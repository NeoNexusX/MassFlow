from typing import Literal, Optional, Tuple
import numpy as np
from numpy.typing import NDArray
from massflow.tools.logger import get_logger
from massflow.module.spectrum import Spectrum
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.data_manager.ms_data_manager import MSDataManager
import massflow.preprocess.numba.peak_align_numba as compute

logger = get_logger("peak_alignment")

MIN_RELATIVE_RES = 5e-7  # 0.5 ppm
MIN_ABSOLUTE_RES = 1e-4  # 0.0001 Da

def _normalize_units(units: str) -> Literal["relative", "absolute"]:
    """
    Normalize the unit parameter (internal helper).

    Args:
        units: Input unit string ('ppm'/'relative' or 'Da'/'absolute').

    Returns:
        str: 'relative' or 'absolute'.
    """
    if units in ("ppm", "relative"):
        return "relative"
    return "absolute"

def get_method_code(tol_method: str) -> int:
    """Convert a tolerance/distance method name into an integer code (for JIT functions).

    Args:
        tol_method: Method name - 'x' (relative to x), 'y' (relative to y; common for PPM),
            or 'abs' (absolute difference).

    Returns:
        int: 0='x', 1='y', 2='abs'
    """
    if tol_method == "x":
        return 0
    if tol_method == "y":
        return 1
    return 2

def mad(x: NDArray, constant: float = 1.4826) -> float:
    """Compute Median Absolute Deviation (MAD)."""
    if x.size == 0:
        return np.nan
    median = np.median(x)
    diff = np.abs(x - median)
    return constant * float(np.median(diff))

def estimate_resolution(
    x: NDArray, tolerance: Optional[float] = None, method: Optional[str] = None
) -> float:
    """Estimate the minimum resolution (minimum spacing) of data points."""
    if x.size <= 1:
        return np.nan

    xs = np.sort(x)

    dx = np.diff(xs)

    a, b = xs[:-1], xs[1:]
    rx = 2.0 * (b - a) / (b + a)

    eps = np.finfo(float).eps
    dx = dx[dx > eps]
    rx = rx[np.abs(rx) > eps]

    chosen_method = method
    if chosen_method is None:
        if dx.size and rx.size:
            range_span = xs[-1] - xs[0]
            lhs = mad(dx / range_span) if range_span > 0 else np.nan
            rhs = mad(rx)
            chosen_method = "abs" if lhs < rhs else "x"
        elif dx.size:
            chosen_method = "abs"
        elif rx.size:
            chosen_method = "x"
        else:
            return np.nan

    target_arr = dx if chosen_method == "abs" else rx
    if target_arr.size == 0:
        return np.nan

    res = float(np.min(target_arr))

    if tolerance is not None:
        residuals = np.mod(target_arr, res)
        if not np.all(residuals <= tolerance):
            return np.nan

    return res

def generate_relative_sequence(start: float, end: float, step: float) -> NDArray:
    """Generate a geometric sequence (for PPM/relative scale)."""
    if start <= 0 or end <= 0:
        logger.warning("Start and end values must be positive for relative sequence generation.")
        return np.array([], dtype=np.float64)

    half = step / 2.0
    ratio = (1.0 + half) / (1.0 - half)

    count = int(np.floor(1.0 + (np.log(end) - np.log(start)) / np.log(ratio)))
    indices = np.arange(count, dtype=np.float64)

    return start * np.power(ratio, indices)

def calc_diff(x: NDArray, y: Optional[NDArray] = None, method: str = "y") -> NDArray:
    """Calculate relative or absolute differences between array elements."""
    if y is None:
        if x.size <= 1:
            return np.array([], dtype=np.float64)
        y = x[:-1]
        x = x[1:]

    x = x.astype(np.float64, copy=False)
    y = y.astype(np.float64, copy=False)

    n = max(x.size, y.size)
    if x.size != n:
        x = np.resize(x, n)
    if y.size != n:
        y = np.resize(y, n)

    if method == "abs":
        return x - y

    denominator = x if method == "x" else y
    return (x - y) / denominator

def scalar_diff(val1: float, val2: float, method: str) -> float:
    """Compute the difference between two scalar values (always non-negative)."""
    diff = abs(val1 - val2)
    if method == "abs":
        return diff

    denominator = val1 if method == "x" else val2
    return diff / abs(denominator) if denominator != 0 else float("inf")

def estimate_domain(
    data_manager: MSDataManager,
    binfun: str = "median",
    binratio: float = 2.0,
    units: str = "relative",
    clear_memory: bool = True,
    batch_size: int = 256,
) -> Tuple[NDArray, float]:
    """
    Estimate a shared m/z reference grid for all spectra.

    Algorithm overview:
    1. Iterate all spectra and estimate each spectrum's resolution.
    2. Aggregate a global resolution using a statistic (median/min/max/mean).
    3. Generate the reference m/z axis and its corresponding tolerance.

    Args:
        data_manager: Mass spectrometry data manager.
        binfun: Aggregation function for resolution ('median', 'min', 'max', 'mean').
        binratio: Tolerance multiplier (tolerance = step * binratio).
        units: Unit type ('relative'/'ppm' or 'absolute'/'Da').
        clear_memory: Whether to clear memory after each batch.
        batch_size: Batch size.

    Returns:
        tuple: (domain, tolerance) - reference m/z axis array and tolerance value.
    """
    ref_method = "x" if units == "relative" else "abs"
    stats = []  # Stores (min_mz, max_mz, resolution)

    for batch in data_manager.batch_generator(batch_size=batch_size):
        for spec in batch:
            # Ensure numpy array for type checkers and np.isfinite
            mz_list = np.asarray(spec.mz_list, dtype=np.float64)
            valid_mz = mz_list[np.isfinite(mz_list)]
            res = estimate_resolution(valid_mz, method=ref_method)

            if valid_mz.size > 0 and np.isfinite(res):
                stats.append((float(np.min(valid_mz)), float(np.max(valid_mz)), res))
            else:
                stats.append((np.nan, np.nan, res))

        if clear_memory:
            data_manager.clear_batch_data_memory(batch)

    stats_arr = np.array(stats, dtype=np.float64)

    # Compute global m/z range
    mz_min = float(np.floor(np.nanmin(stats_arr[:, 0])))
    mz_max = float(np.ceil(np.nanmax(stats_arr[:, 1])))
    resolutions = stats_arr[:, 2]

    # Aggregate resolution using the specified statistic
    if binfun == "median":
        step = float(np.nanmedian(resolutions))  # Median: robust
    elif binfun == "min":
        step = float(np.nanmin(resolutions))  # Min: finest
    elif binfun == "max":
        step = float(np.nanmax(resolutions))  # Max: coarsest
    else:
        step = float(np.nanmean(resolutions))  # Mean: trade-off

    logger.info(
        f"Estimated domain stats: min_mz={mz_min}, max_mz={mz_max}, step={step} (method={binfun})"
    )

    # Generate reference axis based on unit type
    if units == "relative":
        # Relative-error mode: use a geometric sequence
        step = round(2.0 * step, 6) * 0.5  # Round to a multiple of 0.5
        step = max(MIN_RELATIVE_RES, step)  # Ensure not below the minimum resolution
        tolerance = (
            round(2.0 * (step * binratio), 6) * 0.5
        )  # tolerance = step * binratio
        domain = generate_relative_sequence(mz_min, mz_max, step)
    else:
        # Absolute-error mode: use an arithmetic sequence
        step = round(step, 4)  # Round to 4 decimal places
        step = max(MIN_ABSOLUTE_RES, step)  # Ensure not below the minimum resolution
        tolerance = round(step * binratio, 4)  # tolerance = step * binratio
        domain = np.arange(mz_min, mz_max + step, step, dtype=np.float64)

    return domain, tolerance

def bin_peaks(
    data_manager: MSDataManager,
    domain: NDArray,
    tolerance: float,
    tol_method: str = "abs",
    batch_size: int = 256,
) -> NDArray:
    """
    Map peaks from multiple spectra onto a unified reference grid.

    Algorithm overview:
    1. For each raw peak, find the nearest matching point on the reference grid.
    2. Resolve one-to-many conflicts: if multiple peaks map to the same grid point,
       keep the one with the smallest distance.
    3. Accumulate matches and compute a weighted average to refine reference peak positions.
    4. Merge adjacent peaks to generate the final reference axis.

    Args:
        data_manager: Mass spectrometry data manager.
        domain: Initial reference grid.
        tolerance: Tolerance threshold.
        tol_method: Tolerance/distance method ('abs', 'x', 'y').
        batch_size: Batch size.

    Returns:
        NDArray: Final reference m/z axis.
    """
    logger.info(
        f"Binning peaks: {len(data_manager.ms)} spectra, domain size {domain.size}, tolerance {tolerance}"
    )

    code = get_method_code(tol_method)

    peaks_acc = np.zeros(domain.size, dtype=np.float64)  # Peak-position accumulator
    counts = np.zeros(domain.size, dtype=np.int64)  # Counter

    for batch in data_manager.batch_generator(batch_size=batch_size):
        for spec in batch:
            raw_peaks = np.array(spec.mz_list, dtype=np.float64)

            if raw_peaks.size == 0:
                continue

            # Find the nearest reference point for each raw peak
            bin_indices = compute.search_nearest_jit(
                raw_peaks,
                domain,
                tolerance,
                code,
                nomatch_value=-1,
                force_nearest=False,
            )

            valid_mask = bin_indices >= 0  # Filter out unmatched peaks
            if not np.any(valid_mask):
                continue

            valid_bins = bin_indices[valid_mask].astype(np.int64)
            valid_peaks = raw_peaks[valid_mask]

            # Compute distance for each match
            dists = np.abs(
                calc_diff(valid_peaks, domain[valid_bins], method=tol_method)
            )

            # Resolve one-to-many conflicts: multiple peaks mapping to the same reference bin
            # Sort by (bin, dist) and keep only the smallest-distance match per bin
            n_valid = valid_bins.size
            sort_arr = np.empty(
                n_valid, dtype=[("bin", "i8"), ("dist", "f8"), ("idx", "i8")]
            )
            sort_arr["bin"] = valid_bins  # Reference-bin index
            sort_arr["dist"] = dists  # Distance
            sort_arr["idx"] = np.arange(n_valid)  # Index into valid_peaks

            # Sort by bin first, then by distance (stable keeps original order for ties)
            sort_arr.sort(order=["bin", "dist"], kind="stable")

            # For each unique bin, take the first one (smallest distance)
            _, unique_indices = np.unique(sort_arr["bin"], return_index=True)

            # Extract best matches
            best_matches = sort_arr[unique_indices]
            matched_bins = best_matches["bin"]
            matched_peaks = valid_peaks[best_matches["idx"]]

            np.add.at(peaks_acc, matched_bins, matched_peaks)
            np.add.at(counts, matched_bins, 1)

        data_manager.clear_batch_data_memory(batch)

    nonzero = counts != 0
    peaks_acc[nonzero] /= counts[nonzero]
    peaks_acc[~nonzero] = np.nan

    # Extract non_NaN values and before merging
    valid_mask = ~np.isnan(peaks_acc)
    valid_peaks = peaks_acc[valid_mask]
    valid_counts = counts[valid_mask]

    reference = merge_peaks(
        mz_list=valid_peaks, counts=valid_counts, tolerance=tolerance, tol_method=tol_method
    )

    return reference

def merge_peaks(
    mz_list: NDArray,
    counts: NDArray,
    tolerance: float,
    tol_method: str = "abs",
) -> NDArray:
    """
    Merge adjacent peaks.

    Algorithm notes:
    Perform an informed merge based on peak frequency (counts), preserving local maxima
    (mode peaks). Starting from each non-NaN peak, expand left/right until one of the
    following conditions is met:
    - Distance exceeds tolerance
    - Encounter NaN
    - Encounter a higher peak (prevents merging across valleys)

    The merged peak position is a weighted average with weights given by counts.

    Args:
        mz_list: m/z value array (may contain NaN).
        counts: Frequency (counts) for each peak.
        tolerance: Merge tolerance.
        tol_method: Tolerance/distance method.

    Returns:
        NDArray: Merged m/z array (NaNs removed).
    """

    # Sortedness check (ignoring NaNs)
    valid_peaks = mz_list[~np.isnan(mz_list)]
    logger.debug(f"Merging peaks: initial count {valid_peaks.size}")
    if valid_peaks.size > 1 and np.any(np.diff(valid_peaks) < 0):
        raise ValueError("peaks must be sorted")

    n = mz_list.size
    # Find the first non-NaN index k
    k = 0
    while k < n and np.isnan(mz_list[k]):
        k += 1

    while k < n:
        i = k
        j = k

        # Expand to the left
        left_of_mode = False  # Whether we already passed the mode peak (local maximum)
        while (i - 1) >= 0:
            if np.isnan(mz_list[i - 1]):
                break  # Encounter NaN
            # Check whether the distance is within tolerance
            dist = scalar_diff(mz_list[i], mz_list[i - 1], tol_method)
            if dist > tolerance:
                break  # Outside tolerance

            # Continue only while descending/flat; stop if we cross a valley (left higher than current)
            if counts[i - 1] < counts[i]:
                left_of_mode = True  # We are on the right side of the peak (downhill)
            if counts[i - 1] > counts[i] and left_of_mode:
                break  # Crossed a valley; stop merging
            i -= 1

        # Expand to the right
        right_of_mode = False  # Whether we already passed the mode peak
        while (j + 1) < n:
            if np.isnan(mz_list[j + 1]):
                break  # Encounter NaN
            dist = scalar_diff(mz_list[j + 1], mz_list[j], tol_method)
            if dist > tolerance:
                break  # Outside tolerance

            if counts[j + 1] < counts[j]:
                right_of_mode = True  # We are on the left side of the peak (downhill)
            if counts[j + 1] > counts[j] and right_of_mode:
                break  # Crossed a valley; stop merging
            j += 1

        # Merge interval [i, j]
        indices = np.arange(i, j + 1)
        w = counts[indices].astype(np.float64)  # Weights
        w_sum = np.sum(w)

        # Compute merged peak position via weighted average
        if w_sum > 0:
            merged_peak = np.sum(w * mz_list[indices]) / w_sum
        else:
            merged_peak = np.mean(
                mz_list[indices]
            )  # Fall back to simple mean if no weights

        merged_count = np.sum(counts[indices])  # Total count after merge

        # Write into the center position p: floor(mean(i, j))
        p = int(np.floor((i + j) / 2.0))

        mz_list[p] = merged_peak
        counts[p] = int(merged_count)

        # Clear other positions
        mask = np.ones(indices.size, dtype=bool)
        mask[p - i] = False  # Keep p

        # Mark merged-away positions as NaN / 0
        cleared_indices = indices[mask]
        mz_list[cleared_indices] = np.nan
        counts[cleared_indices] = 0

        # Move pointer k to the right of j, skipping NaNs
        k = j + 1
        while k < n and np.isnan(mz_list[k]):
            k += 1

    # Drop NaNs and return the final result
    valid = ~np.isnan(mz_list)
    return mz_list[valid]

def compute_reference(
    data_manager: MSDataManager,
    reference: Optional[NDArray] = None,
    binfun="median",
    binratio=2.0,
    tolerance: Optional[float] = None,
    units="ppm",
    clear_memory=True,
    batch_size: int = 256,
) -> Tuple[NDArray, float]:
    """Calculate reference m/z axis and tolerance."""
    logger.info(
        f"Computing reference: units={units}, binfun={binfun}, binratio={binratio}"
    )
    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs"

    # Estimate domain and tolerance
    domain, estimate_tolerance = estimate_domain(
        data_manager=data_manager,
        binfun=binfun,
        binratio=binratio,
        units=norm_units,
        clear_memory=clear_memory,
        batch_size=batch_size,
    )
    logger.info(
        f"Domain estimated: size={domain.size}, estimate_tolerance={estimate_tolerance}"
    )

    tolerance = estimate_tolerance if tolerance is None else (tolerance * 1e-6 if units == "ppm" else tolerance)

    if reference is not None:
        return reference, tolerance

    reference = bin_peaks(
        data_manager=data_manager,
        domain=domain,
        tolerance=tolerance,
        tol_method=tol_method,
        batch_size=batch_size,
    )
    logger.info(f"Reference computed: final size={reference.size}")

    return reference, tolerance

def align_spectrum(
    spectrum: Spectrum,
    reference: NDArray,
    tolerance: float,
    units: str = "ppm",
) -> SpectrumImzML:
    """Align peaks for a single spectrum."""
    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs"
    code = get_method_code(tol_method)

    mz_list = spectrum.mz_list
    intensity = spectrum.intensity

    aligned_intensity = compute.align_spectrum_jit(
        mz_list=mz_list,
        intensity=intensity,
        reference=reference,
        tolerance=tolerance,
        code=code,
    )

    return SpectrumImzML(
        mz_list=reference,
        intensity=aligned_intensity,
        coordinates=spectrum.coordinate,
    )
