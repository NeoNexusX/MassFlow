from typing import Literal, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from massflow.data_manager import MSDataManager
from massflow.tools.funs import prepare_flat_inputs, infer_shared_mz
import massflow.preprocess.numba.peak_align_numba as compute
from massflow.tools.logger import get_logger

logger = get_logger("massflow.peak_alignment")

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

    eps = np.finfo(float).eps # pylint: disable=no-member
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

def _aggregate_resolution_step(resolutions: NDArray[np.float64], binfun: str) -> float:
    """Aggregate per-spectrum resolution into a global step."""
    if binfun == "median":
        return float(np.nanmedian(resolutions))
    if binfun == "min":
        return float(np.nanmin(resolutions))
    if binfun == "max":
        return float(np.nanmax(resolutions))
    return float(np.nanmean(resolutions))

def estimate_domain_parallel(
    data_manager: MSDataManager,
    binfun: str = "median",
    binratio: float = 2.0,
    units: str = "relative",
    batch_size: int = 256,
    max_threads: int = 0,
) -> Tuple[NDArray, float]:
    """Parallel domain estimation based on flat_generator batches."""
    ref_method = "x" if units == "relative" else "abs"
    method_code = get_method_code(ref_method)

    stats_batches: list[NDArray[np.float64]] = []

    for mz_data, _, lengths, _ in data_manager.flat_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    ):
        if mz_data is None:
            continue

        lengths_arr = np.asarray(lengths, dtype=np.int32)
        if lengths_arr.size == 0:
            continue

        mz_arr = np.asarray(mz_data, dtype=np.float64)

        if mz_arr.ndim != 1:
            raise ValueError(f"Unsupported mz_data ndim={mz_arr.ndim}; expected 1.")

        total_points = int(np.sum(lengths_arr, dtype=np.int64))
        max_len = int(np.max(lengths_arr)) if lengths_arr.size > 0 else 0

        # Continuous mode: shared m/z axis (1-D, not flattened per spectrum)
        if mz_arr.size != total_points:
            if mz_arr.size != max_len:
                raise ValueError(
                    "Incompatible flat batch for shared m/z mode: "
                    f"mz_size={mz_arr.size}, max_len={max_len}, total_points={total_points}."
                )

            n_spec = int(lengths_arr.size)
            batch_stats = np.full((n_spec, 3), np.nan, dtype=np.float64)

            valid_mz = mz_arr[np.isfinite(mz_arr)]
            if valid_mz.size > 0:
                res = estimate_resolution(valid_mz, method=ref_method)
                batch_stats[:, 0] = float(np.min(valid_mz))
                batch_stats[:, 1] = float(np.max(valid_mz))
                batch_stats[:, 2] = float(res)

            stats_batches.append(batch_stats)
            continue

        # Processed mode: per-spectrum flattened m/z array (1-D)
        mins, maxs, resolutions = compute.estimate_domain_stats_flat_parallel_jit(
            mz_flat=mz_arr,
            lengths=lengths_arr,
            method_code=method_code,
        )
        stats_batches.append(np.column_stack((mins, maxs, resolutions)))

    if not stats_batches:
        raise ValueError("No flat batches available for domain estimation.")

    stats_arr = np.vstack(stats_batches)

    mz_min = float(np.floor(np.nanmin(stats_arr[:, 0])))
    mz_max = float(np.ceil(np.nanmax(stats_arr[:, 1])))
    resolutions = stats_arr[:, 2]

    step = _aggregate_resolution_step(resolutions, binfun)

    logger.info(
        f"[parallel] Estimated domain stats: min_mz={mz_min}, max_mz={mz_max}, "
        f"step={step} (method={binfun})"
    )

    if units == "relative":
        step = round(2.0 * step, 6) * 0.5
        step = max(MIN_RELATIVE_RES, step)
        tolerance = round(2.0 * (step * binratio), 6) * 0.5
        domain = generate_relative_sequence(mz_min, mz_max, step)
    else:
        step = round(step, 4)
        step = max(MIN_ABSOLUTE_RES, step)
        tolerance = round(step * binratio, 4)
        domain = np.arange(mz_min, mz_max + step, step, dtype=np.float64)

    return domain, tolerance

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

    method_code = get_method_code(tol_method)
    return compute.merge_peaks_jit(
        mz_list=mz_list,
        counts=counts,
        tolerance=float(tolerance),
        code=method_code,
    )

def bin_peaks_parallel(
    data_manager: MSDataManager,
    domain: NDArray,
    tolerance: float,
    tol_method: str = "abs",
    batch_size: int = 256,
    max_threads: int = 0,
) -> NDArray:
    """Parallel binning using flat_generator + parallel nearest search."""
    logger.info(
        f"[parallel] Binning peaks: {len(data_manager.ms)} spectra, "
        f"domain size {domain.size}, tolerance {tolerance}"
    )

    code = get_method_code(tol_method)
    peaks_acc = np.zeros(domain.size, dtype=np.float64)
    counts = np.zeros(domain.size, dtype=np.int64)

    for mz_data, _, lengths, _ in data_manager.flat_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    ):
        if mz_data is None:
            continue

        lengths_arr = np.asarray(lengths, dtype=np.int32)
        if lengths_arr.size == 0:
            continue

        mz_arr = np.asarray(mz_data, dtype=np.float64)

        if mz_arr.ndim != 1:
            raise ValueError(f"Unsupported mz_data ndim={mz_arr.ndim}; expected 1.")

        total_points = int(np.sum(lengths_arr, dtype=np.int64))
        max_len = int(np.max(lengths_arr)) if lengths_arr.size > 0 else 0
        is_shared_mz = mz_arr.size != total_points

        # Continuous mode: shared m/z for all spectra in this batch
        if is_shared_mz:
            if mz_arr.size != max_len:
                raise ValueError(
                    "Incompatible flat batch for shared m/z mode: "
                    f"mz_size={mz_arr.size}, max_len={max_len}, total_points={total_points}."
                )

            raw_peaks = mz_arr
            if raw_peaks.size == 0:
                continue

            bin_indices = compute.search_nearest_jit(
                queries=raw_peaks,
                targets=domain,
                tolerance=tolerance,
                code=code,
                nomatch_value=-1,
                force_nearest=False,
            )

            if np.all(lengths_arr == raw_peaks.size):
                tmp_peaks_acc = np.zeros_like(peaks_acc)
                tmp_counts = np.zeros_like(counts)
                compute.accumulate_best_matches_for_spectrum_jit(
                    raw_peaks=raw_peaks,
                    bin_indices=bin_indices,
                    domain=domain,
                    method_code=code,
                    peaks_acc=tmp_peaks_acc,
                    counts=tmp_counts,
                )

                n_spec = int(lengths_arr.size)
                if n_spec > 0:
                    peaks_acc += tmp_peaks_acc * n_spec
                    counts += tmp_counts * n_spec
            else:
                for valid_len in lengths_arr:
                    peak_len = int(valid_len)
                    if peak_len <= 0:
                        continue

                    compute.accumulate_best_matches_for_spectrum_jit(
                        raw_peaks=raw_peaks[:peak_len],
                        bin_indices=bin_indices[:peak_len],
                        domain=domain,
                        method_code=code,
                        peaks_acc=peaks_acc,
                        counts=counts,
                    )
            continue

        # Processed mode: per-spectrum flattened m/z array
        nearest_idx = compute.search_nearest_flat_parallel_jit(
            mz_flat=mz_arr,
            lengths=lengths_arr,
            targets=domain,
            tolerance=tolerance,
            code=code,
            nomatch_value=-1,
        )

        compute.accumulate_best_matches_for_flat_jit(
            raw_peaks_flat=mz_arr,
            bin_indices_flat=nearest_idx,
            lengths=lengths_arr,
            domain=domain,
            method_code=code,
            peaks_acc=peaks_acc,
            counts=counts,
        )

    nonzero = counts != 0
    peaks_acc[nonzero] /= counts[nonzero]
    peaks_acc[~nonzero] = np.nan

    valid_mask = ~np.isnan(peaks_acc)
    valid_peaks = peaks_acc[valid_mask]
    valid_counts = counts[valid_mask]

    reference = merge_peaks(
        mz_list=valid_peaks,
        counts=valid_counts,
        tolerance=tolerance,
        tol_method=tol_method,
    )

    return reference

def reference_computer(
    data_manager: MSDataManager,
    reference: Optional[NDArray] = None,
    binfun: str = "median",
    binratio: float = 2.0,
    tolerance: Optional[float] = None,
    units: str = "ppm",
    batch_size: int = 256,
    max_threads: int = 2,
) -> Tuple[NDArray, float]:
    """Parallel version of compute_reference based on flat_generator."""
    logger.info(
        f"[parallel] Computing reference: units={units}, binfun={binfun}, binratio={binratio}"
    )

    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs"

    domain, estimate_tolerance = estimate_domain_parallel(
        data_manager=data_manager,
        binfun=binfun,
        binratio=binratio,
        units=norm_units,
        batch_size=batch_size,
        max_threads=max_threads,
    )

    logger.info(
        f"[parallel] Domain estimated: size={domain.size}, estimate_tolerance={estimate_tolerance}"
    )

    tolerance = estimate_tolerance if tolerance is None else (tolerance * 1e-6 if units == "ppm" else tolerance)

    if reference is not None:
        return reference, tolerance

    reference = bin_peaks_parallel(
        data_manager=data_manager,
        domain=domain,
        tolerance=tolerance,
        tol_method=tol_method,
        batch_size=batch_size,
        max_threads=max_threads,
    )

    logger.info(f"[parallel] Reference computed: final size={reference.size}")
    return reference, tolerance

def peak_aligner(
    mz_data: NDArray[np.float64],
    intensity: NDArray[np.float64],
    lengths: NDArray[np.int32],
    reference: NDArray[np.float64],
    tolerance: float,
    units: str = "ppm",
) -> NDArray[np.float64]:
    """Align a flat batch of spectra to a reference m/z axis."""
    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs"
    code = get_method_code(tol_method)

    mz_arr = np.asarray(mz_data, dtype=np.float64)
    intensity_arr = np.asarray(intensity, dtype=np.float64)
    lengths_arr = np.asarray(lengths, dtype=np.int32)
    reference_arr = np.asarray(reference, dtype=np.float64)

    mz_arr , intensity_arr, lengths_arr = prepare_flat_inputs(mz_arr, intensity_arr, lengths_arr)
    is_shared_mz = infer_shared_mz(mz_arr, lengths_arr)

    return compute.align_spectra_flat_parallel_jit(
        mz_data=mz_arr,
        intensity_flat=intensity_arr,
        lengths=lengths_arr,
        reference=reference_arr,
        tolerance=float(tolerance),
        code=code,
        is_shared_mz=is_shared_mz,
    )
