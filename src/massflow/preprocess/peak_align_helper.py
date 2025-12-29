import numpy as np
from typing import Literal, Optional, Tuple
from numpy.typing import NDArray
from numba import jit
from massflow.logger import get_logger
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML

logger = get_logger("peak_alignment")

def calc_diff(
    x: NDArray[np.float64], 
    y: Optional[NDArray[np.float64]] = None, 
    method: Literal["x", "y", "abs"] = "y"
) -> NDArray[np.float64]:
    """
    Calculate relative or absolute differences between array elements.

    Args:
        x: Target array.
        y: Reference array. If None, calculates the difference between adjacent elements of x (x[1:] - x[:-1]).
        method: Basis for difference calculation.
            - "x": (x - y) / x
            - "y": (x - y) / y  (Default, commonly used for PPM calculation)
            - "abs": x - y

    Returns:
        Array of calculated differences.
    """
    if y is None:
        if x.size <= 1:
            return np.array([], dtype=np.float64)
        y = x[:-1]
        x = x[1:]

    x = x.astype(np.float64, copy=False)
    y = y.astype(np.float64, copy=False)

    n = max(x.size, y.size)
    if x.size != n: x = np.resize(x, n)
    if y.size != n: y = np.resize(y, n)

    if method == "abs":
        return x - y
    
    denominator = x if method == "x" else y
    return (x - y) / denominator

@jit(nopython=True, fastmath=True, cache=True)
def scalar_diff_jit(val1, val2, method_code):
    """JIT optimized scalar difference calculation."""
    # method_code: 0='x', 1='y', 2='abs'
    diff = abs(val1 - val2)
    if method_code == 2: return diff
    denom = val1 if method_code == 0 else val2
    if denom == 0: return np.nan
    return diff / abs(denom)

def scalar_diff(val1: float, val2: float, method: Literal["x", "y", "abs"]) -> float:
    """
    Calculate the absolute difference between two scalar values (result is always positive).
    """
    diff = abs(val1 - val2)
    if method == "abs":
        return diff
    
    denominator = val1 if method == "x" else val2
    return diff / abs(denominator) if denominator != 0 else float('inf')

@jit(nopython=True, fastmath=True, cache=True)
def search_nearest_jit(
    queries: NDArray[np.float64], 
    targets: NDArray[np.float64],
    tolerance: float,
    code: int = 2,
    nomatch_value=-1,
    force_nearest=False):

    n_queries = queries.size
    n_targets = targets.size
    result_indices = np.full(n_queries, nomatch_value, dtype=np.int64)

    if n_targets == 0: return result_indices

    # Binary search (Numba supports np.searchsorted)
    insert_indices = np.searchsorted(targets, queries)

    for i in range(n_queries):
        idx = insert_indices[i]
        q_val = queries[i]
        
        best_idx = nomatch_value
        min_diff = np.nan

        # Check left neighbor
        if idx - 1 >= 0:
            diff_left = scalar_diff_jit(q_val, targets[idx-1], code)
            if force_nearest or diff_left <= tolerance:
                min_diff = diff_left
                best_idx = idx - 1

        # Check right neighbor (current insertion point)
        if idx < n_targets:
            diff_right = scalar_diff_jit(q_val, targets[idx], code)
            is_closer = diff_right < min_diff
            if force_nearest or diff_right <= tolerance:
                if best_idx == nomatch_value or is_closer:
                    best_idx = idx

        result_indices[i] = best_idx
    return result_indices

def search_nearest(queries: NDArray[np.float64],
                   targets: NDArray[np.float64],
                   tolerance: float,
                   tol_method: str,
                   nomatch_value=-1,
                   force_nearest=False):
    code = 0 if tol_method == "x" else (1 if tol_method == "y" else 2)

    return search_nearest_jit(
        queries=queries,
        targets=targets,
        tolerance=tolerance,
        code=code,
        nomatch_value=nomatch_value,
        force_nearest=force_nearest
    )

def generate_relative_sequence(start: float, end: float, step: float) -> NDArray[np.float64]:
    """
    Generate a geometric sequence (for PPM scale).
    """
    if start <= 0 or end <= 0:
        return np.array([], dtype=np.float64)
        
    half = step / 2.0
    ratio = (1.0 + half) / (1.0 - half)
    
    # Calculate sequence length
    count = int(np.floor(1.0 + (np.log(end) - np.log(start)) / np.log(ratio)))
    indices = np.arange(count, dtype=np.float64)
    
    return start * np.power(ratio, indices)

def mad(x: NDArray[np.float64], constant: float = 1.4826) -> float:
    """Calculate Median Absolute Deviation (MAD)."""
    if x.size == 0:
        return np.nan
    median = np.median(x)
    diff = np.abs(x - median)
    return constant * float(np.median(diff))

def estimate_resolution(
    x: NDArray[np.float64], 
    tolerance: Optional[float] = None, 
    method: Optional[Literal["x", "abs"]] = None
) -> float:
    """
    Estimate the minimum resolution (minimum interval) of data points.
    Automatically determines whether absolute error (Da) or relative error (PPM) better describes the data distribution.
    """
    if x.size <= 1:
        return np.nan

    xs = np.sort(x)
    
    # Calculate absolute differences
    dx = np.diff(xs)
    
    # Calculate relative differences (2 * (b-a)/(b+a))
    a, b = xs[:-1], xs[1:]
    rx = 2.0 * (b - a) / (b + a)
    
    # Filter invalid values below machine precision
    eps = np.finfo(float).eps
    dx = dx[dx > eps]
    rx = rx[np.abs(rx) > eps]

    # Automatic selection strategy: Compare MAD to see which distribution is more compact
    chosen_method = method
    if chosen_method is None:
        if dx.size and rx.size:
            # Normalize dx for comparison
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

    # Return estimated value
    target_arr = dx if chosen_method == "abs" else rx
    if target_arr.size == 0:
            return np.nan

    res = float(np.min(target_arr))

    # Optional: Consistency check
    if tolerance is not None:
        residuals = np.mod(target_arr, res)
        if not np.all(residuals <= tolerance):
            return np.nan

    return res

def _normalize_units(units: str) -> Literal["relative", "absolute"]:
    """Normalize unit parameters (Internal helper)"""
    if units in ("ppm", "relative"):
        return "relative"
    return "absolute"

def estimate_domain(
        data_manager: MSDataManagerImzML,
        binfun: str = "median",
        binratio: float = 2.0,
        units: str = "relative",
        clear_memory: bool = True,
        batch_size: int = 256,
) -> Tuple[NDArray[np.float64], float]:
    """
    Estimate the m/z reference grid shared by all spectra.
    """
    ref_method = "x" if units == "relative" else "abs"
    stats = []

    for batch in data_manager.get_batch_generator(batch_size=batch_size):
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

    mz_min = float(np.floor(np.nanmin(stats_arr[:, 0])))
    mz_max = float(np.ceil(np.nanmax(stats_arr[:, 1])))
    resolutions = stats_arr[:, 2]

    if binfun == "median":
        step = float(np.nanmedian(resolutions))
    elif binfun == "min":
        step = float(np.nanmin(resolutions))
    elif binfun == "max":
        step = float(np.nanmax(resolutions))
    else:
        step = float(np.nanmean(resolutions))

    logger.info(f"Estimated domain stats: min_mz={mz_min}, max_mz={mz_max}, step={step} (method={binfun})")

    if units == "relative":
        tolerance = round(2.0 * (step * binratio), 6) * 0.5
        step = round(2.0 * step, 6) * 0.5
        domain = generate_relative_sequence(mz_min, mz_max, step)
    else:
        tolerance = round(step * binratio, 4)
        step = round(step, 4)
        domain = np.arange(mz_min, mz_max + step, step, dtype=np.float64)

    return domain, tolerance

def bin_peaks(
        data_manager: MSDataManagerImzML,
        domain: NDArray[np.float64],
        tolerance: float,
        tol_method: Literal["abs", "x", "y"] = "abs",
        batch_size: int = 256,
) -> NDArray[np.float64]:
    """
    Map peaks from multiple spectra to a unified domain.
    """
    logger.info(f"Binning peaks: {len(data_manager.ms)} spectra, domain size {domain.size}")

    peaks_acc = np.zeros(domain.size, dtype=np.float64)
    counts = np.zeros(domain.size, dtype=np.int64)

    for batch in data_manager.get_batch_generator(batch_size=batch_size):
        for spec in batch:
            raw_peaks = np.array(spec.mz_list, dtype=np.float64)

            if raw_peaks.size == 0:
                continue

            bin_indices = search_nearest(
                raw_peaks, 
                domain, 
                tolerance, 
                tol_method, 
                nomatch_value=-1, 
                force_nearest=False)
            
            valid_mask = bin_indices >= 0
            if not np.any(valid_mask):
                continue
                
            valid_bins = bin_indices[valid_mask].astype(np.int64)
            valid_peaks = raw_peaks[valid_mask]
            valid_indices = np.arange(raw_peaks.size)[valid_mask]

            dists = np.abs(calc_diff(valid_peaks, domain[valid_bins], method=tol_method))

            n_valid = valid_bins.size
            sort_arr = np.empty(n_valid, dtype=[('bin', 'i8'), ('dist', 'f8'), ('idx', 'i8')])
            sort_arr['bin'] = valid_bins
            sort_arr['dist'] = dists
            sort_arr['idx'] = np.arange(n_valid)
            

            sort_arr.sort(order=['bin', 'dist'], kind='stable')

            _, unique_indices = np.unique(sort_arr['bin'], return_index=True)
            

            best_matches = sort_arr[unique_indices]
            matched_bins = best_matches['bin']
            matched_peaks = valid_peaks[best_matches['idx']]

            np.add.at(peaks_acc, matched_bins, matched_peaks)
            np.add.at(counts, matched_bins, 1)

        data_manager.clear_batch_data_memory(batch)

    nonzero = counts != 0
    peaks_acc[nonzero] /= counts[nonzero]
    peaks_acc[~nonzero] = np.nan

    logger.info("Merging peaks...")
    reference = merge_peaks(mz_list=peaks_acc, counts=counts, tolerance=tolerance, tol_method=tol_method)

    valid = ~np.isnan(reference)
    reference = reference[valid]
    return reference

def merge_peaks(
    mz_list: NDArray[np.float64],
    counts: NDArray[np.int64],
    tolerance: float,
    tol_method: Literal["abs", "x", "y"] = "abs",
) -> NDArray[np.float64]:
    """
    Merge adjacent peaks.
    """
    
    # Sort check (ignore NaN)
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
        
        # Left Expansion
        left_of_mode = False
        while (i - 1) >= 0:
            if np.isnan(mz_list[i - 1]):
                break
            # Check if distance is within tolerance
            dist = scalar_diff(mz_list[i], mz_list[i - 1], tol_method)
            if dist > tolerance:
                break
            
            # Continue only if monotonically decreasing or flat; stop if it becomes a valley (left is higher than current)
            if counts[i - 1] < counts[i]:
                left_of_mode = True
            if counts[i - 1] > counts[i] and left_of_mode:
                break
            i -= 1

        # Right Expansion
        right_of_mode = False
        while (j + 1) < n:
            if np.isnan(mz_list[j + 1]):
                break
            dist = scalar_diff(mz_list[j + 1], mz_list[j], tol_method)
            if dist > tolerance:
                break
            
            if counts[j + 1] < counts[j]:
                right_of_mode = True
            if counts[j + 1] > counts[j] and right_of_mode:
                break
            j += 1
        
        # Merge interval [i, j]
        indices = np.arange(i, j + 1)
        w = counts[indices].astype(np.float64)
        w_sum = np.sum(w)
        
        # Weighted average
        if w_sum > 0:
            merged_peak = np.sum(w * mz_list[indices]) / w_sum
        else:
            merged_peak = np.mean(mz_list[indices])
        
        merged_count = np.sum(counts[indices])
        
        # Write to center position p
        p = int(np.floor((i + j) / 2.0)) # Original logic: floor(mean(ij))
        
        mz_list[p] = merged_peak
        counts[p] = int(merged_count)
        
        # Clear other positions
        mask = np.ones(indices.size, dtype=bool)
        mask[p - i] = False # Keep p
        
        # Set merged positions to NaN / 0
        cleared_indices = indices[mask]
        mz_list[cleared_indices] = np.nan
        counts[cleared_indices] = 0
        
        # Move pointer k to the right of j, and skip NaN
        k = j + 1
        while k < n and np.isnan(mz_list[k]):
            k += 1

    # Drop NA
    valid = ~np.isnan(mz_list)
    return mz_list[valid]

def compute_reference(
        data_manager: MSDataManagerImzML,
        binfun="median",
        binratio=2.0,
        units="ppm",
        clear_memory=True,
        batch_size:int=256,
        ) -> Tuple[NDArray[np.float64], float]:
    """Calculate reference m/z axis and tolerance."""
    logger.info(f"Computing reference: units={units}, binfun={binfun}, binratio={binratio}")
    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs"

    domain, tolerance = estimate_domain(
        data_manager=data_manager, 
        binfun=binfun, binratio=binratio, 
        units=norm_units, 
        clear_memory=clear_memory,
        batch_size=batch_size,
        )
    logger.info(f"Domain estimated: size={domain.size}, tolerance={tolerance}")

    reference = bin_peaks(
        data_manager=data_manager, 
        domain=domain, 
        tolerance=tolerance, 
        tol_method=tol_method,
        batch_size=batch_size,
        )
    logger.info(f"Reference computed: final size={reference.size}")

    return reference, tolerance

@jit(nopython=True, fastmath=True, cache=True)
def align_spectrum_jit(mz_list, intensity, reference, tolerance, code):
    """
    JIT-compiled core for spectrum alignment.
    Strictly preserves the original Downsampling/Upsampling logic and Secondary Confirmation.
    """
    n_ref = reference.size
    n_peak = mz_list.size
    aligned = np.zeros(n_ref, dtype=np.float64)
    
    if n_peak == 0:
        return aligned

    if n_ref <= 2 * n_peak:
        pos = search_nearest_jit(reference, mz_list, tolerance, code, -1, False)
        for i in range(n_ref):
            idx = pos[i]
            if idx != -1:
                aligned[i] = intensity[idx]
        return aligned

    upsample_code = 1 if code == 0 else (0 if code == 1 else 2)

    start_pos = search_nearest_jit(mz_list, reference, tolerance, code, -1, False)
    processed_mask = np.zeros(n_ref, dtype=np.bool_)

    for j in range(n_peak):
        pj = start_pos[j]
        if pj == -1: continue
        
        curr_mz = mz_list[j]
        pj_int = int(pj)

        def check_bin(i):
            if processed_mask[i]: return 1 # Blocked
            

            if scalar_diff_jit(reference[i], curr_mz, upsample_code) > tolerance:
                return 1 
            
            q_arr = np.array([reference[i]], dtype=np.float64)
            best_match_arr = search_nearest_jit(q_arr, mz_list, tolerance, code, -1, True)
            best_match_idx = best_match_arr[0]

            if best_match_idx != -1:
                aligned[i] = intensity[best_match_idx]
                processed_mask[i] = True
            
            return 0

        # Expand Right
        for i in range(pj_int, n_ref):
            if check_bin(i) == 1: break

        # Expand Left
        for i in range(pj_int - 1, -1, -1):
            if check_bin(i) == 1: break

    return aligned


def align_spectrum(
    spectrum: SpectrumImzML,
    reference: NDArray[np.float64],
    tolerance: float,
    units: str = "ppm",
) -> SpectrumImzML:
    """Align peaks for a single spectrum."""
    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs" 
    code = 0 if tol_method == "x" else (1 if tol_method == "y" else 2)

    mz_list = spectrum.mz_list
    intensity = spectrum.intensity

    aligned_intensity = align_spectrum_jit(
        mz_list=mz_list,
        intensity=intensity,
        reference=reference,
        tolerance=tolerance,
        code=code
    )

    spectrumed = SpectrumImzML(
        index=spectrum._index,
        mz_list=reference,
        intensity=aligned_intensity,
        coordinates=spectrum.coordinate
    )

    spectrum.clear_data()
    return spectrumed

def align_massdata(
        data_manager: MSDataManagerImzML,
        reference: Optional[NDArray[np.float64]] = None,
        tolerance: Optional[float] = None,
        units: str = "ppm",
        binfun: str = "max",
        binratio: float = 2.0,
        batch_size: int = 256,
        clear_memory: bool = True
    ) -> MSDataManagerImzML:
    """Align all spectra in MSDataManagerImzML."""

    estimate_tolerance: float = tolerance if tolerance is not None else np.nan

    if reference is None or tolerance is None:
        reference, estimate_tolerance = compute_reference(
            data_manager=data_manager,
            binfun=binfun,
            binratio=binratio,
            units=units,
            batch_size=batch_size,
            clear_memory=clear_memory
        )
    tolerance = estimate_tolerance if tolerance is None else tolerance
    logger.info(f"Starting mass data alignment. Reference size: {reference.size}, Tolerance: {tolerance}")

    aligned_ms = MassSpectrumSet()
    aligned_manager = MSDataManagerImzML(aligned_ms, temp_dir="./temp_align_data")
    aligned_manager.copy_meta(data_manager)
    writer = aligned_manager.writer

    total_spectra = len(data_manager.ms)
    processed = 0

    for batch_idx, batch in enumerate(data_manager.get_batch_generator(batch_size=batch_size), start=1):
        for spectrum in batch:
            aligned_spectrum = align_spectrum(
                spectrum=spectrum,
                reference=reference,
                tolerance=tolerance,
                units=units,
            )
            aligned_spectrum.swap_out2disk(writer)
        processed += len(batch)
        logger.info(f"Aligned batch {batch_idx}: {len(batch)} spectra (cumulative {processed}/{total_spectra})")
        data_manager.clear_batch_data_memory(batch)

    aligned_manager.close_writer()
    aligned_manager.load_full_data_from_file()

    logger.info("Mass data alignment completed.")
    if aligned_manager.ms.meta is not None:
        aligned_manager.ms.meta.continuous = True

    return aligned_manager