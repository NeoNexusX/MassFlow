import numpy as np
from typing import Dict, Literal, Optional, Tuple
from numpy.typing import NDArray
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

def scalar_diff(val1: float, val2: float, method: Literal["x", "y", "abs"]) -> float:
    """
    Calculate the absolute difference between two scalar values (result is always positive).
    """
    diff = abs(val1 - val2)
    if method == "abs":
        return diff
    
    denominator = val1 if method == "x" else val2
    return diff / abs(denominator) if denominator != 0 else float('inf')

def search_nearest(
    queries: NDArray[np.float64],
    targets: NDArray[np.float64],
    tolerance: float = 0.0,
    tol_method: Literal["x", "y", "abs"] = "abs",
    nomatch_value: int = -1,
    force_nearest: bool = False,
) -> NDArray[np.int64]:
    """
    Find matching indices for queries in the sorted array targets.
    
    Logic:
    1. Use binary search to locate the insertion point.
    2. Check the left and right neighbors of the insertion point.
    3. Determine the best match based on tolerance and method.

    Args:
        queries: Values to query (any order).
        targets: Sorted target reference array (must be ascending).
        tolerance: Maximum allowed error radius.
        tol_method: Error calculation method ("x"/"y"/"abs").
        nomatch_value: Return value when no match is found (default -1).
        force_nearest: If True, ignore tolerance and force return of the nearest neighbor.

    Returns:
        Index array, same length as queries.
    """
    if targets.size == 0:
        return np.full(queries.shape, nomatch_value, dtype=np.int64)

    if targets.size > 1 and targets[0] > targets[-1]:
            raise ValueError("Targets array must be sorted in ascending order.")

    # 1. Initial binary search position
    # side='left': a[i-1] < v <= a[i]
    insert_indices = np.searchsorted(targets, queries)
    
    result_indices = np.full(queries.shape, nomatch_value, dtype=np.int64)
    n_targets = targets.size

    for i, q_val in enumerate(queries):
        idx = int(insert_indices[i])
        
        # Candidate indices: left neighbor and current position
        idx_left = idx - 1
        idx_right = idx
        
        best_idx = nomatch_value
        min_diff = float('inf')

        # Check left neighbor
        if idx_left >= 0:
            diff_left = scalar_diff(q_val, targets[idx_left], tol_method)
            if force_nearest or diff_left <= tolerance:
                min_diff = diff_left
                best_idx = idx_left

        # Check right neighbor
        if idx_right < n_targets:
            diff_right = scalar_diff(q_val, targets[idx_right], tol_method)
            
            # If right is closer, or no match found on left yet
            is_closer = diff_right < min_diff
            is_within_tol = force_nearest or diff_right <= tolerance
            
            # Update only if right is strictly better than left (Tie-break: prefer left/smaller index)
            if is_within_tol:
                if best_idx == nomatch_value or is_closer:
                    best_idx = idx_right

        result_indices[i] = best_idx

    return result_indices

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
            lhs = mad(dx / range_span) if range_span > 0 else np.inf
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
            return np.inf

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
        clear_memory: bool = False,
) -> Tuple[NDArray[np.float64], float]:
    """
    Estimate the m/z reference grid shared by all spectra.
    """
    ref_method = "x" if units == "relative" else "abs"
    stats = []

    for spec in data_manager.ms:
        mz_list = spec.mz_list
        valid_mz = mz_list[np.isfinite(mz_list)]
        res = estimate_resolution(valid_mz, method=ref_method)
        
        if valid_mz.size > 0 and np.isfinite(res):
            stats.append((float(np.min(valid_mz)), float(np.max(valid_mz)), res))
        else:
            stats.append((np.nan, np.nan, res))

        if clear_memory:
            spec.mz_list = None
            spec.intensity = None
            
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
) -> NDArray[np.float64]:
    """
    Map peaks from multiple spectra to a unified domain.
    """
    logger.info(f"Binning peaks: {len(data_manager.ms)} spectra, domain size {domain.size}")

    # Accumulator initialization
    peaks_acc = np.zeros(domain.size, dtype=np.float64)
    counts = np.zeros(domain.size, dtype=np.int64)

    for spec in data_manager.ms:
        raw_peaks = np.array(spec.mz_list, dtype=np.float64)
        spec.mz_list = None
        spec.intensity = None
        
        # 1. Initial mapping (do not force nearest, only find within tolerance)
        bin_indices = search_nearest(
            raw_peaks, domain, tolerance, tol_method, nomatch_value=-1, force_nearest=False
        )
        
        # 2. Conflict resolution loop
        while True:
            dup_mask = np.zeros(bin_indices.shape, dtype=bool)
            seen_bins: Dict[int, int] = {}
            
            # Mark duplicates
            for idx, bin_idx in enumerate(bin_indices):
                if bin_idx < 0: continue
                b_int = int(bin_idx)
                if b_int in seen_bins:
                    dup_mask[idx] = True
                else:
                    seen_bins[b_int] = idx
            
            if not np.any(dup_mask):
                break # No duplicates, exit loop

            # Handle the first found duplicate Bin
            first_dup_idx = int(np.nonzero(dup_mask)[0][0])
            target_bin = int(bin_indices[first_dup_idx])
            
            # Find all peak indices mapped to this Bin
            conflict_indices = np.where(bin_indices == target_bin)[0]
            
            # Calculate distances of these peaks to the Bin center
            diffs = np.abs(calc_diff(
                raw_peaks[conflict_indices], 
                np.full(conflict_indices.size, domain[target_bin], dtype=np.float64), 
                method=tol_method
            ))

            # Find candidates with distance >= minimum distance, remove the first one
            min_diff = np.min(diffs)
            candidates = conflict_indices[diffs >= min_diff]
            
            # Mark the first candidate in the conflict list as unmatched (-1)
            bin_indices[candidates[0]] = -1
        
        # 3. Accumulate valid matches
        valid_mask = bin_indices >= 0
        matched_bins = bin_indices[valid_mask].astype(np.int64)
        matched_peaks = raw_peaks[valid_mask]
        
        np.add.at(peaks_acc, matched_bins, matched_peaks)
        np.add.at(counts, matched_bins, 1)

    # 4. Calculate mean
    nonzero = counts != 0
    peaks_acc[nonzero] /= counts[nonzero]
    peaks_acc[~nonzero] = np.nan

    # 5. Merge Peaks
    logger.info("Merging peaks...")
    reference = merge_peaks(mz_list=peaks_acc, counts=counts, tolerance=tolerance, tol_method=tol_method)

    # 6. Clean NaN (Drop NA) - Logic from original binpeaks
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
        clear_memory=False) -> Tuple[NDArray[np.float64], float]:
    """Calculate reference m/z axis and tolerance."""
    logger.info(f"Computing reference: units={units}, binfun={binfun}, binratio={binratio}")
    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs"

    domain, tolerance = estimate_domain(data_manager=data_manager, binfun=binfun, binratio=binratio, units=norm_units, clear_memory=clear_memory)
    logger.info(f"Domain estimated: size={domain.size}, tolerance={tolerance}")

    reference = bin_peaks(data_manager=data_manager, domain=domain, tolerance=tolerance, tol_method=tol_method)
    logger.info(f"Reference computed: final size={reference.size}")

    return reference, tolerance

def align_spectrum(
    spectrum: SpectrumImzML,
    reference: NDArray[np.float64],
    tolerance: float,
    units: str = "ppm",
) -> SpectrumImzML:
    """Align peaks for a single spectrum."""
    norm_units = _normalize_units(units)
    tol_method = "x" if norm_units == "relative" else "abs" 

    n_ref = reference.size
    aligned_row = np.zeros(n_ref, dtype=np.float64)

    if spectrum.mz_list.size == 0:
        spectrumed = SpectrumImzML(
                index=spectrum._index,
                mz_list=reference,
                intensity=aligned_row,
                coordinates=spectrum.coordinate
        )
        return spectrumed

    # If reference axis is smaller than (or similar size to) spectrum, or reference axis is unsorted (rare), use Downsample
    if (n_ref <= 2 * spectrum.mz_list.size):
        # Downsample Strategy
        pos = search_nearest(
            queries=reference, 
            targets=spectrum.mz_list, 
            tolerance=tolerance, 
            tol_method=tol_method, 
            nomatch_value=-1, 
        )
        hits = pos >= 0
        if np.any(hits):
            aligned_row[hits] = spectrum.intensity[pos[hits]]
            
    else:
        # 1. Find the starting insertion position of each peak in Ref
        upsample_method = "y" if tol_method == "x" else "x" if tol_method == "y" else "abs"
        
        start_pos = search_nearest(
            queries=spectrum.mz_list,
            targets=reference,
            tolerance=tolerance,
            tol_method=tol_method, 
            nomatch_value=-1,
        )
        
        processed_mask = np.zeros(n_ref, dtype=bool)
        
        for j in range(spectrum.mz_list.size):
            pj = int(start_pos[j])
            if pj < 0: continue
            
            curr_mz = spectrum.mz_list[j]

            def process_bin(i):
                if processed_mask[i]: return True # Stop expansion if blocked
                
                # 1. Check tolerance
                if scalar_diff(float(reference[i]), float(curr_mz), upsample_method) > tolerance: 
                    return True # Stop expansion if out of bounds
                
                # 2. Reconfirm nearest neighbor
                best_match_idx = search_nearest(
                    queries=np.array([reference[i]]), 
                    targets=spectrum.mz_list, 
                    tolerance=tolerance,
                    tol_method=tol_method, 
                    nomatch_value=-1,
                    force_nearest=True # Force find nearest, we already know it's within tolerance
                )[0]
                
                if best_match_idx != -1:
                    aligned_row[i] = spectrum.intensity[best_match_idx]
                    processed_mask[i] = True
                
                return False # Continue expansion

            
            # Expand Right
            for i in range(pj, n_ref):
                if process_bin(i):
                    break
                
            # Expand Left
            for i in range(pj - 1, -1, -1):
                if process_bin(i):
                    break

    spectrumed = SpectrumImzML(
        index=spectrum._index,
        mz_list=reference,
        intensity=aligned_row,
        coordinates=spectrum.coordinate
    )

    spectrum.mz_list = None
    spectrum.intensity = None
    aligned_row = None

    return spectrumed

def align_massdata(
        data_manager: MSDataManagerImzML,
        reference: Optional[NDArray[np.float64]] = None,
        tolerance: Optional[float] = None,
        units: str = "ppm",
        binfun: str = "max",
        binratio: float = 2.0,
    ) -> MSDataManagerImzML:
    """Align all spectra in MSDataManagerImzML."""

    if reference is None or tolerance is None:
        reference, estimate_tolerance = compute_reference(
            data_manager=data_manager,
            binfun=binfun,
            binratio=binratio,
            units=units
        )
    tolerance = estimate_tolerance if tolerance is None else tolerance
    logger.info(f"Starting mass data alignment. Reference size: {reference.size}, Tolerance: {tolerance}")

    aligned_ms = MassSpectrumSet()
    aligned_manager = MSDataManagerImzML(aligned_ms, temp_dir="./temp_align_data")
    aligned_manager.copy_meta(data_manager)
    writer = aligned_manager.writer

    for i, spectrum in enumerate(data_manager.ms):
            aligned_spectrum = align_spectrum(spectrum=spectrum,
                                              reference=reference,
                                              tolerance=tolerance,
                                              units=units)
            aligned_spectrum.swap_out2disk(writer)

            if (i + 1) % 5000 == 0:
                logger.info(f"Aligned {i + 1}/{len(data_manager.ms)} spectra")

    aligned_manager.close_writer()
    aligned_manager.load_full_data_from_file()

    logger.info("Mass data alignment completed.")
    if aligned_manager.ms.meta is not None:
        aligned_manager.ms.meta.continuous = True

    return aligned_manager