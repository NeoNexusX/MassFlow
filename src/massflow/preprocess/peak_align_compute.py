import numpy as np
from typing import Optional
from numpy.typing import NDArray
from numba import jit
from massflow.logger import get_logger
logger = get_logger("peak_align_compute")

def get_method_code(tol_method: str) -> int:
    """
    Convert a tolerance/distance method name into an integer code (for JIT functions).

    Args:
        tol_method: Method name - 'x' (relative to x), 'y' (relative to y; common for PPM),
            or 'abs' (absolute difference).

    Returns:
        int: 0='x', 1='y', 2='abs'
    """
    if tol_method == "x": return 0
    if tol_method == "y": return 1
    return 2

@jit(nopython=True, fastmath=True, cache=True)
def scalar_diff_jit(val1, val2, method_code):
    """
    JIT-optimized scalar difference calculation (always returns a non-negative value).

    Args:
        val1: First value.
        val2: Second value.
        method_code: 0='x' (relative to val1), 1='y' (relative to val2), 2='abs' (absolute).

    Returns:
        float: The computed difference. Returns NaN when the denominator is 0.
    """
    diff = abs(val1 - val2)
    if method_code == 2: return diff  # Absolute difference
    denom = val1 if method_code == 0 else val2  # Choose denominator
    if denom == 0: return np.nan  # Avoid division by zero
    return diff / abs(denom)  # Relative difference

@jit(nopython=True, fastmath=True, cache=True)
def search_nearest_jit(
    queries: NDArray[np.float64], 
    targets: NDArray[np.float64],
    tolerance: float,
    code: int = 2,
    nomatch_value=-1,
    force_nearest=False):
    """
    Use binary search on a sorted target array to find the nearest match index for each query.

    Args:
        queries: Query values.
        targets: Sorted target array.
        tolerance: Tolerance threshold.
        code: Distance method code (0/1/2).
        nomatch_value: Value returned when there is no match (default: -1).
        force_nearest: If True, always return the nearest neighbor (ignore tolerance).
            If False, the match must be within tolerance.

    Returns:
        NDArray: Target indices for each query; nomatch_value where no match is found.
    """
    n_queries = queries.size
    n_targets = targets.size
    result_indices = np.full(n_queries, nomatch_value, dtype=np.int64)

    if n_targets == 0: return result_indices

    # Find insertion positions via binary search (Numba supports np.searchsorted)
    insert_indices = np.searchsorted(targets, queries)

    for i in range(n_queries):
        idx = insert_indices[i]
        q_val = queries[i]
        
        best_idx = nomatch_value
        min_diff = np.nan

        # Check left neighbor (idx - 1)
        if idx - 1 >= 0:
            diff_left = scalar_diff_jit(q_val, targets[idx-1], code)
            if force_nearest or diff_left <= tolerance:
                min_diff = diff_left
                best_idx = idx - 1

        # Check right neighbor (current insertion point idx)
        if idx < n_targets:
            diff_right = scalar_diff_jit(q_val, targets[idx], code)
            is_closer = diff_right < min_diff  # Is right neighbor closer
            if force_nearest or diff_right <= tolerance:
                if best_idx == nomatch_value or is_closer:
                    best_idx = idx

        result_indices[i] = best_idx
    return result_indices

@jit(nopython=True, fastmath=True, cache=True)
def align_spectrum_jit(mz_list, intensity, reference, tolerance, code):
    """
    JIT-compiled core algorithm for spectrum peak alignment.
    Preserves the original downsampling/upsampling logic and the secondary confirmation step.

    Algorithm:
    1. Downsampling mode (n_ref <= 2*n_peak): for each reference point, find the nearest raw peak.
    2. Upsampling mode (n_ref > 2*n_peak): expand from each raw peak to nearby reference points,
       with a secondary confirmation.

    Args:
        mz_list: Raw spectrum m/z array.
        intensity: Raw spectrum intensity array.
        reference: Reference m/z axis.
        tolerance: Tolerance threshold.
        code: Distance method code.

    Returns:
        NDArray: Intensity array aligned onto the reference axis.
    """
    n_ref = reference.size
    n_peak = mz_list.size
    aligned = np.zeros(n_ref, dtype=np.float64)
    
    if n_peak == 0:
        return aligned

    # Downsampling: fewer reference points; match each reference to nearest raw peak
    if n_ref <= 2 * n_peak:
        pos = search_nearest_jit(reference, mz_list, tolerance, code, -1, False)
        for i in range(n_ref):
            idx = pos[i]
            if idx != -1:
                aligned[i] = intensity[idx]
        return aligned

    # Upsampling: more reference points; expand from each raw peak to fill nearby bins
    # Switch tolerance basis (x<->y) to maintain numeric symmetry
    upsample_code = 1 if code == 0 else (0 if code == 1 else 2)

    # Find a starting reference position for each raw peak
    start_pos = search_nearest_jit(mz_list, reference, tolerance, code, -1, False)
    processed_mask = np.zeros(n_ref, dtype=np.bool_)  # Mark reference bins that are already processed

    for j in range(n_peak):
        pj = start_pos[j]  # Starting reference position for the j-th raw peak
        if pj == -1: continue  # No match
        
        curr_mz = mz_list[j]
        pj_int = int(pj)

        # Expand to the right
        for i in range(pj_int, n_ref):

            if processed_mask[i] or scalar_diff_jit(reference[i], curr_mz, upsample_code) > tolerance: break

            q_arr = np.array([reference[i]], dtype=np.float64)
            best_match_arr = search_nearest_jit(q_arr, mz_list, tolerance, code, -1, True)
            best_match_idx = best_match_arr[0]

            if best_match_idx != -1:
                aligned[i] = intensity[best_match_idx]
                processed_mask[i] = True

        # Expand to the left
        for i in range(pj_int - 1, -1, -1):
            if processed_mask[i] or scalar_diff_jit(reference[i], curr_mz, upsample_code) > tolerance: break

            q_arr = np.array([reference[i]], dtype=np.float64)
            best_match_arr = search_nearest_jit(q_arr, mz_list, tolerance, code, -1, True)
            best_match_idx = best_match_arr[0]

            if best_match_idx != -1:
                aligned[i] = intensity[best_match_idx]
                processed_mask[i] = True

    return aligned

def calc_diff(
    x: NDArray,
    y: Optional[NDArray] = None,
    method: str = "y"
) -> NDArray:
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

def scalar_diff(val1: float, val2: float, method: str) -> float:
    """
    Compute the difference between two scalar values (always non-negative).

    Args:
        val1: First value.
        val2: Second value.
        method: 'x' (relative to val1), 'y' (relative to val2), or 'abs' (absolute).

    Returns:
        float: Computed difference; returns inf when the denominator is 0.
    """
    diff = abs(val1 - val2)
    if method == "abs":
        return diff
    
    denominator = val1 if method == "x" else val2
    return diff / abs(denominator) if denominator != 0 else float('inf')

def generate_relative_sequence(start: float, end: float, step: float) -> NDArray:
    """
    Generate a geometric sequence (for PPM/relative scale).

    Under relative error (e.g., PPM), an equally spaced relative step corresponds to a geometric
    progression.
    Formula: x[n] = start * ratio^n, where ratio = (1 + step/2) / (1 - step/2)

    Args:
        start: Start value (must be > 0).
        end: End value (must be > 0).
        step: Relative step (e.g., 5e-6 means 5 ppm).

    Returns:
        NDArray: Geometric sequence array.
    """
    if start <= 0 or end <= 0:
        logger.warning("Start and end values must be positive for relative sequence generation.")
        return np.array([], dtype=np.float64)
        
    half = step / 2.0
    ratio = (1.0 + half) / (1.0 - half)  # Common ratio
    
    # Sequence length: count = 1 + log(end/start) / log(ratio)
    count = int(np.floor(1.0 + (np.log(end) - np.log(start)) / np.log(ratio)))
    indices = np.arange(count, dtype=np.float64)
    
    return start * np.power(ratio, indices)  # x[n] = start * ratio^n

def mad(x: NDArray, constant: float = 1.4826) -> float:
    """
    Compute Median Absolute Deviation (MAD).

    MAD is a robust dispersion metric that is insensitive to outliers.
    Formula: MAD = constant * median(|x - median(x)|)

    Args:
        x: Input array.
        constant: Normalization constant (default 1.4826 makes MAD comparable to std for normal data).

    Returns:
        float: MAD value.
    """
    if x.size == 0:
        return np.nan
    median = np.median(x)
    diff = np.abs(x - median)
    return constant * float(np.median(diff))

def estimate_resolution(
    x: NDArray, 
    tolerance: Optional[float] = None, 
    method: Optional[str] = None
) -> float:
    """
    Estimate the minimum resolution (minimum spacing) of data points.
    Automatically chooses whether absolute spacing (Da) or relative spacing better describes the
    distribution.

    Algorithm:
    1. Compute absolute adjacent differences dx and relative adjacent differences rx.
    2. Compare compactness of the two distributions using MAD.
    3. Return the minimum of the chosen distribution as the resolution estimate.

    Args:
        x: Input data points.
        tolerance: Optional consistency-check tolerance.
        method: Force a method ('abs' or 'x'); if None, choose automatically.

    Returns:
        float: Estimated resolution value.
    """
    if x.size <= 1:
        return np.nan

    xs = np.sort(x)
    
    # Absolute differences: dx[i] = xs[i+1] - xs[i]
    dx = np.diff(xs)
    
    # Relative differences: rx[i] = 2*(b-a)/(b+a), approx (b-a)/mean(a,b)
    a, b = xs[:-1], xs[1:]
    rx = 2.0 * (b - a) / (b + a)
    
    # Filter invalid values below machine precision
    eps = np.finfo(float).eps
    dx = dx[dx > eps]
    rx = rx[np.abs(rx) > eps]

    # Auto-selection: compare MAD and choose the tighter distribution
    chosen_method = method
    if chosen_method is None:
        if dx.size and rx.size:
            # Normalize dx for a fair comparison with rx
            range_span = xs[-1] - xs[0]
            lhs = mad(dx / range_span) if range_span > 0 else np.nan  # MAD of normalized absolute diffs
            rhs = mad(rx)  # MAD of relative diffs
            chosen_method = "abs" if lhs < rhs else "x"  # Smaller MAD indicates a tighter distribution
        elif dx.size:
            chosen_method = "abs"
        elif rx.size:
            chosen_method = "x"
        else:
            return np.nan

    # Estimated resolution: minimum spacing in the chosen distribution
    target_arr = dx if chosen_method == "abs" else rx
    if target_arr.size == 0:
            return np.nan

    res = float(np.min(target_arr))

    # Optional: consistency check (verify all intervals are multiples of res)
    if tolerance is not None:
        residuals = np.mod(target_arr, res)
        if not np.all(residuals <= tolerance):
            return np.nan  # Irregular data; cannot estimate

    return res