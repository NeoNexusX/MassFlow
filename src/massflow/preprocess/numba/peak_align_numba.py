import numpy as np
from numpy.typing import NDArray
from numba import jit
from massflow.tools.logger import get_logger

logger = get_logger("massflow.peak_align_numba")

@jit(nopython=True, fastmath=True, cache=True)
def scalar_diff_jit(val1, val2, method_code):
    """JIT-optimized scalar difference calculation (always returns a non-negative value).

    Args:
        val1: First value.
        val2: Second value.
        method_code: 0='x' (relative to val1), 1='y' (relative to val2), 2='abs' (absolute).

    Returns:
        float: The computed difference. Returns NaN when the denominator is 0.
    """
    diff = abs(val1 - val2)
    if method_code == 2:
        return diff  # Absolute difference
    denom = val1 if method_code == 0 else val2  # Choose denominator
    if denom == 0:
        return np.nan  # Avoid division by zero
    return diff / abs(denom)  # Relative difference

@jit(nopython=True, fastmath=True, cache=True)
def search_nearest_jit(
    queries: NDArray[np.float64],
    targets: NDArray[np.float64],
    tolerance: float,
    code: int = 2,
    nomatch_value=-1,
    force_nearest=False,
):
    """Use binary search on a sorted target array to find the nearest match index for each query.

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

    if n_targets == 0:
        return result_indices

    insert_indices = np.searchsorted(targets, queries)

    for i in range(n_queries):
        idx = insert_indices[i]
        q_val = queries[i]

        best_idx = nomatch_value
        min_diff = np.nan

        if idx - 1 >= 0:
            diff_left = scalar_diff_jit(q_val, targets[idx - 1], code)
            if force_nearest or diff_left <= tolerance:
                min_diff = diff_left
                best_idx = idx - 1

        if idx < n_targets:
            diff_right = scalar_diff_jit(q_val, targets[idx], code)
            is_closer = diff_right < min_diff
            if force_nearest or diff_right <= tolerance:
                if best_idx == nomatch_value or is_closer:
                    best_idx = idx

        result_indices[i] = best_idx
    return result_indices

@jit(nopython=True, fastmath=True, cache=True)
def align_spectrum_jit(mz_list, intensity, reference, tolerance, code):
    """JIT-compiled core algorithm for spectrum peak alignment.

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
        if pj == -1:
            continue

        curr_mz = mz_list[j]
        pj_int = int(pj)

        for i in range(pj_int, n_ref):
            if processed_mask[i] or scalar_diff_jit(reference[i], curr_mz, upsample_code) > tolerance:
                break

            q_arr = np.array([reference[i]], dtype=np.float64)
            best_match_arr = search_nearest_jit(q_arr, mz_list, tolerance, code, -1, True)
            best_match_idx = best_match_arr[0]

            if best_match_idx != -1:
                aligned[i] = intensity[best_match_idx]
                processed_mask[i] = True

        for i in range(pj_int - 1, -1, -1):
            if processed_mask[i] or scalar_diff_jit(reference[i], curr_mz, upsample_code) > tolerance:
                break

            q_arr = np.array([reference[i]], dtype=np.float64)
            best_match_arr = search_nearest_jit(q_arr, mz_list, tolerance, code, -1, True)
            best_match_idx = best_match_arr[0]

            if best_match_idx != -1:
                aligned[i] = intensity[best_match_idx]
                processed_mask[i] = True

    return aligned
