import numpy as np
from numpy.typing import NDArray
from numba import jit, prange
from massflow.tools.logger import get_logger

logger = get_logger("peak_align_compute")

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

@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def align_spectra_parallel_jit(
    mz_matrix: NDArray[np.float64],
    intensity_matrix: NDArray[np.float64],
    lengths: NDArray[np.int32],
    reference: NDArray[np.float64],
    tolerance: float,
    code: int,
) -> NDArray[np.float64]:
    """Align a batch of spectra in parallel over the spectrum axis.

    Args:
        mz_matrix: Padded m/z matrix, shape (n_spectra, max_len).
        intensity_matrix: Padded intensity matrix, shape (n_spectra, max_len).
        lengths: Valid lengths for each spectrum.
        reference: Shared reference m/z axis.
        tolerance: Matching tolerance in absolute or relative units.
        code: Distance method code (0/1/2).

    Returns:
        NDArray[np.float64]: Aligned intensities, shape (n_spectra, n_ref).
    """
    n_spectra = mz_matrix.shape[0]
    n_ref = reference.size
    aligned = np.zeros((n_spectra, n_ref), dtype=np.float64)

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        valid_len = int(lengths[s])
        if valid_len <= 0:
            continue

        mz_list = mz_matrix[s, :valid_len]
        intensity = intensity_matrix[s, :valid_len]
        aligned[s] = align_spectrum_jit(
            mz_list=mz_list,
            intensity=intensity,
            reference=reference,
            tolerance=tolerance,
            code=code,
        )

    return aligned

@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def estimate_domain_stats_parallel_jit(
    mz_matrix: NDArray[np.float64],
    lengths: NDArray[np.int32],
    method_code: int,
):
    """Estimate (min_mz, max_mz, resolution) for each spectrum row in parallel.

    Args:
        mz_matrix: Padded m/z matrix, shape (n_spectra, max_len).
        lengths: Valid lengths for each spectrum.
        method_code: 0='x', 1='y', 2='abs'.

    Returns:
        tuple of 1-D arrays:
            mins: shape (n_spectra,)
            maxs: shape (n_spectra,)
            resolutions: shape (n_spectra,)
    """
    n_spectra = mz_matrix.shape[0]
    mins = np.full(n_spectra, np.nan, dtype=np.float64)
    maxs = np.full(n_spectra, np.nan, dtype=np.float64)
    resolutions = np.full(n_spectra, np.nan, dtype=np.float64)
    eps = np.finfo(np.float64).eps

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        valid_len = int(lengths[s])
        if valid_len <= 1:
            continue

        row = mz_matrix[s, :valid_len]

        finite_count = 0
        for i in range(valid_len):
            if np.isfinite(row[i]):
                finite_count += 1

        if finite_count <= 1:
            continue

        vals = np.empty(finite_count, dtype=np.float64)
        k = 0
        for i in range(valid_len):
            v = row[i]
            if np.isfinite(v):
                vals[k] = v
                k += 1

        vals = np.sort(vals)
        mins[s] = vals[0]
        maxs[s] = vals[finite_count - 1]

        best = np.inf
        if method_code == 2:
            # Absolute resolution: min positive diff
            for i in range(1, finite_count):
                d = vals[i] - vals[i - 1]
                if d > eps and d < best:
                    best = d
        else:
            # Relative resolution: min positive |2*(b-a)/(b+a)|
            for i in range(1, finite_count):
                a = vals[i - 1]
                b = vals[i]
                den = b + a
                if den == 0.0:
                    continue
                r = 2.0 * (b - a) / den
                ar = abs(r)
                if ar > eps and ar < best:
                    best = ar

        if best != np.inf:
            resolutions[s] = best

    return mins, maxs, resolutions


@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def search_nearest_matrix_parallel_jit(
    mz_matrix: NDArray[np.float64],
    lengths: NDArray[np.int32],
    targets: NDArray[np.float64],
    tolerance: float,
    code: int,
    nomatch_value: int = -1,
) -> NDArray[np.int64]:
    """Parallel nearest-bin search for a batch of spectra.

    Args:
        mz_matrix: Padded m/z matrix, shape (n_spectra, max_len).
        lengths: Valid lengths for each spectrum.
        targets: Sorted reference/domain axis.
        tolerance: Match tolerance.
        code: Distance method code.
        nomatch_value: Fallback index when no match.

    Returns:
        NDArray[np.int64]: Bin index matrix, shape (n_spectra, max_len).
    """
    n_spectra, max_len = mz_matrix.shape
    out = np.full((n_spectra, max_len), nomatch_value, dtype=np.int64)

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        valid_len = int(lengths[s])
        if valid_len <= 0:
            continue

        out[s, :valid_len] = search_nearest_jit(
            queries=mz_matrix[s, :valid_len],
            targets=targets,
            tolerance=tolerance,
            code=code,
            nomatch_value=nomatch_value,
            force_nearest=False,
        )

    return out
