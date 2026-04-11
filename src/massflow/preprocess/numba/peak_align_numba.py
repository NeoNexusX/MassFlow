import numpy as np
from numpy.typing import NDArray
from numba import jit, prange


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
def align_spectra_flat_parallel_jit(
    mz_data: NDArray[np.float64],
    intensity_flat: NDArray[np.float64],
    lengths: NDArray[np.int32],
    reference: NDArray[np.float64],
    tolerance: float,
    code: int,
    is_shared_mz: bool,
) -> NDArray[np.float64]:
    """Align a flat batch of spectra in parallel."""
    n_spectra = lengths.size
    n_ref = reference.size
    aligned = np.zeros((n_spectra, n_ref), dtype=np.float64)

    n_total_int = intensity_flat.size
    n_total_mz = mz_data.size

    offsets = np.empty(n_spectra, dtype=np.int64)
    curr = 0
    for i in range(n_spectra):
        offsets[i] = curr
        seg_len = int(lengths[i])
        if seg_len > 0:
            curr += seg_len

    shared_len = n_total_mz

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        valid_len = int(lengths[s])
        if valid_len <= 0:
            continue

        start = int(offsets[s])
        if start >= n_total_int:
            continue

        end = start + valid_len
        if end > n_total_int:
            end = n_total_int

        seg_len = end - start
        if seg_len <= 0:
            continue

        if is_shared_mz:
            mz_len = seg_len if seg_len <= shared_len else shared_len
            if mz_len <= 0:
                continue

            mz_list = mz_data[:mz_len]
            intensity = intensity_flat[start : start + mz_len]
        else:
            if start >= n_total_mz:
                continue

            mz_end = start + seg_len
            if mz_end > n_total_mz:
                mz_end = n_total_mz

            mz_len = mz_end - start
            if mz_len <= 0:
                continue

            mz_list = mz_data[start : start + mz_len]
            intensity = intensity_flat[start : start + mz_len]

        aligned[s] = align_spectrum_jit(
            mz_list=mz_list,
            intensity=intensity,
            reference=reference,
            tolerance=tolerance,
            code=code,
        )

    return aligned

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

@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def estimate_domain_stats_flat_parallel_jit(
    mz_flat: NDArray[np.float64],
    lengths: NDArray[np.int32],
    method_code: int,
):
    """Estimate (min_mz, max_mz, resolution) for each spectrum in flat mode."""
    n_spectra = lengths.size
    mins = np.full(n_spectra, np.nan, dtype=np.float64)
    maxs = np.full(n_spectra, np.nan, dtype=np.float64)
    resolutions = np.full(n_spectra, np.nan, dtype=np.float64)
    eps = np.finfo(np.float64).eps  # pylint: disable=no-member

    n_total = mz_flat.size
    offsets = np.empty(n_spectra, dtype=np.int64)
    curr = 0
    for i in range(n_spectra):
        offsets[i] = curr
        seg_len = int(lengths[i])
        if seg_len > 0:
            curr += seg_len

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        valid_len = int(lengths[s])
        if valid_len <= 1:
            continue

        start = int(offsets[s])
        if start >= n_total:
            continue

        end = start + valid_len
        if end > n_total:
            end = n_total

        seg_len = end - start
        if seg_len <= 1:
            continue

        row = mz_flat[start:end]

        finite_count = 0
        for i in range(seg_len):
            if np.isfinite(row[i]):
                finite_count += 1

        if finite_count <= 1:
            continue

        vals = np.empty(finite_count, dtype=np.float64)
        k = 0
        for i in range(seg_len):
            v = row[i]
            if np.isfinite(v):
                vals[k] = v
                k += 1

        vals = np.sort(vals)
        mins[s] = vals[0]
        maxs[s] = vals[finite_count - 1]

        best = np.inf
        if method_code == 2:
            for i in range(1, finite_count):
                d = vals[i] - vals[i - 1]
                if d > eps and d < best:
                    best = d
        else:
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
def search_nearest_flat_parallel_jit(
    mz_flat: NDArray[np.float64],
    lengths: NDArray[np.int32],
    targets: NDArray[np.float64],
    tolerance: float,
    code: int,
    nomatch_value: int = -1,
) -> NDArray[np.int64]:
    """Parallel nearest-bin search for flat-mode batches."""
    n_total = mz_flat.size
    out = np.full(n_total, nomatch_value, dtype=np.int64)

    n_spectra = lengths.size
    offsets = np.empty(n_spectra, dtype=np.int64)
    curr = 0
    for i in range(n_spectra):
        offsets[i] = curr
        seg_len = int(lengths[i])
        if seg_len > 0:
            curr += seg_len

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        valid_len = int(lengths[s])
        if valid_len <= 0:
            continue

        start = int(offsets[s])
        if start >= n_total:
            continue

        end = start + valid_len
        if end > n_total:
            end = n_total
        if end <= start:
            continue

        out[start:end] = search_nearest_jit(
            queries=mz_flat[start:end],
            targets=targets,
            tolerance=tolerance,
            code=code,
            nomatch_value=nomatch_value,
            force_nearest=False,
        )

    return out

@jit(nopython=True, cache=True)
def _scalar_diff_inf_jit(val1: float, val2: float, method_code: int) -> float:
    """Scalar distance with deterministic inf fallback for zero denominator."""
    diff = abs(val1 - val2)
    if method_code == 2:
        return diff

    denom = val1 if method_code == 0 else val2
    if denom == 0.0:
        return np.inf
    return diff / abs(denom)

@jit(nopython=True, cache=True)
def accumulate_best_matches_for_spectrum_jit(
    raw_peaks: NDArray[np.float64],
    bin_indices: NDArray[np.int64],
    domain: NDArray[np.float64],
    method_code: int,
    peaks_acc: NDArray[np.float64],
    counts: NDArray[np.int64],
) -> None:
    """Resolve one-to-many conflicts per spectrum and update accumulators."""
    n = bin_indices.size

    valid_count = 0
    for i in range(n):
        if bin_indices[i] >= 0:
            valid_count += 1

    if valid_count == 0:
        return

    valid_bins = np.empty(valid_count, dtype=np.int64)
    valid_peaks = np.empty(valid_count, dtype=np.float64)
    valid_dists = np.empty(valid_count, dtype=np.float64)
    valid_order = np.empty(valid_count, dtype=np.int64)

    k = 0
    for i in range(n):
        b = bin_indices[i]
        if b < 0:
            continue

        peak = raw_peaks[i]
        target = domain[b]
        dist = _scalar_diff_inf_jit(peak, target, method_code) # type: ignore
        if np.isnan(dist):
            dist = np.inf

        valid_bins[k] = b
        valid_peaks[k] = peak
        valid_dists[k] = dist
        valid_order[k] = k
        k += 1

    order = np.argsort(valid_bins)

    first = order[0]
    curr_bin = valid_bins[first]
    best_dist = valid_dists[first]
    best_peak = valid_peaks[first]
    best_order = valid_order[first]

    for t in range(1, order.size):
        idx = order[t]
        curr = valid_bins[idx]

        if curr != curr_bin:
            peaks_acc[curr_bin] += best_peak
            counts[curr_bin] += 1

            curr_bin = curr
            best_dist = valid_dists[idx]
            best_peak = valid_peaks[idx]
            best_order = valid_order[idx]
            continue

        dist = valid_dists[idx]
        ord_idx = valid_order[idx]
        if dist < best_dist or (dist == best_dist and ord_idx < best_order):
            best_dist = dist
            best_peak = valid_peaks[idx]
            best_order = ord_idx

    peaks_acc[curr_bin] += best_peak
    counts[curr_bin] += 1

@jit(nopython=True, cache=True)
def accumulate_best_matches_for_flat_jit(
    raw_peaks_flat: NDArray[np.float64],
    bin_indices_flat: NDArray[np.int64],
    lengths: NDArray[np.int32],
    domain: NDArray[np.float64],
    method_code: int,
    peaks_acc: NDArray[np.float64],
    counts: NDArray[np.int64],
) -> None:
    """Resolve one-to-many conflicts for each flat spectrum and update accumulators."""
    n_total_raw = raw_peaks_flat.size
    n_total_idx = bin_indices_flat.size

    offset = 0
    for s in range(lengths.size):
        valid_len = int(lengths[s])
        if valid_len <= 0:
            continue

        start = offset
        end = start + valid_len
        offset = end

        if start >= n_total_raw or start >= n_total_idx:
            continue

        if end > n_total_raw:
            end = n_total_raw
        if end > n_total_idx:
            end = n_total_idx
        if end <= start:
            continue

        accumulate_best_matches_for_spectrum_jit(
            raw_peaks=raw_peaks_flat[start:end],
            bin_indices=bin_indices_flat[start:end],
            domain=domain,
            method_code=method_code,
            peaks_acc=peaks_acc,
            counts=counts,
        )

@jit(nopython=True, cache=True)
def merge_peaks_jit(
    mz_list: NDArray[np.float64],
    counts: NDArray[np.int64],
    tolerance: float,
    code: int,
) -> NDArray[np.float64]:
    """Merge adjacent peaks while preserving local mode barriers."""
    n = mz_list.size

    k = 0
    while k < n and np.isnan(mz_list[k]):
        k += 1

    while k < n:
        i = k
        j = k

        left_of_mode = False
        while (i - 1) >= 0:
            if np.isnan(mz_list[i - 1]):
                break

            dist = _scalar_diff_inf_jit(mz_list[i], mz_list[i - 1], code)
            if dist > tolerance:
                break

            if counts[i - 1] < counts[i]:
                left_of_mode = True
            if counts[i - 1] > counts[i] and left_of_mode:
                break
            i -= 1

        right_of_mode = False
        while (j + 1) < n:
            if np.isnan(mz_list[j + 1]):
                break

            dist = _scalar_diff_inf_jit(mz_list[j + 1], mz_list[j], code)
            if dist > tolerance:
                break

            if counts[j + 1] < counts[j]:
                right_of_mode = True
            if counts[j + 1] > counts[j] and right_of_mode:
                break
            j += 1

        merged_count = 0
        w_sum = 0.0
        weighted_sum = 0.0
        plain_sum = 0.0
        span = j - i + 1

        for idx in range(i, j + 1):
            c = counts[idx]
            merged_count += c
            w = float(c)
            w_sum += w
            weighted_sum += w * mz_list[idx]
            plain_sum += mz_list[idx]

        if w_sum > 0.0:
            merged_peak = weighted_sum / w_sum
        else:
            merged_peak = plain_sum / float(span)

        p = (i + j) // 2
        mz_list[p] = merged_peak
        counts[p] = merged_count

        for idx in range(i, j + 1):
            if idx != p:
                mz_list[idx] = np.nan
                counts[idx] = 0

        k = j + 1
        while k < n and np.isnan(mz_list[k]):
            k += 1

    valid_count = 0
    for idx in range(n):
        if not np.isnan(mz_list[idx]):
            valid_count += 1

    out = np.empty(valid_count, dtype=np.float64)
    out_k = 0
    for idx in range(n):
        if not np.isnan(mz_list[idx]):
            out[out_k] = mz_list[idx]
            out_k += 1

    return out
