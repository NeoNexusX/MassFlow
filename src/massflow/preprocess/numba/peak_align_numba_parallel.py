import numpy as np
from numpy.typing import NDArray
from numba import jit, prange

from massflow.preprocess.numba.peak_align_numba import (
    align_spectrum_jit,
    search_nearest_jit,
)


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
