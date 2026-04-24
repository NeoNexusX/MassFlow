from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from numba import jit, prange

from massflow.tools.funs import lengths_to_offsets
from massflow.tools.logger import get_logger

logger = get_logger("massflow.peak_pick_numba_parallel")


@jit(nopython=True, cache=True)
def _indices_from_bool_mask_jit(mask: NDArray[np.bool_]) -> NDArray[np.int64]:
    """Convert a boolean mask to a compact index array to avoid Python-level boolean indexing overhead."""
    count = 0
    for i in range(mask.size):
        if mask[i]:
            count += 1

    out = np.empty(count, dtype=np.int64)
    k = 0
    for i in range(mask.size):
        if mask[i]:
            out[k] = i
            k += 1
    return out


@jit(nopython=True, cache=True)
def _local_maxima_mask_jit(x: NDArray[np.float64], width: int) -> NDArray[np.bool_]:
    """Detect local maxima."""
    n = x.size
    out = np.zeros(n, dtype=np.bool_)
    if n == 0:
        return out

    width_val = width if width > 0 else 1
    radius = abs(width_val // 2)

    for i in range(n):
        if i < radius or i > n - radius:
            continue

        xi = x[i]
        if not np.isfinite(xi):
            continue

        left = i - radius
        if left < 0:
            left = 0
        right = i + radius
        if right >= n:
            right = n - 1

        is_peak = False
        for j in range(left, right + 1):
            xj = x[j]
            if not np.isfinite(xj):
                continue

            if xi > xj:
                is_peak = True
            if j < i and xj >= xi:
                is_peak = False
                break
            if j > i and xj > xi:
                is_peak = False
                break

        if is_peak:
            out[i] = True

    return out


@jit(nopython=True, cache=True)
def _peak_lbound_jit(x: NDArray[np.float64], peak: int) -> int:
    """Search for the left boundary starting from the peak apex."""
    n = x.size
    lbound = peak
    is_left_of_peak = False

    i = peak - 1
    while i >= 0:
        if x[i] < x[lbound]:
            lbound = i
            is_left_of_peak = True
        elif x[i] > x[lbound] and is_left_of_peak:
            cand = lbound
            lwindow = cand - 2
            if lwindow < 0:
                lwindow = 0

            i -= 1
            while i >= lwindow:
                if x[i] < x[cand]:
                    lbound = i
                    break
                i -= 1

            if cand == lbound:
                break

        i -= 1

    if lbound < 0:
        return 0
    if lbound >= n:
        return n - 1
    return lbound


@jit(nopython=True, cache=True)
def _peak_rbound_jit(x: NDArray[np.float64], peak: int) -> int:
    """Search for the right boundary starting from the peak apex."""
    n = x.size
    rbound = peak
    is_right_of_peak = False

    i = peak + 1
    while i < n:
        if x[i] < x[rbound]:
            rbound = i
            is_right_of_peak = True
        elif x[i] > x[rbound] and is_right_of_peak:
            cand = rbound
            rwindow = cand + 2
            if rwindow >= n:
                rwindow = n - 1

            i += 1
            while i <= rwindow and i < n:
                if x[i] < x[cand]:
                    rbound = i
                    break
                i += 1

            if cand == rbound:
                break

        i += 1

    if rbound < 0:
        return 0
    if rbound >= n:
        return n - 1
    return rbound


@jit(nopython=True, cache=True)
def _peak_boundaries_jit(
    x: NDArray[np.float64], peaks: NDArray[np.int64]
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Compute left and right boundaries for each peak."""
    n = peaks.size
    left = np.empty(n, dtype=np.int64)
    right = np.empty(n, dtype=np.int64)

    for i in range(n):
        p = int(peaks[i])
        left[i] = _peak_lbound_jit(x, p)
        right[i] = _peak_rbound_jit(x, p)

    return left, right


@jit(nopython=True, cache=True)
def _peak_bases_jit(
    x: NDArray[np.float64], peaks: NDArray[np.int64]
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Compute left/right base positions (the local minima on both sides of each peak)."""
    n = x.size
    m = peaks.size

    left = np.empty(m, dtype=np.int64)
    right = np.empty(m, dtype=np.int64)

    for i in range(m):
        p = int(peaks[i])

        l_idx = p
        j = p - 1
        while j >= 0:
            if x[j] > x[p]:
                break
            if x[j] < x[l_idx]:
                l_idx = j
            j -= 1

        r_idx = p
        j = p + 1
        while j < n:
            if x[j] > x[p]:
                break
            if x[j] < x[r_idx]:
                r_idx = j
            j += 1

        left[i] = l_idx
        right[i] = r_idx

    return left, right


@jit(nopython=True, cache=True)
def _trapz_jit(
    domain: NDArray[np.float64],
    signal: NDArray[np.float64],
    lower: int,
    upper: int,
) -> float:
    """Compute peak area using the trapezoidal integration rule."""
    if upper <= lower:
        return 0.0

    total = 0.0
    for i in range(lower + 1, upper + 1):
        dx = domain[i] - domain[i - 1]
        total += 0.5 * (signal[i] + signal[i - 1]) * dx
    return total


@jit(nopython=True, cache=True)
def _select_peaks_by_keep_jit(
    peaks: NDArray[np.int64], keep: NDArray[np.bool_]
) -> NDArray[np.int64]:
    """Select peak indices using the keep mask and return a compact index array."""
    count = 0
    for i in range(keep.size):
        if keep[i]:
            count += 1

    out = np.empty(count, dtype=np.int64)
    k = 0
    for i in range(keep.size):
        if keep[i]:
            out[k] = peaks[i]
            k += 1
    return out


@jit(nopython=True, cache=True)
def _findpeaks_filtered_indices_safe_jit(
    x: NDArray[np.float64],
    width: int,
    use_prominence: bool,
    prominence: float,
    use_relheight: bool,
    relheight: float,
    use_snr: bool,
    snr: float,
    noise_x: NDArray[np.float64],
) -> NDArray[np.int64]:
    """Run peak detection and multi-criteria filtering in numba, returning final peak indices."""
    peaks = _indices_from_bool_mask_jit(_local_maxima_mask_jit(x, width))
    m = peaks.size
    if m == 0:
        return peaks

    keep = np.ones(m, dtype=np.bool_)

    if use_prominence:
        left_bases, right_bases = _peak_bases_jit(x, peaks)
        for i in range(m):
            p = int(peaks[i])
            contour = x[left_bases[i]]
            rb = x[right_bases[i]]
            if rb > contour:
                contour = rb
            prom = x[p] - contour
            if (not np.isfinite(prom)) or prom < prominence:
                keep[i] = False

    if use_relheight:
        max_peak = -np.inf
        for i in range(m):
            if keep[i]:
                v = x[int(peaks[i])]
                if v > max_peak:
                    max_peak = v
        if (not np.isfinite(max_peak)) or max_peak <= 0.0:
            for i in range(m):
                keep[i] = False
        else:
            for i in range(m):
                if keep[i]:
                    rel = x[int(peaks[i])] / max_peak
                    if (not np.isfinite(rel)) or rel < relheight:
                        keep[i] = False

    # SNR filtering: compute peak-to-noise ratios and filter by the SNR threshold.
    if use_snr:
        if noise_x.size != x.size:
            for i in range(m):
                keep[i] = False
        else:
            for i in range(m):
                if not keep[i]:
                    continue

                p = int(peaks[i])
                nv = noise_x[p]
                if (not np.isfinite(nv)) or nv <= 0.0:
                    keep[i] = False
                    continue

                snr_val = x[p] / nv
                if (not np.isfinite(snr_val)) or snr_val < snr:
                    keep[i] = False

    return _select_peaks_by_keep_jit(peaks, keep)


@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def peak_pick_flat_parallel_core_jit(
    mz_arr: NDArray[np.float64],
    intensity_arr: NDArray[np.float64],
    lengths_arr: NDArray[np.int64],
    width: int,
    use_prominence: bool,
    prominence: float,
    use_relheight: bool,
    relheight: float,
    use_snr: bool,
    snr: float,
    return_area: bool,
    is_shared_mz: bool,
    noise_arr: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.int32]]:
    """
    Core parallel routine: detect/filter peaks for each spectrum and write
    results into compact output arrays.
    """
    n_spectra = lengths_arr.size
    input_offsets = lengths_to_offsets(lengths_arr)

    out_lengths = np.zeros(n_spectra, dtype=np.int32)

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        start = int(input_offsets[s])
        end = int(input_offsets[s + 1])
        if end <= start:
            continue

        seg = intensity_arr[start:end]
        noise_seg = noise_arr[start:end] if use_snr else np.empty(0, dtype=np.float64)
        peaks = _findpeaks_filtered_indices_safe_jit(
            x=seg,
            width=width,
            use_prominence=use_prominence,
            prominence=prominence,
            use_relheight=use_relheight,
            relheight=relheight,
            use_snr=use_snr,
            snr=snr,
            noise_x=noise_seg,
        )
        out_lengths[s] = np.int32(peaks.size)

    output_offsets = lengths_to_offsets(out_lengths.astype(np.int64))
    total_output = int(output_offsets[n_spectra])

    out_mz = np.empty(total_output, dtype=np.float64)
    out_values = np.empty(total_output, dtype=np.float64)

    for s in prange(n_spectra):  # pylint: disable=not-an-iterable
        in_start = int(input_offsets[s])
        in_end = int(input_offsets[s + 1])
        if in_end <= in_start:
            continue

        out_start = int(output_offsets[s])
        out_end = int(output_offsets[s + 1])
        if out_end <= out_start:
            continue

        segment = intensity_arr[in_start:in_end]
        if is_shared_mz:
            # continuous: all spectra in a batch share the same m/z axis.
            mz_segment = mz_arr[: (in_end - in_start)]
        else:
            # processed: each spectrum has its own m/z slice in the flat buffer.
            mz_segment = mz_arr[in_start:in_end]

        noise_segment = noise_arr[in_start:in_end] if use_snr else np.empty(0, dtype=np.float64)

        peaks = _findpeaks_filtered_indices_safe_jit(
            x=segment,
            width=width,
            use_prominence=use_prominence,
            prominence=prominence,
            use_relheight=use_relheight,
            relheight=relheight,
            use_snr=use_snr,
            snr=snr,
            noise_x=noise_segment,
        )

        if peaks.size == 0:
            continue

        if return_area:
            left, right = _peak_boundaries_jit(segment, peaks)

        for i in range(peaks.size):
            p = int(peaks[i])
            out_idx = out_start + i
            out_mz[out_idx] = mz_segment[p]
            if return_area:
                out_values[out_idx] = _trapz_jit(mz_segment, segment, int(left[i]), int(right[i])) # type: ignore
            else:
                out_values[out_idx] = segment[p]

    return out_mz, out_values, out_lengths
