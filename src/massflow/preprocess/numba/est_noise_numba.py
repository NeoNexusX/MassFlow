"""Numba-accelerated helpers for noise estimation."""

import numpy as np
from numba import njit, prange
from massflow.tools.funs import prepare_flat_inputs, lengths_to_offsets


@njit(cache=True)
def _finite_compact_numba(data: np.ndarray) -> np.ndarray:
    """Return original data when no NaN, else return NaN-filtered float64 values."""
    n = data.size
    count = 0
    for i in range(n):
        v = data[i]
        if not np.isnan(v):
            count += 1

    if count == n:
        return data

    out = np.empty(count, dtype=np.float64)
    j = 0
    for i in range(n):
        v = data[i]
        if not np.isnan(v):
            out[j] = float(v)
            j += 1
    return out


@njit(cache=True)
def _nanmean_numba(data: np.ndarray) -> float:
    vals = _finite_compact_numba(data)
    n_vals = vals.size
    if n_vals == 0:
        return np.nan
    return np.sum(vals) / n_vals


@njit(cache=True)
def _nanstd_numba(data: np.ndarray) -> float:
    vals = _finite_compact_numba(data)
    n_vals = vals.size
    if n_vals == 0:
        return np.nan
    mu = np.sum(vals) / n_vals
    s = 0.0
    for i in range(n_vals):
        d = vals[i] - mu
        s += d * d
    return np.sqrt(s / n_vals)


@njit(cache=True)
def _nanmedian_numba(data: np.ndarray) -> float:
    vals = _finite_compact_numba(data)
    if vals.size == 0:
        return np.nan
    vals_sorted = np.sort(vals)
    m = vals_sorted.size
    mid = m // 2
    if m % 2 == 1:
        return vals_sorted[mid]
    return 0.5 * (vals_sorted[mid - 1] + vals_sorted[mid])


@njit(cache=True)
def _nanquantile_numba(data: np.ndarray, q: float) -> float:
    """Linear-interpolated quantile over finite values."""
    vals = _finite_compact_numba(data)
    n = vals.size
    if n == 0:
        return np.nan

    if q <= 0.0:
        return np.min(vals)
    if q >= 1.0:
        return np.max(vals)

    vals_sorted = np.sort(vals)
    pos = (n - 1) * q
    lo = int(np.floor(pos))
    hi = int(np.ceil(pos))
    if lo == hi:
        return vals_sorted[lo]
    w = pos - lo
    return vals_sorted[lo] * (1.0 - w) + vals_sorted[hi] * w


@njit(cache=True)
def _mad_numba(data: np.ndarray) -> float:
    vals = _finite_compact_numba(data)
    if vals.size == 0:
        return np.nan
    if vals.size != data.size:
        return np.nan

    med = _nanmedian_numba(vals)
    if np.isnan(med):
        return np.nan
    dev = np.empty(vals.size, dtype=np.float64)
    for i in range(vals.size):
        dev[i] = abs(vals[i] - med)
    return _nanmedian_numba(dev)


@njit(cache=True)
def estimation_fun_numba(data: np.ndarray, method_code: int = 0) -> float:
    """Numba version of estimation_fun.

    method_code:
    - 0: sd
    - 1: mad
    - 2: quantile(0.95)
    - 3: diff (mean absolute deviation from mean)
    """
    if method_code == 0:
        return _nanstd_numba(data)
    if method_code == 1:
        return _mad_numba(data)
    if method_code == 2:
        return _nanquantile_numba(data, 0.95)
    if method_code == 3:
        mu = _nanmean_numba(data)
        if np.isnan(mu):
            return np.nan
        vals = _finite_compact_numba(data)
        if vals.size == 0:
            return np.nan
        s = 0.0
        for i in range(vals.size):
            s += abs(vals[i] - mu)
        return s / vals.size

    raise ValueError("Unknown noise estimation method code")


@njit(cache=True)
def linear_sse_numba(data: np.ndarray, i: int, j: int) -> float:
    """Compute SSE of a first-degree polynomial fit in interval [i, j].

    This is the Numba counterpart of ``_linear_sse`` in
    ``preprocess/helper/est_noise_helper.py``.
    """
    if i > j:
        raise ValueError("Invalid interval: i > j")

    n = data.size
    if n == 0:
        return 0.0

    i0 = 0 if i < 0 else i
    j0 = n - 1 if j >= n else j
    if i0 > j0:
        return 0.0

    # First pass: gather means on valid (non-NaN) samples.
    count = 0
    sum_x = 0.0
    sum_y = 0.0
    for idx in range(i0, j0 + 1):
        yi = data[idx]
        if not np.isnan(yi):
            count += 1
            sum_x += float(idx)
            sum_y += float(yi)

    if count == 0:
        return 0.0

    x_mean = sum_x / count
    y_mean = sum_y / count

    # Second pass: centered sums for slope denominator/numerator.
    denom = 0.0
    numer = 0.0
    for idx in range(i0, j0 + 1):
        yi = data[idx]
        if not np.isnan(yi):
            dx = float(idx) - x_mean
            dy = float(yi) - y_mean
            denom += dx * dx
            numer += dx * dy

    # Degenerate case: one-point (or all x collapsed) -> fit mean only.
    if denom == 0.0:
        sse = 0.0
        for idx in range(i0, j0 + 1):
            yi = data[idx]
            if not np.isnan(yi):
                r = float(yi) - y_mean
                sse += r * r
        return sse

    slope = numer / denom
    intercept = y_mean - slope * x_mean

    sse = 0.0
    for idx in range(i0, j0 + 1):
        yi = data[idx]
        if not np.isnan(yi):
            y_hat = slope * float(idx) + intercept
            r = float(yi) - y_hat
            sse += r * r
    return sse


@njit(cache=True)
def bins_sse_numba(data: np.ndarray, lowers: np.ndarray, uppers: np.ndarray) -> np.ndarray:
    """Compute per-bin SSE in float64 by calling ``linear_sse_numba``.

    Expected dtypes:
    - data: float64 (1D)
    - lowers/uppers: int64 (1D)

    Returns:
    - float64 array with one SSE value per bin.
    """
    n_bins = lowers.size
    out = np.empty(n_bins, dtype=np.float64)
    for k in range(n_bins):
        out[k] = linear_sse_numba(data, int(lowers[k]), int(uppers[k]))
    return out


@njit(cache=True)
def best_split_index_numba(data: np.ndarray, i: int, j: int) -> int:
    """Numba version of ``_best_split_index``.

    Finds split t in [i, j] minimizing
    linear_sse_numba(data, i, t - 1) + linear_sse_numba(data, t, j).
    Returns -1 when no valid split exists.
    """
    if j - i + 1 < 2:
        return -1

    best_t = -1
    best_score = np.inf
    for t in range(i + 1, j + 1):
        s_left = linear_sse_numba(data, i, t - 1)
        s_right = linear_sse_numba(data, t, j)
        score = s_left + s_right
        if score < best_score:
            best_score = score
            best_t = t
    return best_t


@njit(cache=True)
def merge_best_pair_numba(
    lowers: np.ndarray,
    uppers: np.ndarray,
    sse: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Numba version of merge-best adjacent bins."""
    m = lowers.size
    if m <= 1:
        return lowers.copy(), uppers.copy()

    p = 0
    best = sse[0] + sse[1]
    for idx in range(1, m - 1):
        score = sse[idx] + sse[idx + 1]
        if score < best:
            best = score
            p = idx

    new_lowers = np.empty(m - 1, dtype=np.int64)
    new_uppers = np.empty(m - 1, dtype=np.int64)

    for idx in range(p):
        new_lowers[idx] = lowers[idx]
        new_uppers[idx] = uppers[idx]

    new_lowers[p] = lowers[p]
    new_uppers[p] = uppers[p + 1]

    for idx in range(p + 1, m - 1):
        new_lowers[idx] = lowers[idx + 1]
        new_uppers[idx] = uppers[idx + 1]

    return new_lowers, new_uppers


@njit(cache=True)
def split_worst_bin_numba(
    data: np.ndarray,
    lowers: np.ndarray,
    uppers: np.ndarray,
    sse: np.ndarray,
) -> tuple[bool, np.ndarray, np.ndarray]:
    """Numba version of split-worst-bin.

    Returns (ok, new_lowers, new_uppers). If split fails, ok=False and
    original bounds are returned.
    """
    m = lowers.size
    if m == 0:
        return False, lowers, uppers

    order = np.argsort(-sse)
    for oi in range(order.size):
        k = int(order[oi])
        i = int(lowers[k])
        j = int(uppers[k])
        t = best_split_index_numba(data, i, j)
        if t != -1:
            new_lowers = np.empty(m + 1, dtype=np.int64)
            new_uppers = np.empty(m + 1, dtype=np.int64)

            for idx in range(k):
                new_lowers[idx] = lowers[idx]
                new_uppers[idx] = uppers[idx]

            new_lowers[k] = i
            new_uppers[k] = t - 1
            new_lowers[k + 1] = t
            new_uppers[k + 1] = j

            for idx in range(k + 1, m):
                new_lowers[idx + 1] = lowers[idx]
                new_uppers[idx + 1] = uppers[idx]

            return True, new_lowers, new_uppers

    return False, lowers, uppers


@njit(cache=True)
def findbins_numba(
    data: np.ndarray,
    nbins: int = 1,
    dynamic: bool = False,
    niter: int = 10,
    overlap: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Numba version of ``_findbins``.

    Returns Numba-friendly tuple:
    (lower, upper, size, sse, trace)

    Notes:
    - data should be float64 1D for best compatibility.
    - lowers/uppers/size are int64 arrays.
    - sse/trace are float64 arrays; empty when not applicable.
    """
    if data.ndim != 1:
        raise ValueError("Input data must be a 1D array")

    if nbins < 1:
        nbins = 1
    if overlap < 0.0 or overlap >= 0.99:
        overlap = 0.5

    n = data.size
    if n == 0:
        empty_i64 = np.empty(0, dtype=np.int64)
        empty_f64 = np.empty(0, dtype=np.float64)
        return empty_i64, empty_i64, empty_i64, empty_f64, empty_f64

    nbins = nbins if nbins < n else n

    sse = np.empty(0, dtype=np.float64)
    trace = np.empty(0, dtype=np.float64)

    if overlap > 0.0 and not dynamic:
        width = n / (nbins * (1.0 - overlap) + overlap)
        step = width * (1.0 - overlap)
        starts = np.arange(nbins, dtype=np.float64) * step
        lower = np.floor(starts).astype(np.int64)
        upper = (np.ceil(starts + width).astype(np.int64) - 1)
        lower = np.clip(lower, 0, n - 1)
        upper = np.clip(upper, 0, n - 1)

    elif overlap == 0.0:
        edges = np.floor(np.linspace(0.0, float(n), nbins + 1)).astype(np.int64)
        lower = np.clip(edges[:-1], 0, n - 1)
        upper = np.clip(edges[1:] - 1, 0, n - 1)

    elif dynamic:
        if nbins < 3:
            nbins = 3

        edges = np.floor(np.linspace(0.0, float(n), nbins + 1)).astype(np.int64)
        lower = np.clip(edges[:-1], 0, n - 1)
        upper = np.clip(edges[1:] - 1, 0, n - 1)

        sse = bins_sse_numba(data, lower, upper)
        max_steps = int(max(1, niter))
        trace_buf = np.empty(max_steps + 1, dtype=np.float64)
        trace_len = 0
        trace_buf[trace_len] = np.sum(sse)
        trace_len += 1

        for _it in range(max_steps):
            merged_lower, merged_upper = merge_best_pair_numba(lower, upper, sse)
            merged_sse = bins_sse_numba(data, merged_lower, merged_upper)

            ok, split_lower, split_upper = split_worst_bin_numba(
                data, merged_lower, merged_upper, merged_sse
            )
            if not ok:
                break

            new_sse = bins_sse_numba(data, split_lower, split_upper)
            new_score = np.sum(new_sse)

            if new_score < trace_buf[trace_len - 1]:
                lower = split_lower
                upper = split_upper
                sse = new_sse
                trace_buf[trace_len] = new_score
                trace_len += 1
            else:
                trace_buf[trace_len] = trace_buf[trace_len - 1]
                trace_len += 1
                break

        trace = trace_buf[:trace_len]

    else:
        # Should not be reached after overlap clipping, keep a safe fallback.
        edges = np.floor(np.linspace(0.0, float(n), nbins + 1)).astype(np.int64)
        lower = np.clip(edges[:-1], 0, n - 1)
        upper = np.clip(edges[1:] - 1, 0, n - 1)

    size = (upper - lower + 1).astype(np.int64)
    return lower, upper, size, sse, trace


@njit(cache=True)
def _compress_midpoints_numba(midpoints: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sort by x and merge duplicated x by averaging y."""
    n = midpoints.size
    if n == 0:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)

    order = np.argsort(midpoints)
    xs = np.empty(n, dtype=np.float64)
    ys = np.empty(n, dtype=np.float64)
    for i in range(n):
        idx = int(order[i])
        xs[i] = midpoints[idx]
        ys[i] = values[idx]

    out_x = np.empty(n, dtype=np.float64)
    out_y = np.empty(n, dtype=np.float64)

    out_n = 0
    cur_x = xs[0]
    sum_y = ys[0]
    cnt = 1

    for i in range(1, n):
        if xs[i] == cur_x:
            sum_y += ys[i]
            cnt += 1
        else:
            out_x[out_n] = cur_x
            out_y[out_n] = sum_y / cnt
            out_n += 1
            cur_x = xs[i]
            sum_y = ys[i]
            cnt = 1

    out_x[out_n] = cur_x
    out_y[out_n] = sum_y / cnt
    out_n += 1
    return out_x[:out_n], out_y[:out_n]


@njit(cache=True)
def _lagrange_eval_numba(xp: np.ndarray, yp: np.ndarray, xq: float) -> float:
    m = xp.size
    out = 0.0
    for i in range(m):
        term = yp[i]
        xi = xp[i]
        for j in range(m):
            if j != i:
                denom = xi - xp[j]
                if denom == 0.0:
                    return np.nan
                term *= (xq - xp[j]) / denom
        out += term
    return out


@njit(cache=True)
def _local_poly_interp_numba(x: np.ndarray, y: np.ndarray, xq: np.ndarray, k: int) -> np.ndarray:
    """Evaluate interpolating local polynomial of degree k on query points."""
    n = x.size
    out = np.empty(xq.size, dtype=np.float64)

    if n == 0:
        out[:] = np.nan
        return out
    if n == 1:
        out[:] = y[0]
        return out

    degree = k
    if degree < 1:
        degree = 1
    if degree > n - 1:
        degree = n - 1
    span = degree + 1

    for qi in range(xq.size):
        xv = xq[qi]

        # insertion index in sorted x
        pos = n
        for p in range(n):
            if x[p] >= xv:
                pos = p
                break

        left = pos - span // 2
        if left < 0:
            left = 0
        if left > n - span:
            left = n - span

        xp = x[left:left + span]
        yp = y[left:left + span]
        out[qi] = _lagrange_eval_numba(xp, yp, xv)

    return out


@njit(cache=True)
def estimate_noise_curve_numba(
    residuals: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    method_code: int,
    max_degree: int = 3,
) -> np.ndarray:
    """Full numba path: per-bin estimation + interpolated noise curve."""
    m = lower.size
    mids = np.empty(m, dtype=np.float64)
    vals = np.empty(m, dtype=np.float64)

    for i in range(m):
        li = int(lower[i])
        ui = int(upper[i])
        mids[i] = 0.5 * (li + ui)
        vals[i] = estimation_fun_numba(residuals[li:ui + 1], method_code)

    xk, yk = _compress_midpoints_numba(mids, vals)
    if xk.size == 0:
        return np.zeros(residuals.size, dtype=np.float64)
    if xk.size == 1:
        out = np.empty(residuals.size, dtype=np.float64)
        out[:] = yk[0]
        return out

    degree = max_degree
    if degree > xk.size - 1:
        degree = xk.size - 1
    if degree < 1:
        degree = 1

    xq = np.arange(residuals.size, dtype=np.float64)
    return _local_poly_interp_numba(xk, yk, xq, degree)


@njit(cache=True, parallel=True)
def _estimate_flat_core_numba(
    flat: np.ndarray,
    lengths: np.ndarray,
    nbins: int,
    dynamic: bool,
    niter: int,
    overlap: float,
    method_code: int,
    floor_value: float,
) -> np.ndarray:
    out = np.empty(flat.size, dtype=np.float64)
    offsets = lengths_to_offsets(lengths)

    for p in prange(lengths.size):  # pylint: disable=not-an-iterable
        start = int(offsets[p])
        end = int(offsets[p + 1])
        seg_len = end - start
        if seg_len <= 0:
            continue

        segment = flat[start:end]
        if nbins > 1:
            lower, upper, _size, _sse, _trace = findbins_numba(
                segment,
                nbins=nbins,
                dynamic=dynamic,
                niter=niter,
                overlap=overlap,
            )
            curve = estimate_noise_curve_numba(segment, lower, upper, method_code, max_degree=3)
            for i in range(seg_len):
                v = curve[i]
                out[start + i] = v if v >= floor_value else floor_value
        else:
            val = estimation_fun_numba(segment, method_code)
            if np.isnan(val) or val < floor_value:
                val = floor_value
            for i in range(seg_len):
                out[start + i] = val

    return out

def estimate_flat_numba(
    intensity: np.ndarray,
    lengths: np.ndarray | None = None,
    nbins: int = 1,
    dynamic: bool = False,
    niter: int = 10,
    overlap: float = 0.5,
    method_code: int  = 0,
    floor_value: float = 0.001,
) -> np.ndarray:
    """Flat-mode parallel noise estimation.

    - Input is a flat 1D signal array.
    - Segmentation must be provided by ``lengths``.
    - ``nbins > 1`` returns piecewise-interpolated noise curve per spectrum.
    - ``nbins == 1`` returns per-spectrum scalar noise repeated on that spectrum.
    """
    if lengths is None:
        raise ValueError("lengths must be provided for flat-mode noise estimation")

    intensity_arr, lengths_arr = prepare_flat_inputs(
        np.asarray(intensity),
        np.asarray(lengths, dtype=np.int64),
    )

    return _estimate_flat_core_numba(
        intensity_arr,
        lengths_arr,
        int(nbins),
        bool(dynamic),
        int(niter),
        float(overlap),
        int(method_code),
        float(floor_value),
    )
