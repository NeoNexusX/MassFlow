from typing import Optional
import numpy as np
from numba import njit, prange
from scipy.signal import savgol_coeffs
from scipy.spatial import cKDTree  # type: ignore
from scipy import stats
from massflow.tools.funs import prepare_flat_inputs, lengths_to_offsets
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.noise_reduction_numba")


# ==================== JIT-accelerated kernels ====================


@njit(fastmath=True, cache=True)
def _ma_core(
    signal: np.ndarray,
    window: int,
    valid_len: int,
    out: np.ndarray,
) -> None:

    if valid_len <= 0:
        return

    radius = window // 2
    inv_window = 1.0 / window
    left_edge = signal[0]
    right_edge = signal[valid_len - 1]

    # Initialize running sum for the first element
    current_sum = left_edge * np.float64(radius + 1)
    for j in range(1, radius + 1):
        current_sum += signal[j] if j < valid_len else right_edge

    out[0] = current_sum * inv_window

    # Slide the window
    for i in range(1, valid_len):
        leaving_idx = i - 1 - radius
        entering_idx = i + radius

        leaving_val = left_edge if leaving_idx < 0 else np.float64(signal[leaving_idx])
        entering_val = right_edge if entering_idx >= valid_len else np.float64(signal[entering_idx])

        current_sum += entering_val - leaving_val
        out[i] = current_sum * inv_window


@njit(parallel=True, cache=True, fastmath=True)
def _ma_flat_jit(flat: np.ndarray, window: int, lengths: np.ndarray) -> np.ndarray:
    """Flat batch entry point for ma_numba."""
    res = np.empty(flat.size, dtype=flat.dtype)
    offsets = lengths_to_offsets(lengths)

    for p in prange(lengths.size):  # pylint: disable=not-an-iterable
        start = offsets[p]
        end = offsets[p + 1]
        valid_len = end - start
        if valid_len > 0:
            _ma_core(flat[start:end], window, valid_len, res[start:end])

    return res


@njit(fastmath=True, cache=True)
def _savgol_1d_core(signal: np.ndarray,
                    kernels: np.ndarray
) -> np.ndarray:

    """Core 1D Savitzky-Golay routine handling position-dependent convolution kernels."""
    n = signal.size
    window = kernels.shape[0]
    half = window // 2
    res = np.empty(n, dtype=signal.dtype)
    for i in range(n):
        if i < half:
            pos = i
        elif i >= n - half:
            pos = window - (n - i)
        else:
            pos = half

        k_weights = kernels[pos]
        start = i - pos
        val = np.float64(0.0)
        for j in range(window):
            val += signal[start + j] * k_weights[j]
        res[i] = val
    return res


@njit(fastmath=True, cache=True, parallel=True)
def savgol_flat_jit(flat: np.ndarray, kernels: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Flat Savitzky-Golay entry point. Parallel over spectra segments."""
    res = np.empty(flat.size, dtype=flat.dtype)
    offsets = lengths_to_offsets(lengths)
    for p in prange(lengths.size):  # pylint: disable=not-an-iterable
        start = offsets[p]
        end = offsets[p + 1]
        if end > start:
            res[start:end] = _savgol_1d_core(flat[start:end], kernels)
    return res


@njit(fastmath=True, cache=True)
def _convolve1d_core_edge(signal: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Core 1D convolution routine with edge padding."""
    n = signal.size
    k_len = kernel.size
    half = k_len // 2
    res = np.empty(n, dtype=signal.dtype)

    for i in range(n):
        val = np.float64(0.0)
        for j in range(k_len):
            idx = i - half + j
            if idx < 0:
                idx = 0
            elif idx >= n:
                idx = n - 1
            val += signal[idx] * kernel[j]
        res[i] = val
    return res


@njit(fastmath=True, cache=True, parallel=True)
def convolve1d_flat_jit(flat: np.ndarray, kernel: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Flat convolution entry point. Parallel over spectra segments."""
    res = np.empty(flat.size, dtype=flat.dtype)
    offsets = lengths_to_offsets(lengths)
    for p in prange(lengths.size):  # pylint: disable=not-an-iterable
        start = offsets[p]
        end = offsets[p + 1]
        if end > start:
            res[start:end] = _convolve1d_core_edge(flat[start:end], kernel)
    return res


@njit(fastmath=True, cache=True, parallel=True)
def _ns_gaussian_kernel(
    neigh_intensity: np.ndarray, dists: np.ndarray, sd: float
) -> np.ndarray:
    """Numba kernel for NS Gaussian smoothing over k-nearest neighbours (float64 compute)."""
    n, k = neigh_intensity.shape
    out = np.empty(n, dtype=np.float64)
    sd64 = np.float64(sd)
    for i in prange(n): # pylint: disable=not-an-iterable
        row_sum = 0.0
        s = 0.0
        for j in range(k):
            exponent = -0.5 * (dists[i, j] / sd64) ** 2
            exponent = min(max(exponent, -88.0), 0.0)
            w = np.exp(exponent)
            row_sum += w
            s += neigh_intensity[i, j] * w
        if row_sum == 0.0:
            out[i] = 0.0
        else:
            out[i] = s / row_sum
    return out


@njit(fastmath=True, cache=True, parallel=True)
def _ns_bilateral_kernel(
    neigh_intensity: np.ndarray,
    dists: np.ndarray,
    center_intensity: np.ndarray,
    sd_dist: float,
    sd_intensity: float,
) -> np.ndarray:
    """Numba kernel for NS bilateral filtering (float64 compute)."""
    n, k = neigh_intensity.shape
    out = np.empty(n, dtype=np.float64)
    sd_dist64 = np.float64(sd_dist)
    sd_int64 = np.float64(sd_intensity)
    for i in prange(n):# pylint: disable=not-an-iterable
        ci = center_intensity[i]
        row_sum = 0.0
        s = 0.0
        for j in range(k):
            exponent_d = -0.5 * (dists[i, j] / sd_dist64) ** 2
            exponent_d = min(max(exponent_d, -88.0), 0.0)
            wd = np.exp(exponent_d)
            diff = neigh_intensity[i, j] - ci
            exponent_i = -0.5 * (diff / sd_int64) ** 2
            exponent_i = min(max(exponent_i, -88.0), 0.0)
            wi = np.exp(exponent_i)
            w = wd * wi
            row_sum += w
            s += neigh_intensity[i, j] * w
        if row_sum == 0.0:
            out[i] = 0.0
        else:
            out[i] = s / row_sum
    return out


@njit(cache=True, fastmath=True)
def _lowest_numba(
    x: np.ndarray,
    y: np.ndarray,
    xs: float,
    nleft: int,
    nright: int,
    w: np.ndarray,
    userw: bool,
    rw: np.ndarray,
    range_x: float,
) -> tuple[float, bool, int]:
    """Local weighted fit at a single x position (LOWESS lowest step)."""
    n = x.size
    h = max(xs - x[nleft], x[nright] - xs)
    h9 = 0.999 * h
    h1 = 0.001 * h

    weight_sum = 0.0
    j = nleft
    while j < n:
        wj = 0.0
        r = abs(x[j] - xs)
        if r <= h9:
            if r <= h1:
                wj = 1.0
            else:
                t = 1.0 - (r / h) ** 3
                wj = t * t * t
            if userw:
                wj *= rw[j]
            weight_sum += wj
        elif x[j] > xs:
            break
        w[j] = wj
        j += 1

    nrt = j - 1
    if weight_sum <= 0.0:
        return 0.0, False, nrt

    for j in range(nleft, nrt + 1):
        w[j] /= weight_sum

    if h > 0.0:
        x_center = 0.0
        for j in range(nleft, nrt + 1):
            x_center += w[j] * x[j]
        b = xs - x_center
        c = 0.0
        for j in range(nleft, nrt + 1):
            dx = x[j] - x_center
            c += w[j] * dx * dx
        if np.sqrt(c) > 0.001 * range_x:
            b /= c
            for j in range(nleft, nrt + 1):
                w[j] *= b * (x[j] - x_center) + 1.0

    ys = 0.0
    for j in range(nleft, nrt + 1):
        ys += w[j] * y[j]
    return ys, True, nrt


@njit(cache=True, fastmath=True)
def _lowess_core(
    x: np.ndarray,
    y: np.ndarray,
    f: float,
    nsteps: int,
    delta: float,
) -> np.ndarray:
    """Numba LOWESS core (lowest + clowess pipeline)."""
    n = x.size
    ys = np.empty(n, dtype=np.float64)
    if n == 0:
        return ys
    if n == 1:
        ys[0] = y[0]
        return ys

    ns = max(2, min(n, int(f * n + 1.0e-7)))
    rw = np.ones(n, dtype=np.float64)
    res = np.empty(n, dtype=np.float64)
    w = np.zeros(n, dtype=np.float64)
    range_x = x[n - 1] - x[0]

    it = 0
    while it <= nsteps:
        nleft = 0
        nright = ns - 1
        last = -1
        i = 0

        while True:
            if nright < n - 1:
                d1 = x[i] - x[nleft]
                d2 = x[nright + 1] - x[i]
                if d1 > d2:
                    nleft += 1
                    nright += 1
                    continue

            yi, ok, _ = _lowest_numba(
                x, y, x[i], nleft, nright, w, it > 0, rw, range_x
            )
            ys[i] = yi if ok else y[i]

            if last < i - 1:
                denom = x[i] - x[last]
                if denom != 0.0:
                    for j in range(last + 1, i):
                        alpha = (x[j] - x[last]) / denom
                        ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last]
                else:
                    for j in range(last + 1, i):
                        ys[j] = ys[i]

            last = i
            cut = x[last] + delta
            i_scan = last + 1
            while i_scan < n:
                if x[i_scan] > cut:
                    break
                if x[i_scan] == x[last]:
                    ys[i_scan] = ys[last]
                    last = i_scan
                i_scan += 1

            next_i = i_scan - 1
            i = last + 1 if (last + 1) > next_i else next_i
            if last >= n - 1:
                break

        for k in range(n):
            res[k] = y[k] - ys[k]

        if it >= nsteps:
            break

        sc = 0.0
        for k in range(n):
            rk = abs(res[k])
            rw[k] = rk
            sc += rk
        sc /= n

        cmad = 6.0 * np.median(rw)
        if cmad < 1.0e-7 * sc:
            break

        c9 = 0.999 * cmad
        c1 = 0.001 * cmad
        for k in range(n):
            rk = abs(res[k])
            if rk <= c1:
                rw[k] = 1.0
            elif rk <= c9:
                t = 1.0 - (rk / cmad) ** 2
                rw[k] = t * t
            else:
                rw[k] = 0.0

        it += 1

    return ys



# ==================== Non-JIT wrappers and preprocessing ====================

def smooth_lowess_numba(
    x: Optional[np.ndarray] = None,
    y: Optional[np.ndarray] = None,
    f: float = 2.0 / 3.0,
    iter_count: int = 3,
    delta: Optional[float] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    LOWESS wrapper compatible with R-style arguments.
    Returns sorted x and smoothed y aligned to sorted x.
    
    Parameters:
        x: x-coordinates. If None, automatically generated as np.arange(y.size).
        y: y-coordinates (required).
        f: smoothing parameter (0 < f <= 1).
        iter_count: number of robustness iterations.
        delta: distance to consider points as tied in x.
    """
    if y is None:
        raise ValueError("y must not be None")

    # Auto-generate x if not provided
    if x is None:
        x = np.arange(y.size, dtype=np.float64)

    x_dtype = x.dtype
    y_dtype = y.dtype

    if x.size == 0:
        return x, y
    if (x.ndim != 1
        or y.ndim != 1
        or x.size != y.size
        or not np.isfinite(f)
        or f <= 0.0
        or not isinstance(iter_count, (int, np.integer))
        or iter_count < 0
    ):
        raise ValueError(
            "x and y must be 1D arrays with same length; "
            "f must be finite and > 0; "
            "iter_count (R iter) must be an integer >= 0"
        )

    # Keep original input dtypes for output, but run LOWESS in float64.
    x64 = np.ascontiguousarray(x.astype(np.float64, copy=False))
    y64 = np.ascontiguousarray(y.astype(np.float64, copy=False))

    # order = np.argsort(x_arr)
    # x_arr = np.ascontiguousarray(x_arr[order])
    # y_arr = np.ascontiguousarray(y_arr[order])

    if delta is None:
        delta_val = 0.01 * (x64[-1] - x64[0]) if x64.size > 1 else 0.0
    else:
        delta_val = float(delta)
        if not np.isfinite(delta_val) or delta_val < 0.0:
            raise ValueError("delta must be finite and >= 0")

    y_fit64 = _lowess_core(x64, y64, float(f), int(iter_count), float(delta_val))
    return x.astype(x_dtype, copy=False), y_fit64.astype(y_dtype, copy=False)


def smooth_signal_ma_numba(
    intensity: np.ndarray,
    window: int = 5,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Flat-mode moving-average smoothing."""
    if window < 1:
        raise ValueError("window must be >= 1")

    window = window + 1 if window % 2 == 0 else window
    intensity_arr, lengths_arr = prepare_flat_inputs(intensity, lengths)

    return _ma_flat_jit(intensity_arr, window, lengths_arr)


def smooth_signal_savgol_numba(
    intensity: np.ndarray,
    window: int = 5,
    polyorder: int = 3,
    deriv: int = 0,
    delta: float = 1.0,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Flat-mode Savitzky-Golay smoothing."""
    actual_window = max(3, window + (1 - window % 2))
    if polyorder >= actual_window:
        raise ValueError("polyorder must be < window")

    intensity_arr, lengths_arr = prepare_flat_inputs(intensity, lengths)

    kernels = np.zeros((actual_window, actual_window), dtype=np.float64)
    for pos in range(actual_window):
        kernels[pos] = savgol_coeffs(
            actual_window, polyorder, deriv=deriv, delta=delta, pos=pos, use="dot"
        ).astype(np.float64)

    return savgol_flat_jit(intensity_arr, kernels, lengths_arr)


def smooth_signal_gaussian_numba(
    intensity: np.ndarray,
    window: int = 5,
    sd: Optional[float] = None,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Flat-mode Gaussian smoothing."""
    if window < 1:
        raise ValueError("window must be >= 1")

    # Ensure window is odd
    if window % 2 == 0:
        window += 1

    if sd is None:
        sd = window / 4.0
    if sd <= 0:
        raise ValueError("sd must be positive")

    intensity_arr, lengths_arr = prepare_flat_inputs(intensity, lengths)

    # Generate Gaussian kernel
    x = np.arange(-(window // 2), window // 2 + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (x / sd) ** 2).astype(np.float64)
    kernel /= kernel.sum()

    return convolve1d_flat_jit(intensity_arr, kernel, lengths_arr)


def smooth_ns_signal_gaussian_numba(
    intensity: np.ndarray,
    index: np.ndarray,
    window: int = 5,
    p: int = 1,
    sd: float | None = None,
) -> np.ndarray:
    """Neighbourhood-smoothing Gaussian (Numba) entry point.

    Estimates a default spatial scale from neighbour distances when sd is
    None, then applies the parallel Gaussian kernel.
    """
    input_dtype = np.asarray(intensity).dtype
    neigh_intensity, dists, _ = ns_signal_pre(intensity, index, k=window, p=p)
    dists_max = np.max(dists, axis=1)
    sd_val = np.median(dists_max) / 2.0 if sd is None else sd
    if sd_val <= 0:
        raise ValueError("sd must be positive")
    result64 = _ns_gaussian_kernel(neigh_intensity, dists, float(sd_val))
    return result64.astype(input_dtype)


def smooth_ns_signal_bi_numba(
    intensity: np.ndarray,
    index: np.ndarray,
    window: int = 5,
    p: int = 2,
    sd: float | None = None,
    sd_intensity: float | None = None,
) -> np.ndarray:
    """Neighbourhood-smoothing bilateral (Numba) entry point.

    Uses spatial distances to derive sd_dist (if not provided) and robust
    MAD-based estimate for sd_intensity, then applies the parallel bilateral
    kernel.
    """
    input_dtype = np.asarray(intensity).dtype
    neigh_intensity, dists, _ = ns_signal_pre(intensity, index, k=window, p=p)
    dists_max = np.max(dists, axis=1)
    sd_dist_val = np.median(dists_max) / 2.0 if sd is None else sd
    if sd_dist_val <= 0:
        raise ValueError("sd_dist must be positive")
    center = np.asarray(intensity, dtype=np.float64)
    sd_int_val = stats.median_abs_deviation(center, nan_policy="omit", scale="normal") if sd_intensity is None else sd_intensity  # type: ignore
    if sd_int_val <= 0:
        raise ValueError("sd_intensity must be positive")
    result64 = _ns_bilateral_kernel(
        neigh_intensity, dists, center, float(sd_dist_val), float(sd_int_val)
    )
    return result64.astype(input_dtype)


def ns_signal_pre(
    intensity: np.ndarray,
    index: np.ndarray,
    k: int = 5,
    p: int = 1,
):
    """
    Prepare neighborhood search (kNN over `index`) for NS-based smoothing.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity` (e.g., m/z).
        k (int): Number of nearest neighbors used for smoothing. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
            - neigh_intensity: shape (N, k) neighbor intensities
            - dists: shape (N, k) neighbor distances (Minkowski p)
            - idxs: shape (N, k) neighbor indices

    Raises:
        ValueError: If `k` or `p` is not positive.
        ValueError: If `index` is None or length mismatches `intensity`.
    """

    # Enforce positivity
    if k < 1 or p < 1:
        logger.error("k and p must be positive integers")
        raise ValueError("k and p must be positive integers")

    if index is None or len(index) != len(intensity):
        logger.warning("mz_list must be provided and match intensity length; using np.arange as fallback mz_list")
        index = np.arange(len(intensity))
    else:
        logger.debug(f"using provided mz_list for neighborhood search with k={k} and p={p}\r\n{index[:5]}")

    if not isinstance(k, int):
        logger.error("k must be an integer")
        raise ValueError("k must be an integer")
    if not isinstance(p, int):
        logger.error("p must be an integer")
        raise ValueError("p must be an integer")

    tree = cKDTree(index.reshape(-1, 1))

    if len(intensity) < k:
        logger.warning("spectrum length must be greater than k")
        k = len(intensity)

    dists, idxs = tree.query(index.reshape(-1, 1), k=k, p=p)

    # Ensure idxs is (N, k) for consistent neighbor aggregation
    if np.ndim(idxs) == 1:
        idxs = idxs.reshape(-1, 1)

    # Impute NaNs in intensity using mean of its k neighbors
    nan_idx = np.where(np.isnan(intensity))[0]
    if nan_idx.size > 0:
        # shape: (M, k)
        neigh = intensity[idxs[nan_idx]]
        fill_vals = np.nanmean(neigh, axis=1)
        # Fallback for rows whose neighbors are all NaN
        global_median = np.nanmedian(intensity)
        if np.isnan(global_median):
            global_median = 0.0
        fill_vals = np.where(np.isnan(fill_vals), global_median, fill_vals)
        intensity[nan_idx] = fill_vals

    neigh_intensity = intensity[idxs]  # shape: (N, k)

    return  neigh_intensity, dists, idxs
