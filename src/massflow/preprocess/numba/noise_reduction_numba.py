from typing import Optional
import numpy as np
from numba import njit, prange
from scipy.signal import savgol_coeffs
from massflow.tools.funs import prepare_flat_inputs, lengths_to_offsets


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



# ==================== Flat numba wrappers ====================


def smooth_signal_ma_numba(
    intensity: np.ndarray,
    window: int = 5,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Flat-mode moving-average smoothing."""
    if window < 1:
        raise ValueError("window must be >= 1")

    window = window + 1 if window % 2 == 0 else window
    _, intensity_arr, lengths_arr = prepare_flat_inputs(None, intensity, lengths)

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

    _, intensity_arr, lengths_arr = prepare_flat_inputs(None, intensity, lengths)

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

    _, intensity_arr, lengths_arr = prepare_flat_inputs(None, intensity, lengths)

    # Generate Gaussian kernel
    x = np.arange(-(window // 2), window // 2 + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (x / sd) ** 2)
    kernel /= kernel.sum()

    return convolve1d_flat_jit(intensity_arr, kernel, lengths_arr)
