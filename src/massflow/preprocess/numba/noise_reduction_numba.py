import numpy as np
from numba import jit, prange, set_num_threads
from scipy.signal import savgol_coeffs
from scipy.spatial import cKDTree
from scipy import stats
from massflow.tools.logger import get_logger

logger = get_logger("noise_reduction_numba")

# ---------------------------------------------------------------------
# Moving Average (ma_loop)
# ---------------------------------------------------------------------

@jit(nopython=True, fastmath=True, cache=True)
def _ma_loop_core_edge(signal: np.ndarray, window: int) -> np.ndarray:
    """
    Core Moving Average routine using an explicit sliding window loop (O(N)).
    Simulates mode='edge' padding manually.
    """
    n = signal.size
    res = np.empty(n, dtype=np.float32)
    radius = window // 2

    # Initialize running sum for the first element (index 0)
    current_sum = 0.0

    # Calculate initial sum centered at 0: range [-radius, radius]
    for j in range(-radius, radius + 1):
        idx = j
        if idx < 0:
            val = signal[0]
        elif idx >= n:
            val = signal[n - 1]
        else:
            val = signal[idx]
        current_sum += val

    res[0] = current_sum / window

    # Slide the window
    for i in range(1, n):
        leaving_idx = i - 1 - radius
        entering_idx = i + radius

        if leaving_idx < 0:
            leaving_val = signal[0]
        else:
            leaving_val = signal[leaving_idx]

        if entering_idx >= n:
            entering_val = signal[n - 1]
        else:
            entering_val = signal[entering_idx]

        current_sum = current_sum - leaving_val + entering_val
        res[i] = current_sum / window

    return res


@jit(nopython=True, parallel=True, cache=True)
def _ma_loop_batch_jit(intensity: np.ndarray, window: int, lengths: np.ndarray = None) -> np.ndarray:
    """Batch entry point for ma_loop."""
    if intensity.ndim == 1:
        return _ma_loop_core_edge(intensity, window)

    n_pixels, n_mz = intensity.shape
    res = np.empty((n_pixels, n_mz), dtype=np.float32)

    if lengths is None:
        for p in prange(n_pixels):
            res[p] = _ma_loop_core_edge(intensity[p], window)
    else:
        for p in prange(n_pixels):
            valid_len = lengths[p]
            processed = _ma_loop_core_edge(intensity[p, :valid_len], window)
            res[p, :valid_len] = processed

    return res


def smooth_signal_ma_loop(
    intensity: np.ndarray,
    window: int = 5,
    lengths: np.ndarray = None,
    numba_max_threads: int = None,
) -> np.ndarray:
    """
    Numba implementation using explicit loops (primitive operations).
    Performance is O(N) per spectrum, similar to the numpy cumsum trick,
    but fully compiled without objmode.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)

    if window < 1:
        raise ValueError("window must be >= 1")
    if window % 2 == 0:
        window += 1

    if intensity.dtype != np.float32:
        intensity = intensity.astype(np.float32)

    return _ma_loop_batch_jit(intensity, window, lengths)

@jit(nopython=True, fastmath=True, cache=True)
def _savgol_1d_core(signal: np.ndarray, kernels: np.ndarray) -> np.ndarray:
    """Core 1D Savitzky-Golay routine handling position-dependent convolution kernels."""
    n = signal.size
    window = kernels.shape[0]
    half = window // 2
    res = np.empty(n, dtype=np.float32)
    for i in range(n):
        if i < half:
            pos = i
        elif i >= n - half:
            pos = window - (n - i)
        else:
            pos = half
        
        k_weights = kernels[pos]
        start = i - pos
        val = 0.0
        for j in range(window):
            val += signal[start + j] * k_weights[j]
        res[i] = val
    return res

@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def savgol_batch_jit(data: np.ndarray, kernels: np.ndarray, lengths: np.ndarray = None) -> np.ndarray:
    """Batch Savitzky-Golay entry point. Kernels is a 2D array of position-dependent convolution kernels. Parallel over spectra."""
    if data.ndim == 1:
        return _savgol_1d_core(data, kernels)
    else:
        n_pixels, n_mz = data.shape
        res = np.empty((n_pixels, n_mz), dtype=np.float32)
        if lengths is None:
            for p in prange(n_pixels):
                res[p] = _savgol_1d_core(data[p], kernels)
        else:
            for p in prange(n_pixels):
                valid_len = lengths[p]
                # Slice the valid part, process it, and write back
                # Note: This writes valid result to res[p, :valid_len]
                # The rest of res[p] (padding area) remains uninitialized/garbage, which is fine as it will be trimmed later
                processed = _savgol_1d_core(data[p, :valid_len], kernels)
                res[p, :valid_len] = processed
        return res

@jit(nopython=True, fastmath=True, cache=True)
def _convolve1d_core_edge(signal: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Core 1D convolution routine with edge padding."""
    n = signal.size
    k_len = kernel.size
    half = k_len // 2
    res = np.empty(n, dtype=np.float32)

    for i in range(n):
        val = 0.0
        for j in range(k_len):
            idx = i - half + j
            if idx < 0:
                idx = 0
            elif idx >= n:
                idx = n - 1
            val += signal[idx] * kernel[j]
        res[i] = val
    return res

@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def convolve1d_batch_jit(data: np.ndarray, kernel: np.ndarray, lengths: np.ndarray = None) -> np.ndarray:
    """Batch convolution entry point. Parallel over spectra."""
    if data.ndim == 1:
        return _convolve1d_core_edge(data, kernel)
    else:
        n_pixels, n_mz = data.shape
        res = np.empty((n_pixels, n_mz), dtype=np.float32)
        if lengths is None:
            for p in prange(n_pixels):
                res[p] = _convolve1d_core_edge(data[p], kernel)
        else:
            for p in prange(n_pixels):
                valid_len = lengths[p]
                # Process only the valid length
                processed = _convolve1d_core_edge(data[p, :valid_len], kernel)
                res[p, :valid_len] = processed
        return res

def smooth_signal_ma_numba(intensity: np.ndarray, window: int = 5, lengths: np.ndarray = None, numba_max_threads: int = None) -> np.ndarray:
    """
    Numba-accelerated Moving Average smoothing.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)

    if window < 1:
        raise ValueError("window must be >= 1")
    
    # Ensure window is odd for centered convolution
    if window % 2 == 0:
        window += 1
        
    kernel = np.ones(window, dtype=np.float32) / window
    
    # intensity can be 1D or 2D
    # Ensure float32 for Numba
    if intensity.dtype != np.float32:
        intensity = intensity.astype(np.float32)
        
    return convolve1d_batch_jit(intensity, kernel, lengths)

def smooth_signal_gaussian_numba(intensity: np.ndarray, window: int = 5, sd: float = None, lengths: np.ndarray = None, numba_max_threads: int = None) -> np.ndarray:
    """
    Numba-accelerated Gaussian smoothing.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)

    if window < 1:
        raise ValueError("window must be >= 1")

    # Ensure window is odd
    if window % 2 == 0:
        window += 1

    if sd is None:
        sd = window / 4.0
        
    # Generate Gaussian kernel
    x = np.arange(-(window // 2), window // 2 + 1, dtype=np.float32)
    kernel = np.exp(-0.5 * (x / sd) ** 2).astype(np.float32)
    kernel /= kernel.sum()
    
    if intensity.dtype != np.float32:
        intensity = intensity.astype(np.float32)

    return convolve1d_batch_jit(intensity, kernel, lengths)

# --- Python wrapper layer (handles kernel pre‑computation) ---
def smooth_signal_savgol_numba(intensity: np.ndarray, window: int = 5, polyorder: int = 3, deriv: int = 0, delta: float = 1.0, lengths: np.ndarray = None, numba_max_threads: int = None) -> np.ndarray:
    """
    This function wraps the JIT-compiled Savitzky-Golay implementation and
    handles pre-computation of position-dependent convolution kernels.
    Parameters
    ----------
    intensity : numpy.ndarray
        Input intensity array. Either a 1D array of shape ``(n_points,)`` for
        a single spectrum or a 2D array of shape ``(n_spectra, n_points)``
        for a batch of spectra.
    window : int, optional
        Window length of the Savitzky-Golay filter. The actual window used is
        forced to be an odd integer and at least 3. Defaults to 5.
    polyorder : int, optional
        Order of the polynomial used to fit the samples. Must be strictly
        less than the effective window length. Defaults to 3.
    deriv : int, optional
        Order of the derivative to compute. ``0`` means only smoothing is
        performed. Defaults to 0.
    delta : float, optional
        Spacing of the samples along the x-axis. This is passed to
        :func:`scipy.signal.savgol_coeffs` when computing derivative
        coefficients. Defaults to 1.0.
    lengths : np.ndarray, optional
        Array of valid lengths for each spectrum in a 2D batch. If provided,
        smoothing is only applied up to the valid length for each spectrum.
    Returns
    -------
    numpy.ndarray
        Smoothed (or differentiated) intensity array of type ``float32`` with
        the same shape as ``intensity``.
    Raises
    ------
    ValueError
        If ``polyorder`` is greater than or equal to the effective window
        length.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)
        
    actual_window = max(3, window + (1 - window % 2))
    if polyorder >= actual_window:
        raise ValueError("polyorder must be < window")

    # Internal computation now uses float32 to reduce memory usage
    kernels = np.zeros((actual_window, actual_window), dtype=np.float32)
    for pos in range(actual_window):
        kernels[pos] = savgol_coeffs(actual_window, polyorder, deriv=deriv, delta=delta, pos=pos, use="dot").astype(np.float32)

    result32 = savgol_batch_jit(intensity, kernels, lengths)
    # Results are already float32; no further casting required
    return result32


@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def _ns_ma_kernel(neigh_intensity: np.ndarray) -> np.ndarray:
    """Numba kernel for NS moving-average smoothing (float64 compute).

    neigh_intensity has shape (n_points, k) and stores k-nearest-neighbour
    intensities per point. The loop is parallelized over points via prange.
    """
    n, k = neigh_intensity.shape
    out = np.empty(n, dtype=np.float64)
    inv_k = np.float64(1.0 / k)
    for i in prange(n):
        s = 0.0
        for j in range(k):
            s += neigh_intensity[i, j]
        out[i] = s * inv_k
    return out


@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def _ns_gaussian_kernel(neigh_intensity: np.ndarray, dists: np.ndarray, sd: float) -> np.ndarray:
    """Numba kernel for NS Gaussian smoothing over k-nearest neighbours (float64 compute)."""
    n, k = neigh_intensity.shape
    out = np.empty(n, dtype=np.float64)
    sd64 = np.float64(sd)
    for i in prange(n):
        row_sum = 0.0
        s = 0.0
        for j in range(k):
            exponent = -0.5 * (dists[i, j] / sd64) ** 2
            if exponent < -88.0:
                exponent = -88.0
            elif exponent > 0.0:
                exponent = 0.0
            w = np.exp(exponent)
            row_sum += w
            s += neigh_intensity[i, j] * w
        if row_sum == 0.0:
            out[i] = 0.0
        else:
            out[i] = s / row_sum
    return out


@jit(nopython=True, fastmath=True, cache=True, parallel=True)
def _ns_bilateral_kernel(neigh_intensity: np.ndarray, dists: np.ndarray, center_intensity: np.ndarray, sd_dist: float, sd_intensity: float) -> np.ndarray:
    """Numba kernel for NS bilateral filtering (float64 compute)."""
    n, k = neigh_intensity.shape
    out = np.empty(n, dtype=np.float64)
    sd_dist64 = np.float64(sd_dist)
    sd_int64 = np.float64(sd_intensity)
    for i in prange(n):
        ci = center_intensity[i]
        row_sum = 0.0
        s = 0.0
        for j in range(k):
            exponent_d = -0.5 * (dists[i, j] / sd_dist64) ** 2
            if exponent_d < -88.0:
                exponent_d = -88.0
            elif exponent_d > 0.0:
                exponent_d = 0.0
            wd = np.exp(exponent_d)
            diff = neigh_intensity[i, j] - ci
            exponent_i = -0.5 * (diff / sd_int64) ** 2
            if exponent_i < -88.0:
                exponent_i = -88.0
            elif exponent_i > 0.0:
                exponent_i = 0.0
            wi = np.exp(exponent_i)
            w = wd * wi
            row_sum += w
            s += neigh_intensity[i, j] * w
        if row_sum == 0.0:
            out[i] = 0.0
        else:
            out[i] = s / row_sum
    return out


def _ns_pre(intensity: np.ndarray, index: np.ndarray, k: int, p: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute NS neighbour structures for 1D spectra (float64 compute)."""
    intensity = np.asarray(intensity, dtype=np.float64)
    if index is None or index.shape[0] != intensity.shape[0]:
        index = np.arange(intensity.shape[0], dtype=np.float64)
    else:
        index = np.asarray(index, dtype=np.float64)
    if k < 1 or p < 1:
        raise ValueError("k and p must be positive integers")
    if k > intensity.shape[0]:
        k = intensity.shape[0]
    tree = cKDTree(index.reshape(-1, 1))
    dists, idxs = tree.query(index.reshape(-1, 1), k=k, p=p)
    if dists.ndim == 1:
        dists = dists.reshape(-1, 1)
        idxs = idxs.reshape(-1, 1)
    nan_idx = np.where(np.isnan(intensity))[0]
    if nan_idx.size > 0:
        neigh = intensity[idxs[nan_idx]]
        fill_vals = np.nanmean(neigh, axis=1)
        global_median = np.nanmedian(intensity)
        if np.isnan(global_median):
            global_median = 0.0
        mask_nan = np.isnan(fill_vals)
        if mask_nan.any():
            fill_vals[mask_nan] = global_median
        intensity[nan_idx] = fill_vals
    neigh_intensity = intensity[idxs]
    return neigh_intensity, dists, idxs


def smooth_ns_signal_ma_numba(intensity: np.ndarray, index: np.ndarray, k: int = 5, p: int = 1, numba_max_threads: int = None) -> np.ndarray:
    """Neighbourhood-smoothing MA (Numba) entry point.

    Given a 1D spectrum and its index axis, performs k-NN search via _ns_pre
    and applies the parallel moving-average kernel.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)
    neigh_intensity, _, _ = _ns_pre(intensity, index, k=k, p=p)
    result64 = _ns_ma_kernel(neigh_intensity)
    return result64.astype(np.float32)


def smooth_ns_signal_gaussian_numba(intensity: np.ndarray, index: np.ndarray, k: int = 5, p: int = 1, sd: float | None = None, numba_max_threads: int = None) -> np.ndarray:
    """Neighbourhood-smoothing Gaussian (Numba) entry point.

    Estimates a default spatial scale from neighbour distances when sd is
    None, then applies the parallel Gaussian kernel.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)
    neigh_intensity, dists, _ = _ns_pre(intensity, index, k=k, p=p)
    dists_max = np.max(dists, axis=1)
    sd_val = np.median(dists_max) / 2.0 if sd is None else sd
    if sd_val <= 0:
        raise ValueError("sd must be positive")
    result64 = _ns_gaussian_kernel(neigh_intensity, dists, float(sd_val))
    return result64.astype(np.float32)


def smooth_ns_signal_bi_numba(intensity: np.ndarray, index: np.ndarray, k: int = 5, p: int = 2, sd_dist: float | None = None, sd_intensity: float | None = None, numba_max_threads: int = None) -> np.ndarray:
    """Neighbourhood-smoothing bilateral (Numba) entry point.

    Uses spatial distances to derive sd_dist (if not provided) and robust
    MAD-based estimate for sd_intensity, then applies the parallel bilateral
    kernel.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)
    neigh_intensity, dists, _ = _ns_pre(intensity, index, k=k, p=p)
    dists_max = np.max(dists, axis=1)
    sd_dist_val = np.median(dists_max) / 2.0 if sd_dist is None else sd_dist
    if sd_dist_val <= 0:
        raise ValueError("sd_dist must be positive")
    center = np.asarray(intensity, dtype=np.float64)
    sd_int_val = stats.median_abs_deviation(center, nan_policy="omit", scale="normal") if sd_intensity is None else sd_intensity
    if sd_int_val <= 0:
        raise ValueError("sd_intensity must be positive")
    result64 = _ns_bilateral_kernel(neigh_intensity, dists, center, float(sd_dist_val), float(sd_int_val))
    return result64.astype(np.float32)