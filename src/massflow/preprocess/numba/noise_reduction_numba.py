from typing import Optional
import numpy as np
from numba import njit, prange
from scipy.signal import savgol_coeffs
from scipy.spatial import cKDTree  # type: ignore
from scipy import stats
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
    current_sum = left_edge * (radius + 1)
    for j in range(1, radius + 1):
        current_sum += signal[j] if j < valid_len else right_edge

    out[0] = current_sum * inv_window

    # Slide the window
    for i in range(1, valid_len):
        leaving_idx = i - 1 - radius
        entering_idx = i + radius

        leaving_val = left_edge if leaving_idx < 0 else signal[leaving_idx]
        entering_val = right_edge if entering_idx >= valid_len else signal[entering_idx]

        current_sum += entering_val - leaving_val
        out[i] = current_sum * inv_window


@njit(cache=True)
def _lengths_to_offsets(lengths: np.ndarray) -> np.ndarray:
    n = lengths.size
    offsets = np.empty(n + 1, dtype=np.int64)
    offsets[0] = 0
    for i in range(n):
        offsets[i + 1] = offsets[i] + lengths[i]
    return offsets


@njit(parallel=True, cache=True, fastmath=False)
def _ma_flat_jit(flat: np.ndarray, window: int, lengths: np.ndarray) -> np.ndarray:
    """Flat batch entry point for ma_numba."""
    res = np.empty(flat.size, dtype=np.float32)
    offsets = _lengths_to_offsets(lengths)

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


@njit(fastmath=True, cache=True, parallel=True)
def savgol_flat_jit(flat: np.ndarray, kernels: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Flat Savitzky-Golay entry point. Parallel over spectra segments."""
    res = np.empty(flat.size, dtype=np.float32)
    offsets = _lengths_to_offsets(lengths)
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


@njit(fastmath=True, cache=True, parallel=True)
def convolve1d_flat_jit(flat: np.ndarray, kernel: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Flat convolution entry point. Parallel over spectra segments."""
    res = np.empty(flat.size, dtype=np.float32)
    offsets = _lengths_to_offsets(lengths)
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



# ==================== Non-JIT wrappers and preprocessing ====================


def _prepare_intensity_inputs(
    intensity: np.ndarray,
    lengths: Optional[np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize flat input and lengths for flat-mode kernels."""
    if intensity is None or not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        raise ValueError("intensity must be a 1D numpy array")

    intensity_arr = intensity.astype(np.float32, copy=False)

    if lengths is None:
        lengths_arr = np.array([intensity_arr.size], dtype=np.int64)
    else:
        lengths_arr = np.asarray(lengths, dtype=np.int64)
        if lengths_arr.ndim != 1:
            raise ValueError("lengths must be a 1D array")
        if int(np.sum(lengths_arr)) != intensity_arr.size:
            raise ValueError("sum(lengths) must equal intensity.size")

    return intensity_arr, lengths_arr


def smooth_signal_ma_numba(
    intensity: np.ndarray,
    window: int = 5,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Flat-mode moving-average smoothing."""
    if window < 1:
        raise ValueError("window must be >= 1")

    window = window + 1 if window % 2 == 0 else window
    intensity_arr, lengths_arr = _prepare_intensity_inputs(intensity, lengths)

    return _ma_flat_jit(intensity_arr, window, lengths_arr)


# --- Python wrapper layer (handles kernel pre‑computation) ---
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

    intensity_arr, lengths_arr = _prepare_intensity_inputs(intensity, lengths)

    kernels = np.zeros((actual_window, actual_window), dtype=np.float32)
    for pos in range(actual_window):
        kernels[pos] = savgol_coeffs(
            actual_window, polyorder, deriv=deriv, delta=delta, pos=pos, use="dot"
        ).astype(np.float32)

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

    intensity_arr, lengths_arr = _prepare_intensity_inputs(intensity, lengths)

    # Generate Gaussian kernel
    x = np.arange(-(window // 2), window // 2 + 1, dtype=np.float32)
    kernel = np.exp(-0.5 * (x / sd) ** 2).astype(np.float32)
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
    neigh_intensity, dists, _ = ns_signal_pre(intensity, index, k=window, p=p)
    dists_max = np.max(dists, axis=1)
    sd_val = np.median(dists_max) / 2.0 if sd is None else sd
    if sd_val <= 0:
        raise ValueError("sd must be positive")
    result64 = _ns_gaussian_kernel(neigh_intensity, dists, float(sd_val))
    return result64.astype(np.float32)


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
    return result64.astype(np.float32)


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
