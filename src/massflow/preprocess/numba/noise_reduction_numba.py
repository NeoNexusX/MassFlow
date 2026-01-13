import numpy as np
from numba import jit, prange, set_num_threads
from scipy.signal import savgol_coeffs
from massflow.tools.logger import get_logger

logger = get_logger("noise_reduction_numba")
set_num_threads(2)

@jit(nopython=True, fastmath=True, cache=True)
def _savgol_1d_core(signal: np.ndarray, kernels: np.ndarray) -> np.ndarray:
    """Core 1D Savitzky-Golay routine handling position-dependent convolution kernels."""
    n = signal.size
    window = kernels.shape[0]
    half = window // 2
    res = np.empty(n, dtype=np.float64)
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
def savgol_batch_jit(data: np.ndarray, kernels: np.ndarray) -> np.ndarray:
    """Batch Savitzky-Golay entry point. Kernels is a 2D array of position-dependent convolution kernels. Parallel over spectra."""
    if data.ndim == 1:
        return _savgol_1d_core(data, kernels)
    else:
        n_pixels, n_mz = data.shape
        res = np.empty((n_pixels, n_mz), dtype=np.float64)
        for p in prange(n_pixels):
            res[p] = _savgol_1d_core(data[p], kernels)
        return res

# --- Python wrapper layer (handles kernel pre‑computation) ---
def smooth_signal_savgol_numba(intensity: np.ndarray, window: int = 5, polyorder: int = 3, deriv: int = 0, delta: float = 1.0) -> np.ndarray:
    actual_window = max(3, window + (1 - window % 2))
    if polyorder >= actual_window:
        raise ValueError("polyorder must be < window")
    
    kernels = np.zeros((actual_window, actual_window), dtype=np.float64)
    for pos in range(actual_window):
        kernels[pos] = savgol_coeffs(actual_window, polyorder, deriv=deriv, delta=delta, pos=pos, use="dot")
    
    return savgol_batch_jit(intensity.astype(np.float64), kernels)