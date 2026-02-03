import numpy as np
from numba import jit, prange, set_num_threads
@jit(nopython=True, parallel=True, cache=True)
def _gaussian_approx_loop_batch_jit(intensity: np.ndarray, window: int, lengths: np.ndarray = None) -> np.ndarray:
    """Batch entry point for gaussian_approx_loop."""
    if intensity.ndim == 1:
        # Apply MA 3 times
        temp = _ma_loop_core_edge(intensity, window)
        temp = _ma_loop_core_edge(temp, window)
        return _ma_loop_core_edge(temp, window)
    else:
        n_pixels, n_mz = intensity.shape
        res = np.empty((n_pixels, n_mz), dtype=np.float32)
        
        # Explicit parallel loop over spectra
        if lengths is None:
            for p in prange(n_pixels):
                # Apply MA 3 times per spectrum
                temp = _ma_loop_core_edge(intensity[p], window)
                temp = _ma_loop_core_edge(temp, window)
                res[p] = _ma_loop_core_edge(temp, window)
        else:
            for p in prange(n_pixels):
                valid_len = lengths[p]
                # Process only valid part
                # Note: We must be careful about temp buffers.
                # Since _ma_loop_core_edge allocates new array, this is safe but memory intensive (3x alloc).
                # For optimization, we could use in-place or pre-alloc, but let's keep it simple first.
                
                # Step 1
                temp1 = _ma_loop_core_edge(intensity[p, :valid_len], window)
                # Step 2
                temp2 = _ma_loop_core_edge(temp1, window)
                # Step 3
                processed = _ma_loop_core_edge(temp2, window)
                
                res[p, :valid_len] = processed
                
        return res

def smooth_signal_gaussian_approx_loop(intensity: np.ndarray, window: int = 5, lengths: np.ndarray = None, numba_max_threads: int = None) -> np.ndarray:
    """
    Approximates Gaussian smoothing by applying Moving Average 3 times.
    Performance is O(N) instead of O(N*Window), but the shape is an approximation.
    """
    if numba_max_threads is not None:
        set_num_threads(numba_max_threads)
        
    if window < 1:
        raise ValueError("window must be >= 1")
    if window % 2 == 0:
        window += 1
        
    if intensity.dtype != np.float32:
        intensity = intensity.astype(np.float32)
        
    return _gaussian_approx_loop_batch_jit(intensity, window, lengths)

# --- 4. Pure Numba Implementation with Explicit Loop (O(N) Sliding Window) ---

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
            val = signal[n-1]
        else:
            val = signal[idx]
        current_sum += val
        
    res[0] = current_sum / window
    
    # Slide the window
    for i in range(1, n):
        # Window moves from centered at i-1 to centered at i
        # Leaving index: (i-1) - radius
        # Entering index: i + radius
        
        leaving_idx = i - 1 - radius
        entering_idx = i + radius
        
        # Get leaving value
        if leaving_idx < 0:
            leaving_val = signal[0]
        else:
            leaving_val = signal[leaving_idx]
            
        # Get entering value
        if entering_idx >= n:
            entering_val = signal[n-1]
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
    else:
        n_pixels, n_mz = intensity.shape
        res = np.empty((n_pixels, n_mz), dtype=np.float32)
        
        # Explicit parallel loop over spectra
        if lengths is None:
            for p in prange(n_pixels):
                res[p] = _ma_loop_core_edge(intensity[p], window)
        else:
            for p in prange(n_pixels):
                valid_len = lengths[p]
                # Process only valid part
                processed = _ma_loop_core_edge(intensity[p, :valid_len], window)
                res[p, :valid_len] = processed
                
        return res

def smooth_signal_ma_loop(intensity: np.ndarray, window: int = 5, lengths: np.ndarray = None, numba_max_threads: int = None) -> np.ndarray:
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


