from typing import Optional
import numpy as np
from numba import jit, prange

# -----------------------------------------------------------------------------
# Core Calculation Kernels (Helpers)
# -----------------------------------------------------------------------------


@jit(nopython=True, cache=True)
def _get_tic(arr: np.ndarray) -> float:
    return np.sum(arr)


@jit(nopython=True, cache=True)
def _get_rms(arr: np.ndarray) -> float:
    # RMS = sqrt( sum(x^2) / n )
    s = 0.0
    n = arr.shape[0]
    if n == 0:
        return 0.0
    for i in range(n):
        val = arr[i]
        s += val * val
    return np.sqrt(s / n)


@jit(nopython=True, cache=True)
def _get_median(arr: np.ndarray) -> float:
    return np.median(arr)  # type: ignore


@jit(nopython=True, cache=True)
def _apply_unit_scaling_inplace(arr: np.ndarray):
    """Applies Min-Max scaling to [0, 1] in-place."""
    n = arr.shape[0]
    if n == 0:
        return

    min_val = np.inf
    max_val = -np.inf

    # Pass 1: Find Min/Max
    for i in range(n):
        val = arr[i]
        if val < min_val:
            min_val = val
        if val > max_val:
            max_val = val

    rng = max_val - min_val

    # Pass 2: Scale
    if rng > 0:
        for i in range(n):
            arr[i] = (arr[i] - min_val) / rng
    else:
        for i in range(n):
            arr[i] = 0.0


# -----------------------------------------------------------------------------
# Batch Logic
# -----------------------------------------------------------------------------


@jit(nopython=True, cache=True)
def _process_row(
    in_arr: np.ndarray,
    out_arr: np.ndarray,
    method_enum: int,
    scale_factor: float,
    do_unit_scale: bool,
):
    """
    Process a single row (view) and write to output view.
    """
    # 1. Calculate Base Factor
    base = 0.0
    if method_enum == 0:  # TIC
        base = _get_tic(in_arr)
    elif method_enum == 1:  # RMS
        base = _get_rms(in_arr)
    elif method_enum == 2:  # MEDIAN
        base = _get_median(in_arr)

    # Safety Check
    if base <= 0.0 or np.isnan(base):
        # Assuming out_arr is zeroed by caller, just return
        return

    # 2. Apply Normalization
    factor = scale_factor / base
    n = in_arr.shape[0]

    for i in range(n):
        out_arr[i] = in_arr[i] * factor

    # 3. Apply Unit Scaling if requested
    if do_unit_scale:
        _apply_unit_scaling_inplace(out_arr)


@jit(nopython=True, parallel=True, cache=True)
def _normalize_batch_jit(
    intensity: np.ndarray,
    method_enum: int,
    scale_factor: float,
    do_unit_scale: bool,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:

    n_pixels = intensity.shape[0]
    n_mz = intensity.shape[1]

    # Initialize output (zeros handles padding implicitly)
    res = np.zeros((n_pixels, n_mz), dtype=intensity.dtype)

    # Parallel Loop
    if lengths is None:
        for i in prange(n_pixels):
            _process_row(intensity[i], res[i], method_enum, scale_factor, do_unit_scale)
    else:
        for i in prange(n_pixels):
            valid_len = lengths[i]
            if valid_len > 0:
                # Use slicing (views) for efficiency
                _process_row(
                    intensity[i, :valid_len],
                    res[i, :valid_len],
                    method_enum,
                    scale_factor,
                    do_unit_scale,
                )
    return res


# -----------------------------------------------------------------------------
# Main Dispatcher
# -----------------------------------------------------------------------------


def normalizer_numba(
    intensity: np.ndarray,
    method: str = "tic",
    scale: float = 1.0,
    scale_method: str = "none",
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Numba backend for normalization with 2D batch support.

    Parameters:
        intensity: 2D array (n_pixels, n_mz) or 1D array.
        method: 'tic', 'rms', 'median'.
        scale: Amplitude scaling factor.
        scale_method: 'none' or 'unit'.
        lengths: (Optional) 1D int array of valid lengths per spectrum.
    """

    # Ensure float, but prefer float32 for speed if input is already float32
    if intensity.dtype not in [np.float32, np.float64]:
        intensity = intensity.astype(np.float64)

    # Handle 1D input by treating as batch of 1
    is_1d = intensity.ndim == 1
    if is_1d:
        intensity = intensity[np.newaxis, :]
        if lengths is not None and np.isscalar(lengths):
            lengths = np.array([lengths], dtype=np.int64)

    # Map Method
    method_clean = method.strip().lower()
    if method_clean == "tic":
        m_enum = 0
    elif method_clean == "rms":
        m_enum = 1
    elif method_clean == "median":
        m_enum = 2
    else:
        raise ValueError(f"Unknown method for numba normalizer: {method}")

    # Map Scale Method
    scale_clean = (scale_method or "none").strip().lower()
    do_unit = scale_clean == "unit"

    # Execute Parallel Batch
    result = _normalize_batch_jit(intensity, m_enum, float(scale), do_unit, lengths)

    if is_1d:
        return result[0]
    return result
