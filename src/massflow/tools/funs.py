import os
from inspect import signature
from typing import Any, Callable, Optional
from numba import njit
import numpy as np

def is_valid_file(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

def dispatch_with_supported_kwargs(func: Callable[..., Any], **kwargs: Any) -> np.ndarray:
    """Dispatch to function with only supported kwargs."""
    supported = signature(func).parameters
    filtered_kwargs = {name: value for name, value in kwargs.items() if name in supported}
    result = func(**filtered_kwargs)
    return np.asarray(result)

def prepare_flat_inputs(
    mz_data: Optional[np.ndarray],
    intensity: np.ndarray,
    lengths: Optional[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Validate 1D flat input and optional lengths."""
    if intensity is None or not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        raise ValueError("intensity must be a 1D numpy array")

    if mz_data is not None and mz_data.ndim != 1:
        raise ValueError("mz_data must be a 1D array if provided")

    mz_arr = np.asarray(mz_data, dtype=np.float64) if mz_data is not None else np.array([], dtype=intensity.dtype)
    intensity_arr = intensity

    if lengths is None:
        lengths_arr = np.array([intensity_arr.size], dtype=np.uint32)
    else:
        if lengths.ndim != 1:
            raise ValueError("lengths must be a 1D array")
        if np.any(lengths < 0):
            raise ValueError("lengths must contain non-negative integers")
        if int(np.sum(lengths)) != intensity_arr.size:
            raise ValueError("sum(lengths) must equal intensity.size")
        lengths_arr = np.asarray(lengths, dtype=np.uint32)

    return mz_arr, intensity_arr, lengths_arr

def infer_shared_mz(
    mz_data: np.ndarray,
    lengths: np.ndarray
) -> bool:
    """Infer if mz_data is shared across spectra based on its length and the lengths array."""
    total_points = int(np.sum(lengths, dtype=np.int64))
    max_len = int(np.max(lengths)) if lengths.size > 0 else 0
    is_shared_mz = mz_data.size != total_points

    if is_shared_mz and mz_data.size < max_len:
        raise ValueError(
            "Shared m/z axis is shorter than at least one spectrum length: "
            f"mz_size={mz_data.size}, max_len={max_len}."
        )

    return is_shared_mz

@njit(cache=True)
def lengths_to_offsets(lengths: np.ndarray) -> np.ndarray:
    n = lengths.size
    offsets = np.empty(n + 1, dtype=np.int64)
    offsets[0] = 0
    for i in range(n):
        offsets[i + 1] = offsets[i] + lengths[i]
    return offsets
