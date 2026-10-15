import os
from typing import Optional
from numba import njit
import numpy as np

def is_valid_file(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

def is_zero_based(coordinates: np.ndarray) -> bool:
    """Return whether pixel coordinates are 0-based.

    imzML coordinates are 1-based, but some writers emit 0-based ones; a 0 in
    either the x or the y column marks the whole set as 0-based.
    """
    xy = np.asarray(coordinates)[:, :2]
    return bool(xy.size) and (int(xy[:, 0].min()) == 0 or int(xy[:, 1].min()) == 0)

def zero_based_xy(coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the x and y columns of pixel coordinates as 0-based int64 arrays."""
    xy = np.asarray(coordinates)[:, :2].astype(np.int64, copy=False)
    if not is_zero_based(xy):
        xy = xy - 1
    return xy[:, 0], xy[:, 1]

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

    mz_arr = mz_data if mz_data is not None else np.array([], dtype=intensity.dtype)
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
