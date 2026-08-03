"""Numba kernels for filtering continuous flat spectra."""

import numpy as np
from numba import njit, prange

from massflow.tools.funs import lengths_to_offsets


# ---------------------------------------------------------------------------
# Aggregate peak count and intensity across spectra
# ---------------------------------------------------------------------------

@njit(parallel=True, cache=True)
def aggregate_peak_profile_flat_parallel_jit(
    intensity_flat: np.ndarray,
    lengths: np.ndarray,
    spectrum_size: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Compute occurrence counts and summed intensity for one flat batch."""
    counts = np.zeros(spectrum_size, dtype=np.int64)
    summed = np.zeros(spectrum_size, dtype=np.float64)
    offsets = lengths_to_offsets(lengths)

    for peak_index in prange(spectrum_size):  # pylint: disable=not-an-iterable
        count = 0
        value = 0.0
        for spectrum_index in range(lengths.size):
            if peak_index < lengths[spectrum_index]:
                intensity = intensity_flat[offsets[spectrum_index] + peak_index]
                value += intensity
                if intensity > 0:
                    count += 1
        counts[peak_index] = count
        summed[peak_index] = value

    return counts, summed, int(lengths.size)


@njit(parallel=True, cache=True)
def count_peaks_flat_parallel_jit(
    intensity_flat: np.ndarray,
    lengths: np.ndarray,
    spectrum_size: int,
) -> tuple[np.ndarray, int]:
    """Count spectra with positive intensity at each shared m/z position."""
    counts = np.zeros(spectrum_size, dtype=np.int64)
    offsets = lengths_to_offsets(lengths)

    for peak_index in prange(spectrum_size):  # pylint: disable=not-an-iterable
        count = 0
        for spectrum_index in range(lengths.size):
            if (
                peak_index < lengths[spectrum_index]
                and intensity_flat[offsets[spectrum_index] + peak_index] > 0
            ):
                count += 1
        counts[peak_index] = count

    return counts, int(lengths.size)


# ---------------------------------------------------------------------------
# Aggregate intensity spectrum for peak filtering
# ---------------------------------------------------------------------------

@njit(parallel=True, cache=True)
def sum_intensity_flat_parallel_jit(
    intensity_flat: np.ndarray,
    lengths: np.ndarray,
    spectrum_size: int,
) -> tuple[np.ndarray, int]:
    """Sum a flat batch along its shared m/z axis."""
    summed = np.zeros(spectrum_size, dtype=np.float64)
    offsets = lengths_to_offsets(lengths)

    for peak_index in prange(spectrum_size):  # pylint: disable=not-an-iterable
        value = 0.0
        for spectrum_index in range(lengths.size):
            if peak_index < lengths[spectrum_index]:
                value += intensity_flat[offsets[spectrum_index] + peak_index]
        summed[peak_index] = value

    return summed, int(lengths.size)


# ---------------------------------------------------------------------------
# Apply a shared peak-selection mask to every spectrum
# ---------------------------------------------------------------------------

@njit(parallel=True, cache=True)
def filter_peaks_flat_parallel_jit(
    intensity_flat: np.ndarray,
    lengths: np.ndarray,
    keep_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Filter every spectrum in a flat batch with one shared peak mask."""
    offsets = lengths_to_offsets(lengths)
    filtered_lengths = np.empty(lengths.size, dtype=np.uint32)

    for spectrum_index in prange(lengths.size):  # pylint: disable=not-an-iterable
        valid_length = int(lengths[spectrum_index])
        kept = 0
        for peak_index in range(valid_length):
            if keep_mask[peak_index]:
                kept += 1
        filtered_lengths[spectrum_index] = kept

    filtered_offsets = lengths_to_offsets(filtered_lengths)
    filtered = np.empty(filtered_offsets[-1], dtype=intensity_flat.dtype)

    for spectrum_index in prange(lengths.size):  # pylint: disable=not-an-iterable
        input_start = offsets[spectrum_index]
        output_index = filtered_offsets[spectrum_index]
        for peak_index in range(lengths[spectrum_index]):
            if keep_mask[peak_index]:
                filtered[output_index] = intensity_flat[input_start + peak_index]
                output_index += 1

    return filtered, filtered_lengths
