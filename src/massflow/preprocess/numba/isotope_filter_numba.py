"""Numba kernels for per-spectrum isotope filtering."""

import numpy as np
from numba import njit, prange

from massflow.tools.funs import lengths_to_offsets

ISOTOPE_MASS_DIFFERENCE = 1.003354838


@njit(cache=True)
def _find_nearest_peak_jit(
    mz: np.ndarray,
    target: float,
    tolerance: float,
) -> int | np.intp:
    """Return the nearest peak index within tolerance, or -1."""
    left = np.searchsorted(mz, target - tolerance)
    right = np.searchsorted(mz, target + tolerance, side="right")
    if left >= right:
        return -1

    nearest = left
    nearest_error = abs(mz[left] - target)
    for index in range(left + 1, right):
        error = abs(mz[index] - target)
        if error < nearest_error:
            nearest = index
            nearest_error = error
    return nearest


@njit(cache=True)
def _build_isotope_keep_mask_jit(
    mz: np.ndarray,
    intensity: np.ndarray,
    mz_tolerance: float,
    use_ppm: bool,
    min_charge: int,
    max_charge: int,
    min_isotope_count: int,
    max_isotope_count: int,
    use_decreasing_model: bool,
    start_intensity_check: int,
) -> np.ndarray:
    """Build a keep mask using the OpenMS simple deisotoping model."""
    keep = np.ones(mz.size, dtype=np.bool_)
    assigned = np.zeros(mz.size, dtype=np.bool_)
    extensions = np.empty(max_isotope_count, dtype=np.int64)

    for current_peak in range(mz.size):
        if assigned[current_peak]:
            continue

        current_mz = mz[current_peak]
        tolerance = (
            current_mz * mz_tolerance * 1e-6
            if use_ppm
            else mz_tolerance
        )

        for charge in range(max_charge, min_charge - 1, -1):
            extension_count = 1
            extensions[0] = current_peak

            for isotope_index in range(1, max_isotope_count):
                expected_mz = (
                    current_mz
                    + isotope_index * ISOTOPE_MASS_DIFFERENCE / charge
                )
                candidate = _find_nearest_peak_jit(mz, expected_mz, tolerance)
                if candidate < 0 or assigned[candidate]:
                    break
                if (
                    use_decreasing_model
                    and isotope_index >= start_intensity_check
                    and intensity[candidate] > intensity[extensions[extension_count - 1]]
                ):
                    break

                extensions[extension_count] = candidate
                extension_count += 1

            if extension_count >= min_isotope_count:
                assigned[current_peak] = True
                for extension_index in range(1, extension_count):
                    isotope_peak = extensions[extension_index]
                    assigned[isotope_peak] = True
                    keep[isotope_peak] = False
                break

    return keep


@njit(parallel=True, cache=True)
def isotope_keep_masks_flat_parallel_jit(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    mz_tolerance: float,
    use_ppm: bool,
    min_charge: int,
    max_charge: int,
    min_isotope_count: int,
    max_isotope_count: int,
    use_decreasing_model: bool,
    start_intensity_check: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build per-spectrum isotope keep masks and output lengths."""
    offsets = lengths_to_offsets(lengths)
    keep = np.empty(intensity.size, dtype=np.bool_)
    filtered_lengths = np.empty(lengths.size, dtype=np.uint32)

    for spectrum_index in prange(lengths.size):  # pylint: disable=not-an-iterable
        start = offsets[spectrum_index]
        end = offsets[spectrum_index + 1]
        spectrum_keep = _build_isotope_keep_mask_jit(
            mz_data[start:end],
            intensity[start:end],
            mz_tolerance,
            use_ppm,
            min_charge,
            max_charge,
            min_isotope_count,
            max_isotope_count,
            use_decreasing_model,
            start_intensity_check,
        )
        kept = 0
        for local_index in range(spectrum_keep.size):
            keep[start + local_index] = spectrum_keep[local_index]
            if spectrum_keep[local_index]:
                kept += 1
        filtered_lengths[spectrum_index] = kept

    return keep, filtered_lengths


@njit(parallel=True, cache=True)
def apply_isotope_keep_mask_flat_parallel_jit(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    keep: np.ndarray,
    filtered_lengths: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compact flat m/z and intensity arrays using per-spectrum keep masks."""
    input_offsets = lengths_to_offsets(lengths)
    output_offsets = lengths_to_offsets(filtered_lengths)
    filtered_mz = np.empty(output_offsets[-1], dtype=mz_data.dtype)
    filtered_intensity = np.empty(output_offsets[-1], dtype=intensity.dtype)

    for spectrum_index in prange(lengths.size):  # pylint: disable=not-an-iterable
        output_index = output_offsets[spectrum_index]
        for input_index in range(
            input_offsets[spectrum_index],
            input_offsets[spectrum_index + 1],
        ):
            if keep[input_index]:
                filtered_mz[output_index] = mz_data[input_index]
                filtered_intensity[output_index] = intensity[input_index]
                output_index += 1

    return filtered_mz, filtered_intensity
