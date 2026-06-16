"""Helpers for per-spectrum isotope filtering."""

from typing import Literal

import numpy as np

from massflow.preprocess.numba.isotope_filter_numba import (
    apply_isotope_keep_mask_flat_parallel_jit,
    isotope_keep_masks_flat_parallel_jit,
)
from massflow.tools.funs import infer_shared_mz, prepare_flat_inputs

ToleranceUnit = Literal["da", "ppm"]


def _input_validation(
    mz_arr: np.ndarray,
    lengths_arr: np.ndarray,
    mz_tolerance: float,
    units: ToleranceUnit,
    min_charge: int,
    max_charge: int,
    min_isotope_count: int,
    max_isotope_count: int,
    start_intensity_check: int,
) -> None:
    """Validate inputs for isotope filtering."""
    if infer_shared_mz(mz_arr, lengths_arr):
        raise ValueError("Isotope filtering requires processed + centroided m/z data.")
    if mz_tolerance <= 0 or not np.isfinite(mz_tolerance):
        raise ValueError("mz_tolerance must be a finite positive number.")
    if units not in {"da", "ppm"}:
        raise ValueError("units must be either 'da' or 'ppm'.")
    if (units == "ppm" and mz_tolerance > 100) or (
        units == "da" and mz_tolerance > 0.1
    ):
        raise ValueError("mz_tolerance must not exceed 100 ppm or 0.1 Da.")
    if min_charge < 1 or max_charge < min_charge:
        raise ValueError("Charges must satisfy 1 <= min_charge <= max_charge.")
    if (
        min_isotope_count < 2
        or max_isotope_count < min_isotope_count
    ):
        raise ValueError(
            "Isotope counts must satisfy 2 <= min_isotope_count <= max_isotope_count."
        )
    if start_intensity_check < 0:
        raise ValueError("start_intensity_check must be non-negative.")


def filter_peaks_by_isotopes(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    *,
    mz_tolerance: float = 10.0,
    units: Literal["da", "ppm"] = "ppm",
    min_charge: int = 1,
    max_charge: int = 1,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Remove isotope peaks from each centroided spectrum.

    Isotope clusters are identified using the simple OpenMS deisotoping model:
    peaks are expected at successive C13-C12 mass differences divided by the
    charge state. The lowest-m/z peak of each cluster and all unassigned peaks
    are retained.
    """
    mz_arr, intensity_arr, lengths_arr = prepare_flat_inputs(
        mz_data,
        intensity,
        lengths,
    )
    _input_validation(
        mz_arr,
        lengths_arr,
        mz_tolerance,
        units,
        min_charge,
        max_charge,
        min_isotope_count,
        max_isotope_count,
        start_intensity_check
    )

    keep, filtered_lengths = isotope_keep_masks_flat_parallel_jit(
        mz_arr,
        intensity_arr,
        lengths_arr,
        mz_tolerance,
        units == "ppm",
        min_charge,
        max_charge,
        min_isotope_count,
        max_isotope_count,
        use_decreasing_model,
        start_intensity_check,
    )
    filtered_mz, filtered_intensity = apply_isotope_keep_mask_flat_parallel_jit(
        mz_arr,
        intensity_arr,
        lengths_arr,
        keep,
        filtered_lengths,
    )
    return filtered_mz, filtered_intensity, filtered_lengths
