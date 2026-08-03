from __future__ import annotations

from typing import Any, Literal

import numpy as np

from massflow.tools.logger import get_logger

logger = get_logger("massflow.tools.infer_spectrum_type")

SpectrumType = Literal["profile", "centroid"]


def resolve_spectrum_type(data_manager: Any) -> SpectrumType:
    """Resolve spectrum type from canonical metadata flags."""
    if data_manager.ms.meta is None:
        raise ValueError("spectrum metadata is missing. Please load metadata before resolving spectrum type.")

    profile = getattr(data_manager.ms.meta, "profile_spectrum", None)
    centroid = getattr(data_manager.ms.meta, "centroid_spectrum", None)

    if profile is True and centroid is None:
        return "profile"

    if profile is None and centroid is True:
        return "centroid"

    if (profile is None and centroid is None) or (profile is True and centroid is True):
        logger.info("Spectrum type is all True or all None. Falling back to inference.")
        inferred_spectrum_type = infer_spectrum_type(data_manager)
        logger.info(f"Inferred spectrum type: {inferred_spectrum_type!r}. Setting metadata accordingly.")
        return inferred_spectrum_type

    raise ValueError(
        f"invalid spectrum type metadata in imzML file: profile_spectrum={profile!r}, centroid_spectrum={centroid!r}. "
        "Expected exactly one flag to be True and the other to be None."
    )

def infer_spectrum_type(data_manager: Any, *, sample_size: int = 8) -> SpectrumType:
    """Infer spectrum type directly from spectra stored in a data manager."""
    if data_manager is None or not hasattr(data_manager, "ms"):
        raise ValueError("data_manager with attribute 'ms' is required.")

    ms = data_manager.ms
    if not hasattr(ms, "__len__") or not hasattr(ms, "__getitem__"):
        raise ValueError("data_manager.ms must support len() and index access.")

    spectrum_count = len(ms)
    if spectrum_count == 0:
        raise ValueError("Cannot infer spectrum type from an empty data manager.")

    sample_count = min(max(int(sample_size), 1), spectrum_count)
    sample_indices = np.unique(np.linspace(0, spectrum_count - 1, sample_count, dtype=int))
    scores: list[float] = []

    for index in sample_indices:
        spectrum = ms[int(index)]
        mz_data = getattr(spectrum, "mz_list", None)
        intensity_data = getattr(spectrum, "intensity", None)

        if mz_data is None or intensity_data is None:
            continue

        mz_arr = np.asarray(mz_data, dtype=np.float64)
        intensity_arr = np.abs(np.asarray(intensity_data, dtype=np.float64))

        if mz_arr.ndim != 1 or intensity_arr.ndim != 1 or mz_arr.size != intensity_arr.size:
            continue

        valid_mask = np.isfinite(mz_arr) & np.isfinite(intensity_arr)
        mz_arr = mz_arr[valid_mask]
        intensity_arr = intensity_arr[valid_mask]

        if mz_arr.size == 0:
            continue

        if mz_arr.size < 8:
            scores.append(-2.0)
            continue

        amplitude = float(np.max(intensity_arr))
        if amplitude <= 0.0:
            continue

        significant_mask = intensity_arr > max(amplitude * 1e-6, np.finfo(np.float64).eps) # pylint: disable=E1101
        significant_count = int(np.count_nonzero(significant_mask))
        if significant_count == 0:
            continue

        peaks = np.zeros(intensity_arr.size, dtype=bool)
        if intensity_arr.size == 1:
            peaks[0] = bool(significant_mask[0])
        elif intensity_arr.size == 2:
            peaks[0] = bool(significant_mask[0] and intensity_arr[0] >= intensity_arr[1])
            peaks[1] = bool(significant_mask[1] and intensity_arr[1] > intensity_arr[0])
        else:
            peaks[0] = bool(significant_mask[0] and intensity_arr[0] > intensity_arr[1])
            peaks[-1] = bool(significant_mask[-1] and intensity_arr[-1] > intensity_arr[-2])
            peaks[1:-1] = (
                significant_mask[1:-1]
                & (intensity_arr[1:-1] >= intensity_arr[:-2])
                & (intensity_arr[1:-1] > intensity_arr[2:])
            )

        peak_count = int(np.count_nonzero(peaks))
        points_per_peak = float(significant_count / peak_count) if peak_count > 0 else float(significant_count)

        mz_diff = np.diff(mz_arr)
        mz_diff = mz_diff[np.isfinite(mz_diff) & (mz_diff > 0.0)]
        spacing_cv = (
            float(np.std(mz_diff) / np.mean(mz_diff))
            if mz_diff.size > 0 and float(np.mean(mz_diff)) > 0.0
            else np.inf
        )
        coverage_ratio = significant_count / mz_arr.size

        score = 0.0

        if points_per_peak >= 4.0:
            score += 3.0
        elif points_per_peak <= 2.5:
            score -= 3.0

        if coverage_ratio >= 0.35:
            score += 1.0
        elif coverage_ratio <= 0.10:
            score -= 1.0

        if spacing_cv <= 0.25:
            score += 1.0
        elif spacing_cv >= 0.60:
            score -= 1.0

        scores.append(score)

    if not scores:
        raise ValueError("Unable to infer spectrum type from sampled spectra.")

    return "profile" if float(np.mean(scores)) > 0.0 else "centroid"


def set_spectrum_type_metadata(meta: Any, spectrum_type: SpectrumType) -> None:
    """Write spectrum type back to metadata using canonical flag values."""
    if meta is None:
        raise ValueError("spectrum metadata is missing. Cannot set spectrum type.")

    if spectrum_type == "profile":
        meta.profile_spectrum = True
        meta.centroid_spectrum = None
        return

    if spectrum_type == "centroid":
        meta.profile_spectrum = None
        meta.centroid_spectrum = True
        return

    raise ValueError(f"Unsupported spectrum type: {spectrum_type!r}")
