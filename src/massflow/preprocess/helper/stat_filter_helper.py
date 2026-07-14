"""Peak filtering helper functions for mass spectrometry data preprocessing."""

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal

import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.numba.stat_filter_numba import (
    aggregate_peak_profile_flat_parallel_jit,
    count_peaks_flat_parallel_jit,
    filter_peaks_flat_parallel_jit,
    sum_intensity_flat_parallel_jit,
)
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.peak_filter_helper")

IntensityMethod = Literal["mean", "sum"]


@dataclass(frozen=True)
class PeakFilterProfile:
    """Reusable count and intensity statistics for one continuous dataset."""

    counts: np.ndarray
    summed_intensity: np.ndarray
    spectrum_count: int
    mz_axis: np.ndarray

    @property
    def frequencies(self) -> np.ndarray:
        """Return the fraction of spectra containing each peak."""
        return self.counts / self.spectrum_count

    @property
    def mean_intensity(self) -> np.ndarray:
        """Return mean peak intensity across all spectra, including zeros."""
        return self.summed_intensity / self.spectrum_count

    def intensity(self, method: IntensityMethod) -> np.ndarray:
        """Return the requested aggregate intensity spectrum."""
        if method == "mean":
            return self.mean_intensity
        if method == "sum":
            return self.summed_intensity
        raise ValueError("method must be either 'mean' or 'sum'.")


# ---------------------------------------------------------------
# Shared profile and filtering helpers
# ---------------------------------------------------------------
def _validate_continuous_data(data_manager: MSDataManagerImzML) -> None:
    """Validate that filtering can use one shared m/z mask."""
    if data_manager.ms.meta.continuous is not True:
        raise ValueError("Peak filtering is only applicable to continuous data.")


def _aggregate_batches(
    data_manager: MSDataManagerImzML,
    aggregate_batch: Callable[
        [np.ndarray, np.ndarray, int],
        tuple[np.ndarray, int],
    ],
    *,
    batch_size: int = 256,
) -> tuple[np.ndarray, int]:
    """Aggregate one statistic and the spectrum count across flat batches."""
    aggregate = None
    spectrum_count = 0

    for mz_axis, intensity_flat, lengths, _ in data_manager.flat_generator(
        batch_size=batch_size,
        include_mz=True,
    ):
        if lengths.size and int(np.max(lengths)) > mz_axis.size:
            raise ValueError(
                "continuous spectrum length exceeds the shared m/z axis length."
            )
        batch_aggregate, batch_spectrum_count = aggregate_batch(
            intensity_flat,
            lengths,
            mz_axis.size,
        )
        if aggregate is None:
            aggregate = batch_aggregate
        else:
            aggregate += batch_aggregate
        spectrum_count += batch_spectrum_count

    if aggregate is None or spectrum_count == 0:
        raise ValueError("Cannot aggregate peaks from empty data.")

    return aggregate, spectrum_count


def compute_peak_filter_profile(
    data_manager: MSDataManagerImzML,
    *,
    batch_size: int = 256,
) -> PeakFilterProfile:
    """Compute count and summed-intensity statistics in one dataset scan."""
    _validate_continuous_data(data_manager)

    counts = None
    summed_intensity = None
    spectrum_count = 0

    for mz_axis, intensity_flat, lengths, _ in data_manager.flat_generator(
        batch_size=batch_size,
        include_mz=True,
    ):
        batch_counts, batch_sum, batch_spectrum_count = (
            aggregate_peak_profile_flat_parallel_jit(
                intensity_flat,
                lengths,
                mz_axis.size,
            )
        )
        if counts is None:
            counts = batch_counts
            summed_intensity = batch_sum
        else:
            counts += batch_counts
            summed_intensity += batch_sum # type: ignore
        spectrum_count += batch_spectrum_count

    if counts is None or summed_intensity is None or spectrum_count == 0:
        raise ValueError("Cannot compute a peak filter profile from empty data.")

    return PeakFilterProfile(
        counts=counts,
        summed_intensity=summed_intensity,
        spectrum_count=spectrum_count,
        mz_axis=np.asarray(data_manager.ms.shared_mz_list).copy(),
    )


def _validate_profile(
    data_manager: MSDataManagerImzML,
    profile: PeakFilterProfile,
) -> None:
    """Validate that a profile can be reused with the data manager."""
    if not isinstance(profile, PeakFilterProfile):
        raise TypeError("profile must be a PeakFilterProfile.")
    if profile.counts.ndim != 1 or profile.summed_intensity.ndim != 1:
        raise ValueError("profile counts and summed_intensity must be 1-D arrays.")
    if (
        profile.counts.size != profile.summed_intensity.size
        or profile.counts.size != profile.mz_axis.size
    ):
        raise ValueError("profile arrays must have the same length.")
    if profile.spectrum_count <= 0:
        raise ValueError("profile spectrum_count must be positive.")

    shared_mz = data_manager.ms.shared_mz_list
    if shared_mz is None or not np.array_equal(profile.mz_axis, shared_mz):
        raise ValueError("profile does not match the data manager shared m/z axis.")


def _resolve_profile(
    data_manager: MSDataManagerImzML,
    profile: PeakFilterProfile | None,
    *,
    batch_size: int = 256,
) -> PeakFilterProfile:
    """Return a validated reusable profile, computing it when omitted."""
    if profile is None:
        return compute_peak_filter_profile(data_manager, batch_size=batch_size)
    _validate_continuous_data(data_manager)
    _validate_profile(data_manager, profile)
    return profile


def filter_data_by_mask(
    data_manager: MSDataManagerImzML,
    keep_mask: np.ndarray,
    *,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
    batch_size: int = 256,
) -> MSDataManagerImzML:
    """Apply one shared peak mask and write a filtered continuous dataset.

    Public so other dataset-scope helpers (e.g. ``adduct_filter_helper``)
    can reuse the same shared-mask write-back path without duplicating logic.
    """
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    filtered_data = MSDataManagerImzML(
        filepath=None if output_path is None else str(output_path),
        temp_dir=None if temp_dir is None else str(temp_dir),
        max_threads=data_manager.max_threads,
        mz_dtype=data_manager.mz_dtype,
        intensity_dtype=data_manager.intensity_dtype,
    )
    filtered_data.copy_meta(data_manager)
    filtered_data.ms.meta.continuous = True
    filtered_data.ms.meta.processed = False

    filtered_mz = None
    for mz_axis, intensity_flat, lengths, coordinates in data_manager.flat_generator(
        batch_size=batch_size,
        include_mz=True,
    ):
        if filtered_mz is None:
            filtered_mz = mz_axis[keep_mask]
        filtered_intensity, filtered_lengths = filter_peaks_flat_parallel_jit(
            intensity_flat,
            lengths,
            keep_mask,
        )
        filtered_data.write_flat_batch(
            intensity_flat=filtered_intensity,
            lengths=filtered_lengths,
            coordinates=coordinates,
            mz_axis=filtered_mz,
        )

    filtered_data.close_writer()
    filtered_data.load_head_data()
    return filtered_data


# ---------------------------------------------------------------
# Peak-selection masks
# ---------------------------------------------------------------
def _build_count_mask(
    profile: PeakFilterProfile,
    *,
    min_count: int | None = None,
    min_frequency: float | None = None,
) -> np.ndarray:
    """Build a peak mask from occurrence count and frequency thresholds."""
    if min_count is not None and (
        not isinstance(min_count, (int, np.integer)) or min_count < 1
    ):
        raise ValueError("min_count must be a positive integer when provided.")
    if min_frequency is not None and (
        not np.isfinite(min_frequency) or not 0 <= min_frequency <= 1
    ):
        raise ValueError("min_frequency must be in [0, 1] when provided.")

    keep_mask = np.ones(profile.counts.size, dtype=bool)
    if min_count is not None:
        keep_mask &= profile.counts >= min_count
    if min_frequency is not None:
        keep_mask &= profile.frequencies >= min_frequency
    return keep_mask


def _build_intensity_mask(
    profile: PeakFilterProfile,
    *,
    method: IntensityMethod = "mean",
    min_intensity: float | None = None,
    quantile: float | None = None,
) -> np.ndarray:
    """Build a peak mask from aggregate intensity thresholds."""
    if method not in {"mean", "sum"}:
        raise ValueError("Method must be either 'mean' or 'sum'.")
    if min_intensity is not None and (
        not np.isfinite(min_intensity) or min_intensity < 0
    ):
        raise ValueError("min_intensity must be a finite non-negative number.")
    if quantile is not None and not 0 < quantile < 1:
        raise ValueError("Quantile must be between 0 and 1.")

    intensity = profile.intensity(method)
    keep_mask = np.ones(intensity.size, dtype=bool)
    if min_intensity is not None:
        keep_mask &= intensity >= min_intensity
    if quantile is not None:
        keep_mask &= intensity > np.quantile(intensity, quantile)
    return keep_mask


# ---------------------------------------------------------------
# Public peak filtering operations of statistics-based filtering and profile computation
# ---------------------------------------------------------------
def filter_peaks_by_stats(
    data_manager: MSDataManagerImzML,
    *,
    profile: PeakFilterProfile | None = None,
    min_count: int | None = None,
    min_frequency: float | None = None,
    intensity_method: IntensityMethod = "mean",
    min_intensity: float | None = None,
    intensity_quantile: float | None = None,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
    batch_size: int = 256,
) -> MSDataManagerImzML:
    """Filter peaks using count and intensity masks combined with AND."""
    if (
        min_count is None
        and min_frequency is None
        and min_intensity is None
        and intensity_quantile is None
    ):
        raise ValueError("At least one peak filtering threshold must be provided.")

    profile = _resolve_profile(data_manager, profile, batch_size=batch_size)
    count_mask = _build_count_mask(
        profile,
        min_count=min_count,
        min_frequency=min_frequency,
    )
    intensity_mask = _build_intensity_mask(
        profile,
        method=intensity_method,
        min_intensity=min_intensity,
        quantile=intensity_quantile,
    )
    return filter_data_by_mask(
        data_manager,
        count_mask & intensity_mask,
        output_path=output_path,
        temp_dir=temp_dir,
        batch_size=batch_size,
    )


def filter_peaks_by_count(
    data_manager: MSDataManagerImzML,
    *,
    min_count: int | None = None,
    min_frequency: float | None = None,
    profile: PeakFilterProfile | None = None,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
) -> MSDataManagerImzML:
    """Filter peaks using occurrence count and frequency thresholds."""
    if min_count is None and min_frequency is None:
        raise ValueError("Either min_count or min_frequency must be provided.")

    profile = _resolve_profile(data_manager, profile)
    keep_mask = _build_count_mask(
        profile,
        min_count=min_count,
        min_frequency=min_frequency,
    )
    return filter_data_by_mask(
        data_manager,
        keep_mask,
        output_path=output_path,
        temp_dir=temp_dir,
    )


def filter_peaks_by_intensity(
    data_manager: MSDataManagerImzML,
    *,
    method: IntensityMethod = "mean",
    min_intensity: float | None = None,
    intensity_quantile: float | None = 0.05,
    profile: PeakFilterProfile | None = None,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
) -> MSDataManagerImzML:
    """Filter peaks using aggregate intensity thresholds."""
    if min_intensity is None and intensity_quantile is None:
        raise ValueError("Either min_intensity or intensity_quantile must be provided.")

    profile = _resolve_profile(data_manager, profile)
    keep_mask = _build_intensity_mask(
        profile,
        method=method,
        min_intensity=min_intensity,
        quantile=intensity_quantile,
    )
    return filter_data_by_mask(
        data_manager,
        keep_mask,
        output_path=output_path,
        temp_dir=temp_dir,
    )


def compute_sum_intensity(
    data_manager: MSDataManagerImzML,
    *,
    batch_size: int = 256,
) -> np.ndarray:
    """Compute the sum spectrum across all spectra in the data manager."""
    summed, _ = _aggregate_batches(
        data_manager,
        sum_intensity_flat_parallel_jit,
        batch_size=batch_size,
    )
    return summed


def compute_mean_intensity(data_manager: MSDataManagerImzML) -> np.ndarray:
    """Compute the mean spectrum across all spectra in the data manager."""
    summed, spectrum_count = _aggregate_batches(
        data_manager,
        sum_intensity_flat_parallel_jit,
    )
    return summed / spectrum_count


def compute_peak_counts(data_manager: MSDataManagerImzML) -> np.ndarray:
    """Compute the occurrence count spectrum across all spectra."""
    counts, _ = _aggregate_batches(data_manager, count_peaks_flat_parallel_jit)
    return counts
