from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray


CARBON_ISOTOPE_MASS_DIFF = 1.0033548378


@dataclass(slots=True)
class PeakReferenceStats:
    """Statistics aligned one-to-one with a reference m/z axis."""

    reference: NDArray[np.float64]
    counts: Optional[NDArray[np.int64]] = None
    intensity_sum: Optional[NDArray[np.float64]] = None
    intensity_max: Optional[NDArray[np.float64]] = None
    n_spectra: Optional[int] = None

    @property
    def frequency(self) -> Optional[NDArray[np.float64]]:
        """Return per-reference occurrence frequency when counts and n_spectra exist."""
        if self.counts is None or self.n_spectra is None or self.n_spectra <= 0:
            return None
        capped_counts = np.minimum(self.counts, int(self.n_spectra))
        return capped_counts.astype(np.float64) / float(self.n_spectra)

    @property
    def mean_intensity(self) -> Optional[NDArray[np.float64]]:
        """Return per-reference mean intensity over matched spectra."""
        if self.counts is None or self.intensity_sum is None:
            return None

        out = np.zeros(self.reference.size, dtype=np.float64)
        mask = self.counts > 0
        out[mask] = self.intensity_sum[mask] / self.counts[mask]
        return out


def _as_1d_array(name: str, values: Optional[np.ndarray], size: int) -> Optional[np.ndarray]:
    if values is None:
        return None

    arr = np.asarray(values)
    if arr.ndim != 1 or arr.size != size:
        raise ValueError(f"{name} must be a 1D array with the same length as reference")
    return arr


def _validate_stats(stats: PeakReferenceStats) -> PeakReferenceStats:
    reference = np.asarray(stats.reference, dtype=np.float64)
    if reference.ndim != 1:
        raise ValueError("reference must be a 1D array")

    size = reference.size
    counts = _as_1d_array("counts", stats.counts, size)
    intensity_sum = _as_1d_array("intensity_sum", stats.intensity_sum, size)
    intensity_max = _as_1d_array("intensity_max", stats.intensity_max, size)

    if counts is not None:
        counts = counts.astype(np.int64, copy=False)
        if np.any(counts < 0):
            raise ValueError("counts must contain non-negative values")

    if stats.n_spectra is not None and stats.n_spectra <= 0:
        raise ValueError("n_spectra must be positive when provided")

    return PeakReferenceStats(
        reference=reference,
        counts=counts,
        intensity_sum=None if intensity_sum is None else intensity_sum.astype(np.float64, copy=False),
        intensity_max=None if intensity_max is None else intensity_max.astype(np.float64, copy=False),
        n_spectra=stats.n_spectra,
    )


def _require_array(name: str, values: Optional[np.ndarray]) -> np.ndarray:
    if values is None:
        raise ValueError(f"{name} is required for this peak filter")
    return values


def _score_for_top_k(stats: PeakReferenceStats, top_k_by: str) -> np.ndarray:
    key = top_k_by.strip().lower()
    if key == "count":
        return _require_array("counts", stats.counts).astype(np.float64)
    if key == "frequency":
        return _require_array("frequency", stats.frequency)
    if key == "mean_intensity":
        return _require_array("mean_intensity", stats.mean_intensity)
    if key == "max_intensity":
        return _require_array("intensity_max", stats.intensity_max)
    if key == "total_intensity":
        return _require_array("intensity_sum", stats.intensity_sum)
    raise ValueError(
        "top_k_by must be one of: count, frequency, mean_intensity, "
        "max_intensity, total_intensity"
    )


def _apply_top_k_mask(keep: np.ndarray, score: np.ndarray, top_k: int) -> np.ndarray:
    if top_k < 1:
        raise ValueError("top_k must be a positive integer when provided")

    kept_idx = np.flatnonzero(keep)
    if kept_idx.size <= top_k:
        return keep

    kept_scores = np.asarray(score[kept_idx], dtype=np.float64)
    kept_scores = np.where(np.isfinite(kept_scores), kept_scores, -np.inf)
    selected_local = np.argpartition(kept_scores, -top_k)[-top_k:]
    selected = kept_idx[selected_local]

    out = np.zeros_like(keep, dtype=bool)
    out[selected] = True
    return out


def _isotope_score(stats: PeakReferenceStats, isotope_score: str) -> np.ndarray:
    key = isotope_score.strip().lower()
    if key == "count":
        return _require_array("counts", stats.counts).astype(np.float64)
    if key == "mean_intensity":
        return _require_array("mean_intensity", stats.mean_intensity)
    if key == "max_intensity":
        return _require_array("intensity_max", stats.intensity_max)
    if key == "total_intensity":
        return _require_array("intensity_sum", stats.intensity_sum)
    raise ValueError("isotope_score must be one of: count, mean_intensity, max_intensity, total_intensity")


def _remove_isotope_candidates(
    reference: np.ndarray,
    keep: np.ndarray,
    score: np.ndarray,
    *,
    charge_states: Sequence[int],
    mz_tolerance: float,
) -> np.ndarray:
    """Remove simple +C13 isotope candidates when the lighter peak has >= score."""
    if mz_tolerance < 0.0 or not np.isfinite(mz_tolerance):
        raise ValueError("isotope_mz_tolerance must be a finite non-negative number")

    out = keep.copy()
    active_idx = np.flatnonzero(keep)
    if active_idx.size <= 1:
        return out

    active_mz = reference[active_idx]
    if np.any(np.diff(active_mz) < 0):
        order = np.argsort(active_mz)
        active_idx = active_idx[order]
        active_mz = active_mz[order]

    active_score = np.asarray(score[active_idx], dtype=np.float64)
    active_score = np.where(np.isfinite(active_score), active_score, 0.0)

    for pos, mz_value in enumerate(active_mz):
        if not out[active_idx[pos]]:
            continue

        base_score = active_score[pos]
        for charge in charge_states:
            if charge <= 0:
                raise ValueError("isotope charge states must be positive integers")

            target = mz_value + CARBON_ISOTOPE_MASS_DIFF / float(charge)
            right = np.searchsorted(active_mz, target + mz_tolerance, side="right")
            left = np.searchsorted(active_mz, target - mz_tolerance, side="left")

            for hit in range(left, right):
                if hit == pos or not out[active_idx[hit]]:
                    continue
                if active_score[hit] <= base_score:
                    out[active_idx[hit]] = False

    return out


def filter_reference_peaks(
    stats: PeakReferenceStats,
    *,
    min_count: Optional[int] = None,
    min_frequency: Optional[float] = None,
    min_mean_intensity: Optional[float] = None,
    min_max_intensity: Optional[float] = None,
    top_k: Optional[int] = None,
    top_k_by: str = "mean_intensity",
    remove_isotopes: bool = False,
    isotope_charge_states: Sequence[int] = (1, 2),
    isotope_mz_tolerance: float = 0.01,
    isotope_score: str = "mean_intensity",
) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """Filter a reference m/z axis using precomputed per-reference statistics."""
    stats = _validate_stats(stats)
    keep = np.ones(stats.reference.size, dtype=bool)

    if min_count is not None:
        if min_count < 1:
            raise ValueError("min_count must be a positive integer when provided")
        keep &= _require_array("counts", stats.counts) >= int(min_count)

    if min_frequency is not None:
        if not np.isfinite(min_frequency) or min_frequency < 0.0 or min_frequency > 1.0:
            raise ValueError("min_frequency must be in [0, 1] when provided")
        keep &= _require_array("frequency", stats.frequency) >= float(min_frequency)

    if min_mean_intensity is not None:
        keep &= _require_array("mean_intensity", stats.mean_intensity) >= float(min_mean_intensity)

    if min_max_intensity is not None:
        keep &= _require_array("intensity_max", stats.intensity_max) >= float(min_max_intensity)

    if top_k is not None:
        keep = _apply_top_k_mask(keep, _score_for_top_k(stats, top_k_by), int(top_k))

    if remove_isotopes:
        keep = _remove_isotope_candidates(
            stats.reference,
            keep,
            _isotope_score(stats, isotope_score),
            charge_states=isotope_charge_states,
            mz_tolerance=float(isotope_mz_tolerance),
        )

    return stats.reference[keep], keep
