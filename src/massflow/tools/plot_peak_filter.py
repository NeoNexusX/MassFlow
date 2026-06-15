"""Plotting utilities for inspecting peak-filter statistics."""

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.helper.peak_filter_helper import (
    IntensityMethod,
    PeakFilterProfile,
    compute_peak_filter_profile,
)

CountDisplay = Literal["count", "frequency"]


def plot_peak_filter_histogram(
    data_manager: MSDataManagerImzML,
    *,
    count_display: CountDisplay = "frequency",
    intensity_method: IntensityMethod = "mean",
    min_count: int | None = None,
    min_frequency: float | None = None,
    min_intensity: float | None = None,
    quantile: float | None = None,
    bins: int = 100,
    log_x: bool = False,
    log_y: bool = True,
    show: bool = True,
    save_path: str | Path | None = None,
    dpi: int = 150,
    figsize: tuple[float, float] = (14, 5),
) -> PeakFilterProfile:
    """Plot count/frequency and intensity histograms, returning the reusable profile."""
    if count_display not in {"count", "frequency"}:
        raise ValueError("count_display must be either 'count' or 'frequency'.")
    if intensity_method not in {"mean", "sum"}:
        raise ValueError("intensity_method must be either 'mean' or 'sum'.")
    if not isinstance(bins, int) or bins < 1:
        raise ValueError("bins must be a positive integer.")
    if min_count is not None and (
        not isinstance(min_count, (int, np.integer)) or min_count < 1
    ):
        raise ValueError("min_count must be a positive integer when provided.")
    if min_frequency is not None and (
        not np.isfinite(min_frequency) or not 0 <= min_frequency <= 1
    ):
        raise ValueError("min_frequency must be in [0, 1] when provided.")
    if min_intensity is not None and (
        not np.isfinite(min_intensity) or min_intensity < 0
    ):
        raise ValueError("min_intensity must be a finite non-negative number.")
    if quantile is not None and not 0 < quantile < 1:
        raise ValueError("quantile must be between 0 and 1.")

    profile = compute_peak_filter_profile(data_manager)

    count_values = (
        profile.counts if count_display == "count" else profile.frequencies
    )
    intensity_values = profile.intensity(intensity_method)

    fig, (count_ax, intensity_ax) = plt.subplots(
        1,
        2,
        figsize=figsize,
        dpi=dpi,
    )
    count_ax.hist(_plot_values(count_values, log_x), bins=bins)
    intensity_ax.hist(_plot_values(intensity_values, log_x), bins=bins)

    if min_count is not None:
        count_threshold = (
            float(min_count)
            if count_display == "count"
            else float(min_count) / profile.spectrum_count
        )
        count_ax.axvline(count_threshold, linestyle="--", label=f"min_count={min_count}")
    if min_frequency is not None:
        frequency_threshold = (
            float(min_frequency) * profile.spectrum_count
            if count_display == "count"
            else float(min_frequency)
        )
        count_ax.axvline(
            frequency_threshold,
            linestyle="--",
            label=f"min_frequency={min_frequency}",
        )
    if min_intensity is not None:
        intensity_ax.axvline(
            min_intensity,
            linestyle="--",
            label=f"min_intensity={min_intensity}",
        )
    if quantile is not None:
        intensity_ax.axvline(
            np.quantile(intensity_values, quantile),
            linestyle="--",
            label=f"quantile={quantile}",
        )

    count_ax.set_title(f"Peak occurrence {count_display}")
    count_ax.set_xlabel(count_display)
    count_ax.set_ylabel("Number of m/z positions")
    intensity_ax.set_title(f"Peak {intensity_method} intensity")
    intensity_ax.set_xlabel(f"{intensity_method} intensity")
    intensity_ax.set_ylabel("Number of m/z positions")

    for axis in (count_ax, intensity_ax):
        if log_x:
            axis.set_xscale("log")
        if log_y:
            axis.set_yscale("log")
        if axis.lines:
            axis.legend()

    fig.tight_layout()
    if save_path is not None:
        output_path = Path(save_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return profile


def _plot_values(values: np.ndarray, log_x: bool) -> np.ndarray:
    """Return finite values suitable for histogram rendering."""
    finite_values = values[np.isfinite(values)]
    if log_x:
        finite_values = finite_values[finite_values > 0]
    if finite_values.size == 0:
        raise ValueError("profile contains no finite values available for plotting.")
    return finite_values
