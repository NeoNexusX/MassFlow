"""Plotting utilities for inspecting peak-filter statistics."""

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.helper.stat_filter_helper import (
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
    intensity_log_x: bool = True,
    log_y: bool = False,
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
    count_plot_values = _plot_values(count_values, log_x)
    intensity_plot_values = _plot_values(
        intensity_values,
        log_x or intensity_log_x,
    )
    histogram_style = {
        "alpha": 0.85,
        "edgecolor": "white",
        "linewidth": 0.5,
    }
    count_ax.hist(
        count_plot_values,
        bins=_count_bins(count_display, profile.spectrum_count, bins),
        color="#4C78A8",
        **histogram_style,
    )
    intensity_ax.hist(
        intensity_plot_values,
        bins=_intensity_bins(intensity_plot_values, bins, log_x or intensity_log_x),
        color="#F58518",
        **histogram_style,
    )

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
        if log_y:
            axis.set_yscale("log")
        axis.grid(axis="y", alpha=0.25, linewidth=0.8)
        axis.set_axisbelow(True)
        if axis.lines:
            axis.legend()

    if count_display == "frequency" and not log_x:
        count_ax.set_xlim(0, 1)
    if log_x:
        count_ax.set_xscale("log")
    if log_x or intensity_log_x:
        intensity_ax.set_xscale("log")
        zero_count = int(np.count_nonzero(intensity_values == 0))
        if zero_count:
            intensity_ax.text(
                0.98,
                0.95,
                f"{zero_count:,} zero-intensity positions omitted",
                ha="right",
                va="top",
                transform=intensity_ax.transAxes,
                fontsize="small",
                color="0.35",
            )

    total_positions = profile.counts.size
    fig.suptitle(
        f"Peak filter statistics ({total_positions:,} m/z positions, "
        f"{profile.spectrum_count:,} spectra)",
        fontsize="large",
    )
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


def _count_bins(
    count_display: CountDisplay,
    spectrum_count: int,
    bins: int,
) -> np.ndarray | int:
    """Return bins matching the discrete count or bounded frequency domain."""
    if count_display == "frequency":
        return np.linspace(0, 1, bins + 1)
    if spectrum_count <= bins:
        return np.arange(-0.5, spectrum_count + 1.5)
    return bins


def _intensity_bins(
    values: np.ndarray,
    bins: int,
    log_x: bool,
) -> np.ndarray | int:
    """Return geometric bins when intensity is displayed on a logarithmic axis."""
    if not log_x:
        return bins
    minimum = float(values.min())
    maximum = float(values.max())
    if minimum == maximum:
        return bins
    return np.geomspace(minimum, maximum, bins + 1)
