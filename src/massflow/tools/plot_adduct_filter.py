"""Plot auditable adduct-filter decisions before modifying an MSI dataset."""

from pathlib import Path

from matplotlib.axes import Axes
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.helper.adduct_filter_helper import (
    AdductAnnotationResult,
    AdductRule,
    Polarity,
    _select_summary_annotation_indices,
    _validate_annotation_result,
    compute_adduct_annotations,
    write_adduct_filter_report,
)

# Stable, colorblind-friendly rule colors and representative-rule markers.
_RULE_COLORS = (
    "#4C78A8",
    "#F58518",
    "#54A24B",
    "#E45756",
    "#72B7B2",
    "#B279A2",
    "#FF9DA6",
    "#9D755D",
)
_RULE_MARKERS = ("o", "s", "^", "D", "P", "X", "v", "<")
_PASS_COLOR = "#54A24B"
_REMOVAL_EVIDENCE_COLOR = "#E45756"
_CANDIDATE_COLOR = "#9E9E9E"


# ---------------------------------------------------------------------------
# Public preview operations
# ---------------------------------------------------------------------------


def plot_adduct_filter(  # pylint: disable=too-many-arguments
    data_manager: MSDataManagerImzML,
    *,
    polarity: Polarity = "positive",
    rules: tuple[AdductRule, ...] | None = None,
    ppm: float = 5.0,
    mzabs: float = 0.015,
    require_quasi: bool = True,
    min_spatial_correlation: float = 0.7,
    batch_size: int = 256,
    show: bool = True,
    output_path: str | Path | None = None,
    report_path: str | Path | None = None,
    dpi: int = 150,
    figsize: tuple[float, float] = (15, 7),
) -> AdductAnnotationResult:
    """Compute and plot a reusable adduct-filter preview without removing peaks.

    The returned result can be passed to ``apply_adduct_annotations`` so the
    reviewed decisions are applied without recomputing hypotheses or spatial
    evidence.
    """
    _validate_plot_settings(dpi, figsize)
    result = compute_adduct_annotations(
        data_manager,
        polarity=polarity,
        rules=rules,
        ppm=ppm,
        mzabs=mzabs,
        require_quasi=require_quasi,
        min_spatial_correlation=min_spatial_correlation,
        batch_size=batch_size,
    )
    if report_path is not None:
        write_adduct_filter_report(result, report_path)
    return _render(result, show, output_path, dpi, figsize)


def plot_adduct_filter_result(  # pylint: disable=too-many-arguments
    result: AdductAnnotationResult,
    *,
    show: bool = True,
    output_path: str | Path | None = None,
    report_path: str | Path | None = None,
    dpi: int = 150,
    figsize: tuple[float, float] = (15, 7),
) -> AdductAnnotationResult:
    """Render an existing result using its recorded decision threshold."""
    _validate_plot_settings(dpi, figsize)
    if report_path is not None:
        write_adduct_filter_report(result, report_path)
    return _render(result, show, output_path, dpi, figsize)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _render(
    result: AdductAnnotationResult,
    show: bool,
    output_path: str | Path | None,
    dpi: int,
    figsize: tuple[float, float],
) -> AdductAnnotationResult:
    _validate_annotation_result(result)
    fig, (removal_ax, evidence_ax) = plt.subplots(
        1,
        2,
        figsize=figsize,
        dpi=dpi,
        gridspec_kw={"width_ratios": (1.08, 1.0)},
    )

    _plot_removal_summary(removal_ax, result)
    _plot_correlation_distribution(evidence_ax, result)
    fig.suptitle(_summary_title(result), fontsize="large", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return result


def _plot_removal_summary(ax: Axes, result: AdductAnnotationResult) -> None:
    """Plot one strongest-evidence row for every uniquely removed peak."""
    summary_indices = _select_summary_annotation_indices(result.annotations)
    annotations = [
        annotation
        for index, annotation in enumerate(result.annotations)
        if index in summary_indices
    ]
    rule_colors = _rule_colors(result.rule_names)
    rule_markers = _rule_markers(result.rule_names)

    if not annotations:
        ax.text(
            0.5,
            0.5,
            "No peaks would be removed\nunder the recorded settings",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="0.35",
            fontsize="large",
        )
        _set_empty_mz_limits(ax, result.mz_axis)
    else:
        for representative_rule in result.rule_names:
            for removed_rule in result.rule_names:
                matching = [
                    annotation
                    for annotation in annotations
                    if annotation.representative_rule == representative_rule
                    and annotation.removed_rule == removed_rule
                ]
                if not matching:
                    continue
                ax.scatter(
                    [annotation.removed_mz for annotation in matching],
                    [
                        annotation.removed_mz - annotation.representative_mz
                        for annotation in matching
                    ],
                    s=38,
                    marker=rule_markers[representative_rule],
                    color=rule_colors[removed_rule],
                    alpha=0.72,
                    linewidths=0.5,
                    edgecolors="black",
                    rasterized=True,
                )

        removed_rules = {annotation.removed_rule for annotation in annotations}
        representative_rules = {
            annotation.representative_rule for annotation in annotations
        }
        color_legend = ax.legend(
            handles=[
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="None",
                    markerfacecolor=rule_colors[name],
                    markeredgecolor="none",
                    markersize=7,
                    label=name,
                )
                for name in result.rule_names
                if name in removed_rules
            ],
            loc="upper right",
            fontsize="x-small",
            title="Removed-peak rule (color)",
            title_fontsize="x-small",
        )
        ax.add_artist(color_legend)
        ax.legend(
            handles=[
                Line2D(
                    [],
                    [],
                    marker=rule_markers[name],
                    linestyle="None",
                    markerfacecolor="0.75",
                    markeredgecolor="black",
                    markersize=7,
                    label=name,
                )
                for name in result.rule_names
                if name in representative_rules
            ],
            loc="lower right",
            fontsize="x-small",
            title="Kept-peak rule (marker)",
            title_fontsize="x-small",
        )

    ax.axhline(0.0, color="0.35", linewidth=0.9, linestyle="--")
    ax.set_xlabel("removed peak m/z")
    ax.set_ylabel("removed m/z − kept representative m/z")
    alternative_count = len(result.annotations) - len(annotations)
    ax.set_title(
        f"Unique removed peaks ({len(annotations):,}; "
        f"{alternative_count:,} alternative evidence rows in CSV)"
    )
    ax.grid(alpha=0.25, linewidth=0.8)
    ax.set_axisbelow(True)


def _plot_correlation_distribution(ax: Axes, result: AdductAnnotationResult) -> None:
    """Plot candidate, passed, and final-summary correlation distributions."""
    threshold = result.config.min_spatial_correlation
    candidate_correlations = np.asarray(
        [
            evidence.spatial_correlation
            for evidence in result.pair_evidence
            if np.isfinite(evidence.spatial_correlation)
        ],
        dtype=np.float64,
    )
    passed_correlations = np.asarray(
        [
            evidence.spatial_correlation
            for evidence in result.pair_evidence
            if evidence.passes_threshold
        ],
        dtype=np.float64,
    )
    summary_indices = _select_summary_annotation_indices(result.annotations)
    summary_correlations = np.asarray(
        [
            annotation.spatial_correlation
            for index, annotation in enumerate(result.annotations)
            if index in summary_indices
        ],
        dtype=np.float64,
    )

    ax.axvline(
        threshold,
        linestyle="--",
        color="black",
        linewidth=1.2,
        label=f"Recorded threshold = {threshold:g}",
    )
    if candidate_correlations.size == 0:
        ax.text(
            0.5,
            0.5,
            "No candidate adduct relationships",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="0.35",
            fontsize="large",
        )
    else:
        bins = np.linspace(-1.0, 1.0, 61)
        ax.hist(
            candidate_correlations,
            bins=bins,
            color=_CANDIDATE_COLOR,
            alpha=0.5,
            label=f"All finite candidates ({candidate_correlations.size:,})",
        )
        ax.hist(
            passed_correlations,
            bins=bins,
            histtype="step",
            color=_PASS_COLOR,
            linewidth=1.8,
            label=f"Passed threshold ({passed_correlations.size:,})",
        )
        ax.hist(
            summary_correlations,
            bins=bins,
            histtype="step",
            color=_REMOVAL_EVIDENCE_COLOR,
            linewidth=1.8,
            label=f"Unique-removal evidence ({summary_correlations.size:,})",
        )
        ax.set_yscale("log")

    ax.set_xlim(-1.0, 1.0)
    ax.set_xlabel("spatial Pearson correlation")
    ax.set_ylabel("relationship count (log scale)")
    ax.set_title("Spatial-correlation distributions")
    ax.grid(axis="y", alpha=0.25, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize="x-small")


# ---------------------------------------------------------------------------
# Plot-data formatting
# ---------------------------------------------------------------------------


def _summary_title(result: AdductAnnotationResult) -> str:
    removed_peaks = {annotation.removed_peak for annotation in result.annotations}
    total_peaks = int(result.mz_axis.size)
    removal_fraction = len(removed_peaks) / total_peaks if total_peaks else 0.0
    passed_count = sum(evidence.passes_threshold for evidence in result.pair_evidence)
    config = result.config
    counts = (
        f"Adduct filter preview: {total_peaks:,} peaks | "
        f"{len(result.pair_evidence):,} candidate relationships | "
        f"{passed_count:,} passed threshold | "
        f"{len(removed_peaks):,} unique peaks removed ({removal_fraction:.1%}; "
        f"{len(result.annotations):,} evidence rows)"
    )
    settings = (
        f"polarity={config.polarity}; ppm={config.ppm:g}; mzabs={config.mzabs:g}; "
        f"require_quasi={config.require_quasi}; "
        f"min_spatial_correlation={config.min_spatial_correlation:g}; "
        f"batch_size={config.batch_size}"
    )
    return f"{counts}\n{settings}"


def _rule_colors(rule_names: tuple[str, ...]) -> dict[str, str]:
    return {
        name: _RULE_COLORS[index % len(_RULE_COLORS)]
        for index, name in enumerate(rule_names)
    }


def _rule_markers(rule_names: tuple[str, ...]) -> dict[str, str]:
    return {
        name: _RULE_MARKERS[index % len(_RULE_MARKERS)]
        for index, name in enumerate(rule_names)
    }


def _set_empty_mz_limits(ax: Axes, mz_axis: np.ndarray) -> None:
    finite = np.asarray(mz_axis, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return
    minimum = float(finite.min())
    maximum = float(finite.max())
    if minimum == maximum:
        padding = max(abs(minimum) * 0.01, 1.0)
        ax.set_xlim(minimum - padding, maximum + padding)
    else:
        padding = (maximum - minimum) * 0.02
        ax.set_xlim(minimum - padding, maximum + padding)


def _validate_plot_settings(
    dpi: int,
    figsize: tuple[float, float],
) -> None:
    if not isinstance(dpi, (int, np.integer)) or dpi <= 0:
        raise ValueError("dpi must be a positive integer.")
    if not (len(figsize) == 2 and all(np.isfinite(d) and d > 0 for d in figsize)):
        raise ValueError("figsize must be a (width, height) pair of positive floats.")
