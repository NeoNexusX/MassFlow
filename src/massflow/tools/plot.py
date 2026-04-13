from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from massflow.module import Spectrum
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tools")


@dataclass(frozen=True)
class SeriesStyle:
    color: str
    line_width: float
    mode: str
    label: Optional[str] = None


@dataclass(frozen=True)
class PlotConfig:
    save_path: Optional[str]
    figsize: Tuple[float, float]
    dpi: int
    mz_range: Optional[Tuple[float, float]]
    intensity_range: Optional[Tuple[float, float]]
    metrics_box: bool


class SpectrumPlotter:
    """Unified plotter for single and paired spectrum visualization."""

    DEFAULT_COLORS = ("#EECF31", "#4AE38C")
    ZERO_EPS = 1e-12

    def __init__(self, config: PlotConfig):
        self.config = config

    @staticmethod
    def _normalize_pair(values, default_pair):
        if values is None:
            return [default_pair[0], default_pair[1]]
        if isinstance(values, str):
            return [values, values]

        seq = list(values)
        if len(seq) == 0:
            return [default_pair[0], default_pair[1]]
        if len(seq) == 1:
            return [seq[0], seq[0]]
        return [seq[0], seq[1]]

    @staticmethod
    def _normalize_mode(mode: Optional[str]) -> str:
        value = (mode or "line").strip().lower()
        return value if value in {"line", "stem"} else "line"

    @staticmethod
    def _validate_spectrum(name: str, spectrum: Optional[Spectrum]) -> Spectrum:
        if spectrum is None:
            raise ValueError(f"'{name}' spectrum cannot be None.")
        if len(spectrum) == 0:
            raise ValueError(f"'{name}' spectrum has no data points.")
        return spectrum

    @staticmethod
    def _slice_arrays(
        spectrum: Spectrum,
        mz_range: Optional[Tuple[float, float]],
    ) -> Tuple[np.ndarray, np.ndarray]:
        if mz_range is None:
            mz = np.asarray(spectrum.mz_list)
            intensity = np.asarray(spectrum.intensity)
        else:
            cropped = spectrum.crop_range(mz_range=mz_range, sort_by_mz=True, mode="new")
            mz = np.asarray(cropped.mz_list)
            intensity = np.asarray(cropped.intensity)

        if mz.shape[0] != intensity.shape[0]:
            raise ValueError("m/z array and intensity array length mismatch.")
        if mz.size == 0:
            raise ValueError("No data points in selected 'mz_range'.")

        return mz, intensity

    @staticmethod
    def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
        if a.size <= 1 or b.size <= 1:
            return 0.0
        std_a = float(np.std(a))
        std_b = float(np.std(b))
        if std_a <= SpectrumPlotter.ZERO_EPS or std_b <= SpectrumPlotter.ZERO_EPS:
            return 0.0
        value = float(np.corrcoef(a, b)[0, 1])
        return 0.0 if np.isnan(value) else value

    @staticmethod
    def _safe_snr(spectrum: Spectrum) -> float:
        try:
            value = float(SpectrumPreprocess.calculate_snr_spectrum(spectrum))
            return value if np.isfinite(value) else 0.0
        except Exception as exc:
            logger.warning(f"SNR calculation failed: {exc}")
            return 0.0

    @staticmethod
    def _draw(ax, mz: np.ndarray, intensity: np.ndarray, style: SeriesStyle) -> None:
        if style.mode == "line":
            ax.plot(
                mz,
                intensity,
                color=style.color,
                linewidth=style.line_width,
                alpha=0.85,
                label=style.label,
            )
            return

        markerline, stemlines, baseline = ax.stem(mz, intensity, label=style.label)
        plt.setp(stemlines, color=style.color, linewidth=style.line_width, alpha=1)
        plt.setp(markerline, color=style.color, markersize=3, alpha=1)
        plt.setp(baseline, color="gray", linewidth=0.5, alpha=1)

    @staticmethod
    def _style_axes(ax) -> None:
        # Keep figure clean by hiding top/right spines.
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def _axis_limits(
        self,
        mz_a: np.ndarray,
        int_a: np.ndarray,
        mz_b: Optional[np.ndarray] = None,
        int_b: Optional[np.ndarray] = None,
    ) -> Tuple[float, float, float, float]:
        if self.config.mz_range is not None:
            x_min = float(self.config.mz_range[0])
            x_max = float(self.config.mz_range[1])
        else:
            x_data = [mz_a]
            if mz_b is not None:
                x_data.append(mz_b)
            x_all = np.concatenate(x_data)
            x_min, x_max = float(np.min(x_all)), float(np.max(x_all))

        if self.config.intensity_range is not None:
            y_min = float(self.config.intensity_range[0])
            y_max = float(self.config.intensity_range[1])
        else:
            y_data = [int_a]
            if int_b is not None:
                y_data.append(int_b)
            y_all = np.concatenate(y_data)
            y_min = 0.0
            y_max = max(float(np.max(y_all)) * 1.05, 1.0)

        return x_min, x_max, y_min, y_max

    @staticmethod
    def _apply_axis_limits(axes, limits: Tuple[float, float, float, float]) -> None:
        x_min, x_max, y_min, y_max = limits
        if not isinstance(axes, (list, tuple, np.ndarray)):
            axes = (axes,)

        for ax in axes:
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)

    def _metrics_text(self, base: Spectrum, target: Spectrum) -> str:
        n = min(len(base), len(target))
        base_intensity = np.asarray(base.intensity)
        target_intensity = np.asarray(target.intensity)
        base_mz = np.asarray(base.mz_list)

        if base_intensity.size == 0 or target_intensity.size == 0 or base_mz.size == 0:
            return "Metrics unavailable: empty spectrum window."

        o = base_intensity[:n]
        d = target_intensity[:n]

        corr = self._safe_corr(o, d)

        o_sum = float(np.sum(o))
        d_sum = float(np.sum(d))
        tic_ratio = d_sum / o_sum if abs(o_sum) > 1e-12 else 1.0

        snr_orig = self._safe_snr(base)
        snr_now = self._safe_snr(target)
        snr_gain = snr_now / snr_orig if snr_orig > 1e-12 else 1.0

        return (
            f"Range: {base_mz[0]:.4f} - {base_mz[-1]:.4f}\n"
            f"Correlation: {corr:.4f}\n"
            f"TIC ratio: {tic_ratio:.3f}\n"
            f"SNR orig: {snr_orig:.1f}\n"
            f"SNR den: {snr_now:.1f}\n"
            f"SNR improvement: {snr_gain:.2f}x"
        )

    @staticmethod
    def _save_or_show(fig, save_path: Optional[str], dpi: int) -> None:
        if save_path:
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()

    def add_metrics_box(self, ax, base: Spectrum, target: Spectrum) -> None:
        if len(base) <= 5 or len(target) <= 5:
            return

        text = self._metrics_text(base, target)
        logger.info(text)
        ax.text(
            0.02,
            0.98,
            text,
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
        )

    def plot_single(self, base: Spectrum, style: SeriesStyle) -> None:
        mz, intensity = self._slice_arrays(base, self.config.mz_range)
        fig, ax = plt.subplots(1, 1, figsize=self.config.figsize)

        self._draw(ax, mz, intensity, style)

        limits = self._axis_limits(mz, intensity)
        self._apply_axis_limits(ax, limits)
        ax.set_title("Mass Spectrum", fontweight="bold")
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")
        self._style_axes(ax)
        fig.tight_layout()

        self._save_or_show(fig, self.config.save_path, self.config.dpi)

    def plot_overlay(self, base: Spectrum, target: Spectrum, styles: Sequence[SeriesStyle]) -> None:
        base_mz, base_int = self._slice_arrays(base, self.config.mz_range)
        target_mz, target_int = self._slice_arrays(target, self.config.mz_range)

        fig, ax = plt.subplots(1, 1, figsize=self.config.figsize)
        self._draw(ax, base_mz, base_int, styles[0])
        self._draw(ax, target_mz, target_int, styles[1])

        limits = self._axis_limits(base_mz, base_int, target_mz, target_int)
        self._apply_axis_limits(ax, limits)
        ax.set_title("Original & Preprocessed (Overlay)", fontweight="bold")
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")
        self._style_axes(ax)
        ax.legend(loc="upper right")

        if self.config.metrics_box:
            self.add_metrics_box(ax, base, target)

        fig.tight_layout()
        self._save_or_show(fig, self.config.save_path, self.config.dpi)

    def plot_compare(self, base: Spectrum, target: Spectrum, styles: Sequence[SeriesStyle]) -> None:
        base_mz, base_int = self._slice_arrays(base, self.config.mz_range)
        target_mz, target_int = self._slice_arrays(target, self.config.mz_range)

        fig, (ax_top, ax_bottom) = plt.subplots(
            2,
            1,
            figsize=self.config.figsize,
            sharex=True,
            sharey=True,
        )

        self._draw(ax_top, base_mz, base_int, styles[0])
        self._draw(ax_bottom, target_mz, target_int, styles[1])

        limits = self._axis_limits(base_mz, base_int, target_mz, target_int)
        self._apply_axis_limits((ax_top, ax_bottom), limits)

        ax_top.set_title("Original Spectrum", fontweight="bold")
        ax_bottom.set_title("Preprocessed Spectrum", fontweight="bold")
        ax_bottom.set_xlabel("m/z")
        ax_top.set_ylabel("Intensity")
        ax_bottom.set_ylabel("Intensity")
        self._style_axes(ax_top)
        self._style_axes(ax_bottom)

        if self.config.metrics_box:
            self.add_metrics_box(ax_bottom, base, target)

        fig.tight_layout()
        self._save_or_show(fig, self.config.save_path, self.config.dpi)



def plot_spectrum(
    base: Optional[Spectrum] = None,
    target: Optional[Spectrum] = None,
    save_path=None,
    figsize=(20, 5),
    dpi: int = 300,
    colors: Optional[Sequence[str]] = None,
    line_width: Optional[Sequence[float]] = None,
    plot_mode: Optional[Sequence[str]] = None,
    mz_range: Optional[Tuple[float, float]] = None,
    intensity_range: Optional[Tuple[float, float]] = None,
    metrics_box: bool = True,
    overlay: bool = False,
):
    """Public API: plot one or two spectra with optional overlay and metrics."""
    if base is None and target is None:
        raise ValueError("At least one of 'base' or 'target' spectrum must be provided.")

    if base is None and target is not None:
        base = target
        target = None

    base = SpectrumPlotter._validate_spectrum("base", base)

    color_pair = SpectrumPlotter._normalize_pair(colors, SpectrumPlotter.DEFAULT_COLORS)
    line_pair = [float(v) for v in SpectrumPlotter._normalize_pair(line_width, [1.5, 1.5])]
    mode_pair = [SpectrumPlotter._normalize_mode(m) for m in SpectrumPlotter._normalize_pair(plot_mode, ["line", "line"])]

    resolved_figsize = figsize if (overlay or target is None) else (figsize[0], figsize[1] * 2)

    config = PlotConfig(
        save_path=save_path,
        figsize=(float(resolved_figsize[0]), float(resolved_figsize[1])),
        dpi=int(dpi),
        mz_range=mz_range,
        intensity_range=intensity_range,
        metrics_box=bool(metrics_box),
    )

    logger.info(
        f"Plotting spectrum: overlay={overlay}, "
        f"modes={mode_pair}, linewidth={line_pair}, mz_range={mz_range}, "
        f"intensity_range={intensity_range}, metrics_box={metrics_box}"
    )

    plotter = SpectrumPlotter(config)

    if target is None:
        style = SeriesStyle(color=str(color_pair[0]), line_width=line_pair[0], mode=mode_pair[0])
        plotter.plot_single(base, style)
        return

    target = SpectrumPlotter._validate_spectrum("target", target)
    styles = [
        SeriesStyle(
            color=str(color_pair[0]),
            line_width=line_pair[0],
            mode=mode_pair[0],
            label="Original",
        ),
        SeriesStyle(
            color=str(color_pair[1]),
            line_width=line_pair[1],
            mode=mode_pair[1],
            label="Preprocessed",
        ),
    ]

    if overlay:
        plotter.plot_overlay(base, target, styles)
    else:
        plotter.plot_compare(base, target, styles)