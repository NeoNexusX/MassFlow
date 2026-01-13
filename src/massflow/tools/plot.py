from typing import Sequence, Tuple, Optional
import numpy as np
import matplotlib.pyplot as plt
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.module.spectrum import Spectrum
from massflow.tools.logger import get_logger

logger = get_logger("tools.plot")

def plot_spectrum(
    base: Optional["Spectrum"] = None,
    target: Optional["Spectrum"] = None,
    save_path=None,
    figsize=(20, 5),
    dpi: int = 300,
    colors: Optional[Sequence[str]] = None,
    line_width: Optional[Sequence[float]] = None,
    plot_mode: Optional[Sequence[str]] = None,
    mz_range: Optional[Tuple[float, float]] = None,
    intensity_range: Optional[Tuple[float, float]] = None,
    metrics_box: bool = True,
    title_suffix: Optional[str] = None,
    overlay: bool = False,
):
    """
    Plot one or two spectra with optional overlay and metrics box.

    Parameters:
        base (SpectrumBaseModule | None): Base/original spectrum to plot.
        target (SpectrumBaseModule | None): Target/processed spectrum to plot.
        save_path (str | None): File path to save the figure; display if None.
        figsize (Tuple[int,int]): Figure size.
        dpi (int): Saving DPI.
        colors (List[str]): Colors for base/target.
        plot_mode (List[str]): Per-spectrum mode: 'line' or 'stem'.
        mz_range (Tuple[float,float] | None): X-axis range.
        intensity_range (Tuple[float,float] | None): Y-axis range.
        metrics_box (bool): Whether to overlay metrics.
        title_suffix (str | None): Optional suffix added to title.
        overlay (bool): Plot spectra together if True, otherwise on separate axes.

    Returns:
        None
    """
    # initialize default values
    line_width = [1, 1] if line_width is None else line_width
    colors = ['#5c9dba', '#df4c5b'] if colors is None else colors
    plot_mode = ["line","line"] if plot_mode is None else plot_mode

    logger.info(f"Plotting spectrum with plot_mode={plot_mode},"
                "line_width={line_width}, mz_range={mz_range}, intensity_range={intensity_range},"
                "metrics_box={metrics_box}, title_suffix={title_suffix}, overlay={overlay}")

    figsize=figsize if overlay or target is None else (figsize[0], figsize[1] * 2)

    if base is None and target is None:
        logger.warning("No spectrum data provided for plotting.")
        raise ValueError("At least one of 'base' or 'target' spectrum must be provided.")
    elif target and base is None:
        logger.info("Only target spectrum provided; plotting single spectrum.")
        base = target
        target = None

    if target is None:
        plot_single(
            base=base,
            save_path=save_path,
            figsize=figsize,
            dpi=dpi,
            color=colors[0],
            line_width=line_width[0],
            plot_mode=plot_mode[0],
            mz_range=mz_range,
            intensity_range=intensity_range,
            title_suffix=title_suffix,
        )
    elif overlay:
        plot_two_together(
            target=target,
            base=base,
            save_path=save_path,
            figsize=figsize if overlay else (figsize[0], figsize[1] * 2),
            dpi=dpi,
            color=colors,
            plot_mode=plot_mode,
            line_width=line_width,
            mz_range=mz_range,
            intensity_range=intensity_range,
            metrics_box=metrics_box,
            title_suffix=title_suffix,
        )
    else:
        plot_two_individual(
            target=target,
            base=base,
            save_path=save_path,
            figsize=figsize,
            dpi=dpi,
            color=colors,
            plot_mode=plot_mode,
            line_width=line_width,
            mz_range=mz_range,
            intensity_range=intensity_range,
            metrics_box=metrics_box,
            title_suffix=title_suffix,
        )

def add_metrics_box(ax,
                    base,
                    target,
                    box_loc: Tuple[float, float] = (0.02, 0.98),
                    fontsize: int = 9):
    """
    Overlay a metrics box summarizing correlation, TIC ratio, and SNR.

    Parameters:
        ax (matplotlib.axes.Axes): Target axes.
        base (SpectrumBaseModule): Original spectrum.
        target (SpectrumBaseModule): Processed spectrum.
        box_loc (Tuple[float,float]): Axes-relative location of the box.
        fontsize (int): Font size.

    Returns:
        None
    """
    if base is None or target is None:
        return
    min_len = min(len(base), len(target))
    if min_len <= 1:
        return

    o = base.intensity[:min_len]
    d = target.intensity[:min_len]

    corr = float(np.corrcoef(o, d)[0, 1])
    tic_ratio = float(d.sum() / o.sum()) if o.sum() > 0 else 1.0

    snr_orig = SpectrumPreprocess.calculate_snr_spectrum(base)
    snr_update = SpectrumPreprocess.calculate_snr_spectrum(target)
    snr_improvement = snr_update / snr_orig if snr_orig > 0 else 1.0

    metrics_text = (f"Range: {base.mz_list[0]:.4f} - {base.mz_list[-1]:.4f}\n"
                    f"Correlation: {corr:.4f}\n"
                    f"TIC ratio: {tic_ratio:.3f}\n"
                    f"SNR orig: {snr_orig:.1f}\n"
                    f"SNR den: {snr_update:.1f}\n"
                    f"SNR improvement: {snr_improvement:.2f}x")
    # Log the metrics
    logger.info(metrics_text)
    ax.text(box_loc[0],
            box_loc[1],
            metrics_text,
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=fontsize)

def plot_single(base,
                save_path=None,
                figsize=(20, 5),
                dpi: int = 300,
                color='#5c9dba',
                line_width: float = 1,
                plot_mode: str = "line",
                mz_range: Optional[Tuple[float, float]] = None,
                intensity_range: Optional[Tuple[float, float]] = None,
                title_suffix: Optional[str] = None):
    """
    Plot a single spectrum in either line or stem mode.

    Parameters:
        base (SpectrumBaseModule): Spectrum to plot.
        save_path (str | None): File path to save; display if None.
        figsize (Tuple[int,int]): Figure size.
        dpi (int): Saving DPI.
        color (str): Plot color.
        plot_mode (str): 'line' or 'stem'.
        mz_range (Tuple[float,float] | None): X-axis range.
        intensity_range (Tuple[float,float] | None): Y-axis range.
        title_suffix (str | None): Title suffix.

    Returns:
        None
    """

    plt.figure(figsize=figsize)
    mode = (plot_mode or "stem").lower()
    if mode == "line":
        plt.plot(base.mz_list, base.intensity, color=color, linewidth=line_width, alpha=0.8)
    else:
        markerline, stemlines, baseline = plt.stem(base.mz_list, base.intensity)
        plt.setp(stemlines, linewidth=line_width, color=color, alpha=0.8)
        plt.setp(markerline, markersize=3, color=color, alpha=0.8)
        plt.setp(baseline, linewidth=0.5, color='gray', alpha=0.6)

    # Axis range control
    x_min, x_max = (float(min(base.mz_list)), float(max(base.mz_list))) if mz_range is None else (mz_range[0], mz_range[1])
    y_min, y_max = (0.0, float(max(base.intensity)) * 1.05) if intensity_range is None else (intensity_range[0], intensity_range[1])
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    title = "Mass Spectrum" if not title_suffix else f"Mass Spectrum - {title_suffix}"
    plt.title(title)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi)
    else:
        plt.show()
    return

def plot_two_together(target: Optional['Spectrum'] = None,
                      base: Optional['Spectrum'] = None,
                      save_path=None,
                      figsize=(20, 5),
                      dpi: int = 300,
                      color: Optional[Sequence[str]] = None,
                      line_width: Optional[Sequence[float]] = None,
                      plot_mode: Optional[Sequence[str]] = None,
                      mz_range: Optional[Tuple[float, float]] = None,
                      intensity_range: Optional[Tuple[float, float]] = None,
                      metrics_box: bool = True,
                      title_suffix: Optional[str] = None):
    """
    Overlay base and processed spectra on a single axes.

    Parameters:
        target (SpectrumBaseModule | None): Processed spectrum.
        base (SpectrumBaseModule | None): Original spectrum.
        save_path (str | None): File path to save; display if None.
        figsize (Tuple[int,int]): Figure size.
        dpi (int): Saving DPI.
        color (List[str]): Colors for base/target.
        plot_mode (List[str]): Modes for base/target: 'line' or 'stem'.
        mz_range (Tuple[float,float] | None): X-axis range.
        intensity_range (Tuple[float,float] | None): Y-axis range.
        metrics_box (bool): Whether to overlay metrics.
        title_suffix (str | None): Title suffix.

    Returns:
        None
    """

    _, ax = plt.subplots(1, 1, figsize=figsize if isinstance(figsize, tuple) else (12, 6))
    color = ['#5d7db3', '#d2c3d5'] if color is None else color
    line_width = [1, 1] if line_width is None else line_width
    plot_mode = ["line","line"] if plot_mode is None else plot_mode

    if plot_mode[0] == "line":
        ax.plot(base.mz_list, base.intensity, color=color[0], linewidth=line_width[0], label='Original')

    else:
        m1, s1, b1 = ax.stem(base.mz_list, base.intensity, label='Original')
        plt.setp(s1, linewidth=line_width[0], color=color[0], alpha=0.8)
        plt.setp(m1, markersize=3, color=color[0], alpha=0.8)
        plt.setp(b1, linewidth=0.5, color='gray', alpha=0.6)

    if plot_mode[1] == "line":
        ax.plot(target.mz_list, target.intensity, color=color[1], linewidth=line_width[1], label='Preprocessed' if not title_suffix else f'Preprocessed ({title_suffix})')
    else:
        m2, s2, b2 = ax.stem(target.mz_list, target.intensity, label='Preprocessed' if not title_suffix else f'Preprocessed ({title_suffix})')
        plt.setp(s2, linewidth=line_width[1], color=color[1], alpha=0.8)
        plt.setp(m2, markersize=3, color=color[1], alpha=0.8)
        plt.setp(b2, linewidth=0.5, color='gray', alpha=0.6)

    orig_c = base.crop_range(mz_range) 
    now_c = target.crop_range(mz_range)

    # Axis range settings (combined)
    x_min = float(mz_range[0]) if mz_range is not None else float(min(orig_c.mz_list.min(), now_c.mz_list.min()))
    x_max = float(mz_range[1]) if mz_range is not None else float(max(orig_c.mz_list.max(), now_c.mz_list.max()))
    y_max_comb = float(intensity_range[1]) if intensity_range is not None else float(max(orig_c.intensity.max(), now_c.intensity.max(), 1.0)) * 1.05
    y_min_comb = 0.0 if intensity_range is None else float(intensity_range[0])

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min_comb, y_max_comb)

    # Titles, labels, legend, grid
    ax.set_title('Original & Preprocessed (Overlay)' if not title_suffix else f'Original & Preprocessed (Overlay) - {title_suffix}', fontweight='bold')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Overlay metrics box on the same axes
    if metrics_box and len(orig_c) > 5 and len(now_c) > 5:
        add_metrics_box(ax, base=orig_c, target=now_c)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    else:
        plt.show()

def plot_two_individual(target: Spectrum,
                        base: Spectrum,
                        save_path=None,
                        figsize=(20, 5),
                        dpi: int = 300,
                        color: Optional[Sequence[str]] = None,
                        line_width: Optional[Sequence[float]] = None,
                        plot_mode: Optional[Sequence[str]] = None,
                        mz_range: Optional[Tuple[float, float]] = None,
                        intensity_range: Optional[Tuple[float, float]] = None,
                        metrics_box: bool = True,
                        title_suffix: Optional[str] = None):
    """
    Plot base and processed spectra on separate aligned subplots.

    Parameters:
        target (SpectrumBaseModule): Processed spectrum.
        base (SpectrumBaseModule): Original spectrum.
        save_path (str | None): File path to save; display if None.
        figsize (Tuple[int,int]): Figure size.
        dpi (int): Saving DPI.
        color (List[str]): Colors for base/target.
        plot_mode (List[str]): Modes for base/target: 'line' or 'stem'.
        mz_range (Tuple[float,float] | None): X-axis range.
        intensity_range (Tuple[float,float] | None): Y-axis range.
        metrics_box (bool): Whether to overlay metrics.
        title_suffix (str | None): Title suffix.

    Returns:
        None
    """
    # get two subplots
    _, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=figsize if isinstance(figsize, tuple) else (12, 8), sharex=True, sharey=True)
    color = ['#5d7db3', '#d2c3d5'] if color is None else color
    line_width = [1, 1] if line_width is None else line_width
    plot_mode = ["line","line"] if plot_mode is None else plot_mode

    if plot_mode[0] == "line":
        ax_top.plot(base.mz_list, base.intensity, color=color[0], linewidth=line_width[0])

    else:
        m1, s1, b1 = ax_top.stem(base.mz_list, base.intensity)
        plt.setp(s1, linewidth=line_width[0], color=color[0], alpha=0.8)
        plt.setp(m1, markersize=3, color=color[0], alpha=0.8)
        plt.setp(b1, linewidth=0.5, color='gray', alpha=0.6)

    if plot_mode[1] == "line":
        ax_bottom.plot(target.mz_list, target.intensity, color=color[1], linewidth=line_width[1])
    else:
        m2, s2, b2 = ax_bottom.stem(target.mz_list, target.intensity)
        plt.setp(s2, linewidth=line_width[1], color=color[1], alpha=0.8)
        plt.setp(m2, markersize=3, color=color[1], alpha=0.8)
        plt.setp(b2, linewidth=0.5, color='gray', alpha=0.6)

    # Axis range settings
    x_min = float(mz_range[0]) if mz_range is not None else float(min(base.mz_list.min(), target.mz_list.min()))
    x_max = float(mz_range[1]) if mz_range is not None else float(max(base.mz_list.max(), target.mz_list.max()))
    y_top = float(intensity_range[1]) if intensity_range is not None else float(max(base.intensity.max(), 1.0)) * 1.05
    y_bot = float(intensity_range[1]) if intensity_range is not None else float(max(target.intensity.max(), 1.0)) * 1.05

    ax_top.set_xlim(x_min, x_max)
    ax_bottom.set_xlim(x_min, x_max)
    ax_top.set_ylim(0.0 if intensity_range is None else float(intensity_range[0]), y_top)
    ax_bottom.set_ylim(0.0 if intensity_range is None else float(intensity_range[0]), y_bot)

    # Titles and grid
    ax_top.set_title('Original Spectrum', fontweight='bold')
    den_title = 'Preprocessed Spectrum' if not title_suffix else f'Preprocessed Spectrum ({title_suffix})'
    ax_bottom.set_title(den_title, fontweight='bold')
    ax_bottom.set_xlabel('m/z')
    ax_top.set_ylabel('Intensity')
    ax_bottom.set_ylabel('Intensity')
    ax_top.grid(True, alpha=0.3)
    ax_bottom.grid(True, alpha=0.3)


    base = base.crop_range(mz_range)
    target = target.crop_range(mz_range)

    # Overlay metrics text
    if metrics_box and len(base) > 5 and len(target) > 5:
        add_metrics_box(ax_bottom, base=base, target=target)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    else:
        plt.show()
