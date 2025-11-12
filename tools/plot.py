from typing import List, Tuple, Optional
import numpy as np
import matplotlib.pyplot as plt
from preprocess.ms_preprocess import MSIPreprocessor
from module.ms_module import SpectrumBaseModule
from logger import get_logger

logger = get_logger("tools.plot")

def plot_spectrum(
    base: Optional["SpectrumBaseModule"] = None,
    target:Optional["SpectrumBaseModule"] = None,
    save_path=None,
    figsize=(20, 5),
    dpi: int = 300,
    colors: List[str] = ['#5c9dba', '#df4c5b'],
    plot_mode: str = "line",
    mz_range: Optional[Tuple[float, float]] = None,
    intensity_range: Optional[Tuple[float, float]] = None,
    metrics_box: bool = True,
    title_suffix: Optional[str] = None,
    overlay: bool = False,
):
    logger.info(f"Plotting spectrum with plot_mode={plot_mode}, mz_range={mz_range}, intensity_range={intensity_range}, metrics_box={metrics_box}, title_suffix={title_suffix}, overlay={overlay}")
    figsize=figsize if overlay or target is None else (figsize[0], figsize[1] * 2)

    if target is None:
        plot_single(
            base=base,
            save_path=save_path,
            figsize=figsize,
            dpi=dpi,
            color=colors[0],
            plot_mode=plot_mode,
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
    Overlay metrics text box on the given axes based on original and processed intensity arrays.
    Automatically aligns to the shortest length when arrays differ.

    Metrics:
    - Correlation coefficient
    - TIC ratio (Total Ion Current ratio)
    - SNR orig / SNR den (MAD-based SNR estimate)
    - SNR improvement multiplier
    """
    if base is None or target is None:
        return
    min_len = min(len(base), len(target))
    if min_len <= 1:
        return

    o = base.intensity[:min_len]
    d = target.intensity[:min_len]

    corr = float(np.corrcoef(o, d)[0, 1]) if min_len > 1 else 0.0
    tic_ratio = float(d.sum() / o.sum()) if o.sum() > 0 else 1.0

    snr_orig = MSIPreprocessor.calculate_snr_spectrum(base)
    snr_update = MSIPreprocessor.calculate_snr_spectrum(target)
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
                color='#5d7db3',
                plot_mode: str = "line",
                mz_range: Optional[Tuple[float, float]] = None,
                intensity_range: Optional[Tuple[float, float]] = None,        
                title_suffix: Optional[str] = None):

    plt.figure(figsize=figsize)
    mode = (plot_mode or "stem").lower()
    if mode == "line":
        plt.plot(base.mz_list, base.intensity, color=color, linewidth=0.8, alpha=0.8)
    else:
        markerline, stemlines, baseline = plt.stem(base.mz_list, base.intensity)
        plt.setp(stemlines, linewidth=0.7, color=color, alpha=0.7)
        plt.setp(markerline, markersize=3, color=color, alpha=0.7)
        plt.setp(baseline, linewidth=0.5, color='gray', alpha=0.4)

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

def plot_two_together(target: Optional['SpectrumBaseModule'] = None,
                      base: Optional['SpectrumBaseModule'] = None,
                      save_path=None,
                      figsize=(20, 5),
                      dpi: int = 300,
                      color: List[str] = ['#5d7db3', '#d2c3d5'],
                      plot_mode: str = "line",
                      mz_range: Optional[Tuple[float, float]] = None,
                      intensity_range: Optional[Tuple[float, float]] = None,
                      metrics_box: bool = True,
                      title_suffix: Optional[str] = None):

    _, ax = plt.subplots(1, 1, figsize=figsize if isinstance(figsize, tuple) else (12, 6))

    if plot_mode == "line":
        ax.plot(base.mz_list, base.intensity, color=color[0], linewidth=1, label='Original')
        ax.plot(target.mz_list, target.intensity, color=color[1], linewidth=1, label='Preprocessed' if not title_suffix else f'Preprocessed ({title_suffix})')

    else:
        m1, s1, b1 = ax.stem(base.mz_list, base.intensity, label='Original')
        plt.setp(s1, linewidth=0.7, color=color[0], alpha=0.7)
        plt.setp(m1, markersize=3, color=color[0], alpha=0.7)
        plt.setp(b1, linewidth=0.5, color='gray', alpha=0.4)

        m2, s2, b2 = ax.stem(target.mz_list, target.intensity, label='Preprocessed' if not title_suffix else f'Preprocessed ({title_suffix})')
        plt.setp(s2, linewidth=0.7, color=color[1], alpha=0.7)
        plt.setp(m2, markersize=3, color=color[1], alpha=0.7)
        plt.setp(b2, linewidth=0.5, color='gray', alpha=0.4)

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

def plot_two_individual(target: SpectrumBaseModule,
                        base: SpectrumBaseModule,
                        save_path=None,
                        figsize=(20, 5),
                        dpi: int = 300,
                        color: List[str] = ['#5d7db3', '#d2c3d5'],
                        plot_mode: str = "line",
                        mz_range: Optional[Tuple[float, float]] = None,
                        intensity_range: Optional[Tuple[float, float]] = None,
                        metrics_box: bool = True,
                        title_suffix: Optional[str] = None):
    # get two subplots
    _, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=figsize if isinstance(figsize, tuple) else (12, 8), sharex=True, sharey=True)

    if plot_mode == "line":
        ax_top.plot(base.mz_list, base.intensity, color=color[0], linewidth=1)
        ax_bottom.plot(target.mz_list, target.intensity, color=color[1], linewidth=1)

    else:
        m1, s1, b1 = ax_top.stem(base.mz_list, base.intensity)
        plt.setp(s1, linewidth=0.7, color=color[0], alpha=0.7)
        plt.setp(m1, markersize=3, color=color[0], alpha=0.7)
        plt.setp(b1, linewidth=0.5, color='gray', alpha=0.4)

        m2, s2, b2 = ax_bottom.stem(target.mz_list, target.intensity)
        plt.setp(s2, linewidth=0.7, color=color[1], alpha=0.7)
        plt.setp(m2, markersize=3, color=color[1], alpha=0.7)
        plt.setp(b2, linewidth=0.5, color='gray', alpha=0.4)

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
