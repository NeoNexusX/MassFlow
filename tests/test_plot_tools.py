import numpy as np

import matplotlib

matplotlib.use("Agg")

from massflow.module import Spectrum
from massflow.tools.plot import plot_spectrum



def _build_spectrum(scale: float = 1.0) -> Spectrum:
    mz = np.linspace(100.0, 1000.0, 200)
    intensity = np.sin(np.linspace(0, 10, 200)) * 10 + 50
    intensity = (intensity * scale).astype(float)
    return Spectrum(mz_list=mz, intensity=intensity, coordinate=[1, 1, 1])



def test_plot_spectrum_requires_input():
    try:
        plot_spectrum(base=None, target=None)
    except ValueError as exc:
        assert "At least one" in str(exc)
        return

    assert False, "Expected ValueError when both base and target are None"



def test_plot_spectrum_single_save(tmp_path):
    out_file = tmp_path / "single_plot.png"
    base = _build_spectrum()

    plot_spectrum(
        base=base,
        save_path=str(out_file),
        overlay=False,
        metrics_box=False,
    )

    assert out_file.exists()
    assert out_file.stat().st_size > 0



def test_plot_spectrum_overlay_save(tmp_path):
    out_file = tmp_path / "overlay_plot.png"
    base = _build_spectrum(1.0)
    target = _build_spectrum(1.2)

    plot_spectrum(
        base=base,
        target=target,
        save_path=str(out_file),
        overlay=True,
        metrics_box=True,
        mz_range=(150.0, 900.0),
    )

    assert out_file.exists()
    assert out_file.stat().st_size > 0



def test_plot_spectrum_invalid_mz_range_raises(tmp_path):
    out_file = tmp_path / "invalid_plot.png"
    base = _build_spectrum()

    try:
        plot_spectrum(
            base=base,
            save_path=str(out_file),
            overlay=False,
            mz_range=(900.0, 150.0),
        )
    except ValueError as exc:
        assert "min_mz" in str(exc)
        return

    assert False, "Expected ValueError for invalid mz_range"
