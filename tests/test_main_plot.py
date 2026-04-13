import numpy as np

import main as main_module


class _DummySpectrum:
    def __init__(self, scale: float = 1.0):
        self.mz_list = np.array([100.0, 101.0, 102.0], dtype=float)
        self.intensity = np.array([10.0, 20.0, 30.0], dtype=float) * scale


def test_plot_before_after_calls_plot_with_expected_options(tmp_path, monkeypatch):
    calls = {}

    def _fake_plot_spectrum(**kwargs):
        calls.update(kwargs)

    monkeypatch.setattr(main_module, "plot_spectrum", _fake_plot_spectrum)

    before = _DummySpectrum(scale=1.0)
    after = _DummySpectrum(scale=1.5)
    out_file = tmp_path / "plots" / "before_after.png"

    main_module.plot_before_after(before, after, save_path=str(out_file))

    assert out_file.parent.exists()
    assert calls["base"] is before
    assert calls["target"] is after
    assert calls["save_path"] == str(out_file)
    assert calls["overlay"] is False
    assert calls["metrics_box"] is True


class _DummyProcessedResult:
    def __init__(self, intensity):
        self.intensity = intensity


class _DummyManager:
    def __init__(self, filepath=None, temp_dir=None):
        self.filepath = filepath
        self.temp_dir = temp_dir
        self.ms = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def load_head_data(self):
        if not self.ms:
            self.ms = [_DummySpectrum(1.0), _DummySpectrum(1.2)]

    def copy_meta(self, _other):
        return None

    def flat_generator(self, **_kwargs):
        mz_flat = np.array([100.0, 101.0, 102.0], dtype=float)
        intensity_flat = np.array([10.0, 20.0, 30.0], dtype=float)
        lengths = [3]
        coordinates = [[1, 1, 1]]
        yield mz_flat, intensity_flat, lengths, coordinates

    def swap_flat_data_out2disk(self, mz_flat, intensity_flat, lengths, coordinates):
        self.ms = []
        start = 0
        for length in lengths:
            end = start + length
            spectrum = _DummySpectrum(1.0)
            spectrum.mz_list = np.array(mz_flat[start:end], dtype=float)
            spectrum.intensity = np.array(intensity_flat[start:end], dtype=float)
            self.ms.append(spectrum)
            start = end

    def close_writer(self):
        return None

    def close(self):
        return None


def test_main_triggers_plot_before_after(monkeypatch):
    plot_calls = []

    def _fake_noise_reduction_flat(mz_data, intensity, lengths, method):
        assert method == "gaussian_numba"
        return _DummyProcessedResult(intensity=intensity * 0.8)

    def _fake_snr_details(_spectrum):
        return 1.0, 1.0, 1.0

    def _fake_plot_before_after(before_spectrum, after_spectrum, save_path=None):
        plot_calls.append((before_spectrum, after_spectrum, save_path))

    monkeypatch.setattr(main_module, "MSDataManagerImzML", _DummyManager)
    monkeypatch.setattr(main_module, "MassSpectrumSet", lambda: object())
    monkeypatch.setattr(main_module.FlatPreprocess, "noise_reduction_flat", staticmethod(_fake_noise_reduction_flat))
    monkeypatch.setattr(main_module, "calculate_snr_details", _fake_snr_details)
    monkeypatch.setattr(main_module, "log_snr_details", lambda *args, **kwargs: None)
    monkeypatch.setattr(main_module, "plot_before_after", _fake_plot_before_after)

    main_module.main()

    assert len(plot_calls) == 1
    before, after, save_path = plot_calls[0]
    assert isinstance(before, _DummySpectrum)
    assert isinstance(after, _DummySpectrum)
    assert save_path is None
