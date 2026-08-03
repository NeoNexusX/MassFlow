import main as main_module


class _DummyInputManager:
    def __init__(self, filepath):
        self.filepath = filepath
        self.loaded = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def load_head_data(self):
        self.loaded = True


class _DummyProcessedManager:
    n_ions = 12
    n_pixels = 34
    height = 5
    width = 6

    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True


def test_main_runs_ion_image_h5_pipeline(tmp_path, monkeypatch, capsys):
    calls = {}
    processed_manager = _DummyProcessedManager()

    class _DummyPreprocessor:
        def __init__(self, data_manager, output_format, output_path):
            calls["data_manager"] = data_manager
            calls["output_format"] = output_format
            calls["output_path"] = output_path

        def start(self):
            calls["started"] = True
            return processed_manager

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "example.h5").write_bytes(b"test")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(main_module, "MSDataManagerImzML", _DummyInputManager)
    monkeypatch.setattr(main_module, "Preprocessor", _DummyPreprocessor)

    main_module.main()

    input_manager = calls["data_manager"]
    assert input_manager.filepath == "data/example.imzML"
    assert input_manager.loaded is True
    assert calls["output_format"] == "ion_image_h5"
    assert calls["output_path"] == "data/example.h5"
    assert calls["started"] is True
    assert processed_manager.closed is True

    output = capsys.readouterr().out
    assert "Output: data/example.h5" in output
    assert "n_ions: 12" in output
    assert "n_pixels: 34" in output
    assert "image: 5 x 6" in output
