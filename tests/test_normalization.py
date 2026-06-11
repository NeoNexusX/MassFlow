import time
from pathlib import Path

import numpy as np
import pytest
import matplotlib

matplotlib.use("Agg")

from massflow.data_manager import MSDataManagerImzML
from massflow.module import Spectrum
from massflow.preprocess import BatchPreprocess
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.tools.logger import get_logger
from massflow.tools.plot import plot_spectrum

logger = get_logger("test_normalization")

ROUNDS = 5
BATCH_NORM_METHODS = ["tic", "rms"]
BATCH_NORM_METRIC_METHODS = ["tic", "rms", "rl_norm"]
FLAT_NORM_METHODS = ["tic_numba", "rms_numba", "ref_numba"]
FLAT_SCALE_NORM_METHODS = ["tic_numba", "rms_numba"]


def _normalization_flat_from_flat_batches(
    flat_batches,
    method: str,
    scale: float | None = None,
    ref_tolerance: float = 0.1,
):
    for mz_flat, intensity_flat, lengths, ref in flat_batches:
        kwargs = {
            "mz_data": mz_flat,
            "intensity": intensity_flat,
            "method": method,
            "scale": scale,
            "lengths": lengths,
        }
        if method == "ref_numba":
            kwargs["ref"] = ref
            kwargs["ref_tolerance"] = ref_tolerance

        _ = FlatPreprocess.normalization_flat(**kwargs)


def _make_metric_test_batch(rows: int = 3, cols: int = 3, n_ion: int = 48) -> list[Spectrum]:
    rng = np.random.default_rng(13)
    base = rng.lognormal(mean=1.6, sigma=0.25, size=n_ion)
    mz = np.linspace(100.0, 900.0, n_ion, dtype=float)

    spectra: list[Spectrum] = []
    for y in range(rows):
        for x in range(cols):
            field = float(np.exp(0.22 * np.sin((x + 1) / 2.0) - 0.18 * np.cos((y + 1) / 2.1)))
            intensity = np.asarray(base * field, dtype=float)
            spectra.append(
                Spectrum(
                    mz_list=mz,
                    intensity=intensity,
                    coordinate=(x + 1, y + 1, 1),
                )
            )
    return spectra


class TestNormalizationAPI:
    """Normalization API tests: memory, speed, consistency, and normalization invariants."""

    @pytest.fixture(scope="module", params=["data/example.imzML"])
    def ms_raw_data(self, request) -> MSDataManagerImzML:
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()
        for _ in dm.batch_generator(batch_size=512):
            pass
        return dm

    @pytest.fixture(scope="module", params=["data/example.imzML"])
    def flat_caches(self, request):
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()

        caches = []
        for mz_data, intensity_flat, lengths, _ in dm.flat_generator(
            batch_size=4096,
            include_mz=True,
            max_threads=16,
        ):
            ref = float(mz_data[mz_data.size // 2]) if mz_data is not None and mz_data.size > 0 else None
            caches.append((mz_data, intensity_flat, lengths, ref))

        return caches

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", FLAT_NORM_METHODS)
    def test_norm_flat_speed(
        self,
        benchmark,
        method,
        flat_caches,
    ):
        logger.info("Benchmarking flat normalization method=%s", method)

        benchmark.pedantic(
            _normalization_flat_from_flat_batches,
            args=(flat_caches, method),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

    @pytest.mark.parametrize("method", FLAT_SCALE_NORM_METHODS)
    def test_norm_invariants_default_equals_length(self, method, flat_caches):
        _, intensity_flat, lengths, _ = flat_caches[0]
        intensity_flat = np.asarray(intensity_flat, dtype=np.float64)

        normalized_flat = FlatPreprocess.normalization_flat(
            mz_data=None, # type: ignore
            intensity=intensity_flat,
            method=method,
            scale=None,
            lengths=lengths,
        )

        offset = 0
        for valid_len in lengths:
            end = offset + int(valid_len)
            y = normalized_flat.intensity[offset:end]
            if method == "tic_numba":
                assert np.isclose(np.sum(y), float(valid_len), rtol=1e-6)
            else:
                assert np.isclose(np.sqrt(np.mean(y ** 2)), float(valid_len), rtol=1e-6)
            offset = end

    @pytest.mark.parametrize("method", FLAT_SCALE_NORM_METHODS)
    def test_norm_invariants_scale_equals_param(self, method, flat_caches):
        _, intensity_flat, lengths, _ = flat_caches[0]
        intensity_flat = np.asarray(intensity_flat, dtype=np.float64)

        scale_target = 3.5
        normalized_flat = FlatPreprocess.normalization_flat(
            mz_data=None, # type: ignore
            intensity=intensity_flat,
            method=method,
            scale=scale_target,
            lengths=lengths,
        )

        offset = 0
        for valid_len in lengths:
            end = offset + int(valid_len)
            y = normalized_flat.intensity[offset:end]
            if method == "tic_numba":
                assert np.isclose(np.sum(y), scale_target, rtol=1e-6)
            else:
                assert np.isclose(np.sqrt(np.mean(y ** 2)), scale_target, rtol=1e-6)
            offset = end

    def test_norm_ref_peak_scaled_to_target(self, flat_caches):
        mz_flat, intensity_flat, lengths, ref = flat_caches[0]
        assert mz_flat is not None
        assert ref is not None

        is_shared_mz = (
            mz_flat.ndim == 1
            and mz_flat.size != intensity_flat.size
            and lengths.size > 0
            and mz_flat.size == int(lengths[0])
        )

        scale_target = 3.5
        normalized_flat = FlatPreprocess.normalization_flat(
            mz_data=mz_flat,
            intensity=intensity_flat,
            method="ref_numba",
            scale=scale_target,
            ref=ref,
            ref_tolerance=0.1,
            lengths=lengths,
        )

        offset = 0
        for valid_len in lengths:
            end = offset + int(valid_len)
            mz_slice = mz_flat if is_shared_mz else mz_flat[offset:end]
            y = normalized_flat.intensity[offset:end]

            idx = int(np.argmin(np.abs(mz_slice - ref)))
            if np.abs(mz_slice[idx] - ref) <= 0.1 and intensity_flat[offset + idx] > 0.0:
                logger.info(
                       f"scaling to target={scale_target},y[idx]={y[idx]}"
                )
                assert np.isclose(y[idx], scale_target, rtol=1e-6, atol=1e-6)

            offset = end

    @pytest.mark.parametrize("method", BATCH_NORM_METRIC_METHODS)
    def test_norm_batch_metrics_effective(self, method):
        """Use quantitative metrics to verify each batch normalization method is effective."""
        batch = _make_metric_test_batch()

        before_tic = np.array([float(np.asarray(s.intensity, dtype=float).sum()) for s in batch], dtype=float)

        normalized = BatchPreprocess.normalization_batch(
            batch_spectra=batch,
            method=method,
            scale=1.0,
        )

        assert len(normalized) == len(batch)

        after_tic = []
        for spectrum in normalized:
            y = np.asarray(spectrum.intensity, dtype=float)
            assert y.ndim == 1
            assert y.size > 0
            assert np.all(np.isfinite(y))
            assert float(np.linalg.norm(y)) > 0.0
            after_tic.append(float(np.asarray(y, dtype=float).sum()))

        after_tic_arr = np.asarray(after_tic, dtype=float)

        if method == "tic":
            assert np.allclose(after_tic_arr, 1.0, rtol=1e-6, atol=1e-6)
        elif method == "rms":
            for spectrum in normalized:
                y = np.asarray(spectrum.intensity, dtype=float)
                assert np.isclose(float(np.sqrt(np.mean(y ** 2))), 1.0, rtol=1e-6, atol=1e-6)
        elif method == "rl_norm":
            before_cv = float(np.std(before_tic) / (np.mean(before_tic) + 1e-12))
            after_cv = float(np.std(after_tic_arr) / (np.mean(after_tic_arr) + 1e-12))
            logger.info("RL_NORM_CV before=%f after=%f", before_cv, after_cv)
            assert after_cv < before_cv

    def test_norm_batch_rl_norm_scale_metric(self):
        """RL-Norm should preserve linear scaling behavior."""
        batch = _make_metric_test_batch()

        y1 = BatchPreprocess.normalization_batch(
            batch_spectra=batch,
            method="rl_norm",
            scale=1.0,
        )
        y2 = BatchPreprocess.normalization_batch(
            batch_spectra=batch,
            method="rl_norm",
            scale=2.0,
        )

        for s1, s2 in zip(y1, y2):
            np.testing.assert_allclose(
                np.asarray(s2.intensity),
                2.0 * np.asarray(s1.intensity),
                rtol=1e-6,
                atol=1e-6,
            )

    def test_norm_batch_rl_norm_visualization_saved(self):
        batch = _make_metric_test_batch()

        normalized = BatchPreprocess.normalization_batch(
            batch_spectra=batch,
            method="rl_norm",
            scale=1.0,
        )

        before_cv = float(
            np.std([float(np.asarray(s.intensity, dtype=float).sum()) for s in batch])
            / (np.mean([float(np.asarray(s.intensity, dtype=float).sum()) for s in batch]) + 1e-12)
        )
        after_cv = float(
            np.std([float(np.asarray(s.intensity, dtype=float).sum()) for s in normalized])
            / (np.mean([float(np.asarray(s.intensity, dtype=float).sum()) for s in normalized]) + 1e-12)
        )

        output_path = Path("temp") / "plot_examples" / "normalization" / "batch_rl_norm_before_after.png"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plot_spectrum(
            base=batch[0],
            target=normalized[0],
            save_path=str(output_path),
            overlay=True,
            metrics_box=True,
        )

        logger.info("Saved normalization example plot to %s", output_path)

        assert output_path.exists()
        assert output_path.stat().st_size > 0
        assert after_cv < before_cv
