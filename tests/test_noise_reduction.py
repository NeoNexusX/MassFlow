import time
import pytest
import numpy as np
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import BatchPreprocess
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.tools.dm_process import speed_process
from massflow.tools.logger import get_logger

logger = get_logger("massflow.test.test_noise_reduction")

ROUNDS = 5


def _noise_reduction_flat_from_flat_batches(
    flat_batches,
    method: str,
    window: int,
):
    for intensity_flat, lengths in flat_batches:
        _ = FlatPreprocess.noise_reduction_flat(
            intensity=intensity_flat,
            method=method,
            window=window,
            lengths=lengths,
        )


class TestNoiseReductionAPI:
    """
    Test suite for noise reduction API functionality and performance.
        use :
        uv run  pytest ./tests/test_noise_reduction.py -k "test_nr_speed or test_nr_flat_speed" -q
    """

    @pytest.fixture(scope="module", params=["data/example.imzML"])
    def ms_raw_data(self, request) -> MSDataManagerImzML:
        """Fixture providing MSDataManagerImzML instance with fully initialized spectra for noise reduction tests."""
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
        for _, intensity_flat, lengths, _ in dm.flat_generator(batch_size=4096, include_mz=False, max_threads=16):
            caches.append((intensity_flat, lengths))

        return caches

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", ["ma_numba", "gaussian_numba", "savgol_numba"])
    def test_nr_speed(self,
                      benchmark,
                      method,
                      ms_raw_data):
        """Test noise reduction speed for different methods."""
        logger.info(f"Benchmarking noise reduction method={method}")

        batch_kwargs = {
            "method": method,
            "window": 5,

        }

        benchmark.pedantic(
            speed_process,
            args=(ms_raw_data, 4096, BatchPreprocess.noise_reduction_batch, batch_kwargs,8),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", ["ma_numba", "gaussian_numba", "savgol_numba"])
    def test_nr_flat_speed(self,
                           benchmark,
                           method,
                           flat_caches):
        """Test flat noise reduction speed using pre-generated flat batches."""
        logger.info(f"Benchmarking flat noise reduction method={method}")

        flat_kwargs = {
            "method": method,
            "window": 5,
        }

        benchmark.pedantic(
            _noise_reduction_flat_from_flat_batches,
            args=(flat_caches, flat_kwargs["method"], flat_kwargs["window"]),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

    @pytest.mark.parametrize("method", ["ma_numba", "gaussian_numba", "savgol_numba"])
    def test_nr_flat_batch_intensity_consistency(self, method, ms_raw_data):
        """Test flat and batch noise reduction intensity consistency under identical parameters."""
        batch = next(ms_raw_data.batch_generator(batch_size=256))

        batch_result = BatchPreprocess.noise_reduction_batch(
            batch_spectra=batch,
            method=method,
            window=5,
        )

        lengths = np.array([spectrum.intensity.size for spectrum in batch], dtype=np.int64)
        intensity_flat = np.concatenate(
            [spectrum.intensity.astype(np.float32, copy=False) for spectrum in batch]
        )
        flat_result = FlatPreprocess.noise_reduction_flat(
            intensity=intensity_flat,
            method=method,
            window=5,
            lengths=lengths,
        )

        offset = 0
        for spectrum_batch, valid_len in zip(batch_result, lengths):
            end = offset + int(valid_len)
            flat_slice = flat_result[offset:end]
            assert spectrum_batch.intensity is not None
            np.testing.assert_allclose(
                spectrum_batch.intensity,
                flat_slice,
                rtol=1e-6,
                atol=1e-6,
            )
            offset = end
