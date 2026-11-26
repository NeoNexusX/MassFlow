import time
import numpy as np
import pytest
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.preprocess.numba.noise_reduction_numba import smooth_signal_savgol_numba
from massflow.tools.logger import get_logger

logger = get_logger("massflow.test.test_noise_reduction")

ROUNDS = 5
FLAT_NR_METHODS = ["ma_numba", "gaussian_numba", "savgol_numba"]


def _noise_reduction_flat_from_flat_batches(
    flat_batches,
    method: str,
    window: int,
):
    for intensity_flat, lengths in flat_batches:
        _ = FlatPreprocess.noise_reduction_flat(
            mz_data=None, # type: ignore
            intensity=intensity_flat,
            method=method,
            window=window,
            lengths=lengths,
        )


class TestNoiseReductionAPI:
    """
    Test suite for noise reduction API functionality and performance.
        use :
        uv run pytest ./tests/test_noise_reduction.py -k "test_nr_flat_speed" -q
    """

    @pytest.fixture(scope="module", params=["data/example.imzML"])
    def flat_caches(self, request):
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()

        caches = []
        for _, intensity_flat, lengths, _ in dm.flat_generator(batch_size=4096, include_mz=False):
            caches.append((intensity_flat, lengths))

        return caches

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", FLAT_NR_METHODS)
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


def test_savgol_short_spectrum_uses_edge_padding():
    first_spectrum = np.array([1.0, 2.0, 3.0])
    lengths = np.array([3, 3], dtype=np.int64)

    result_a = smooth_signal_savgol_numba(
        np.concatenate((first_spectrum, np.zeros(3))),
        window=5,
        polyorder=3,
        lengths=lengths,
    )
    result_b = smooth_signal_savgol_numba(
        np.concatenate((first_spectrum, np.full(3, 1000.0))),
        window=5,
        polyorder=3,
        lengths=lengths,
    )

    expected = np.array([41.0 / 35.0, 2.0, 99.0 / 35.0])
    np.testing.assert_allclose(result_a[:3], expected)
    np.testing.assert_allclose(result_b[:3], expected)
