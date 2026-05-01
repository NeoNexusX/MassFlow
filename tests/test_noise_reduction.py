import time
import pytest
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.flat_pre_fun import FlatPreprocess
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
        for _, intensity_flat, lengths, _ in dm.flat_generator(batch_size=4096, include_mz=False, max_threads=16):
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
