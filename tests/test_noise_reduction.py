import time
import pytest
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import BatchPreprocess
from massflow.tools.dm_process import speed_process
from massflow.tools.logger import get_logger

logger = get_logger("massflow.test.test_noise_reduction")

ROUNDS = 2


class TestNoiseReductionAPI:
    """Test suite for noise reduction API functionality and performance."""

    @pytest.fixture(scope="module", params=["data/example.imzML"])
    def ms_raw_data(self, request) -> MSDataManagerImzML:
        """Fixture providing MSDataManagerImzML instance with fully initialized spectra for noise reduction tests."""
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()
        for _ in dm.batch_generator(batch_size=512):
            pass
        return dm

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", ["ma", "ma_numba"])
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
