import time
import pytest

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import BatchPreprocess
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.tools.dm_process import speed_process
from massflow.tools.logger import get_logger

logger = get_logger("test_baseline")

ROUNDS = 5
BATCH_BASELINE_METHODS = ["locmin", "snip"]
FLAT_BASELINE_METHODS = ["locmin_numba", "snip_numba"]


def _baseline_reduction_flat_from_flat_batches(
    flat_batches,
    method: str,
    width: int,
):
    for intensity_flat, lengths in flat_batches:
        _ = FlatPreprocess.baseline_reduction_flat(
            intensity=intensity_flat,
            method=method,
            width=width,
            lengths=lengths,
        )


class TestBaseline:
    """
    Baseline correction benchmark tests.
            use:
            uv run pytest ./tests/test_baseline.py -k "test_baseline_speed or test_baseline_flat_speed" -q
    """

    @pytest.fixture(scope="module", params=["data/baseline_test.ibd"])
    def ms_raw_data(self, request) -> MSDataManagerImzML:
        """Fixture providing batch-readable data manager cache for baseline benchmarks."""
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()
        for _ in dm.batch_generator(batch_size=512):
            pass
        return dm

    @pytest.fixture(scope="module", params=["data/baseline_test.ibd"])
    def flat_caches(self, request):
        """Fixture providing pre-generated flat arrays and lengths for flat benchmarks."""
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()

        caches = []
        for _, intensity_flat, lengths, _ in dm.flat_generator(
            batch_size=4096, include_mz=False, max_threads=16
        ):
            caches.append((intensity_flat, lengths))

        return caches

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", BATCH_BASELINE_METHODS)
    def test_baseline_speed(self, benchmark, method, ms_raw_data):
        """Benchmark batch baseline correction via speed_process."""
        logger.info(f"Benchmarking batch baseline correction method={method}")

        batch_kwargs = {
            "method": method,
            "width": 5,
            "smooth": "none",
        }

        benchmark.pedantic(
            speed_process,
            args=(
                ms_raw_data,
                1024,
                BatchPreprocess.baseline_correction_batch,
                batch_kwargs,
            ),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("method", FLAT_BASELINE_METHODS)
    def test_baseline_flat_speed(self, benchmark, method, flat_caches):
        """Benchmark flat numba baseline reduction via baseline_reduction_flat."""
        logger.info(f"Benchmarking flat baseline reduction method={method}")

        flat_kwargs = {
            "method": method,
            "width": 5,
        }

        benchmark.pedantic(
            _baseline_reduction_flat_from_flat_batches,
            args=(flat_caches, flat_kwargs["method"], flat_kwargs["width"]),
            rounds=ROUNDS,
            iterations=1,
            warmup_rounds=1,
        )
