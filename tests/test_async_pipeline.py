import time
import pytest
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.tools.logger import get_logger

logger = get_logger("test_async_pipeline")
ROUND = 10


@pytest.fixture(scope="function")
def async_data_manager() -> MSDataManagerImzML:  # type: ignore[return]
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath="data/example.imzML")
    dm.load_head_data()
    yield dm  # type: ignore[return]
    dm.close()


class TestAsyncPipelineProfile:
    """Memory and runtime profile tests for async pipeline and each internal stage."""

    @pytest.mark.parametrize("batch_size", [128, 512, 1024])
    @pytest.mark.parametrize("queue_ab_size", [10])
    @pytest.mark.parametrize("queue_bc_size", [10])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_async_pipeline_benchmark(
        self,
        benchmark,
        async_data_manager: MSDataManagerImzML,
        batch_size: int,
        queue_ab_size: int,
        queue_bc_size: int,
    ) -> None:
        """Whole pipeline runtime benchmark."""

        def run_pipeline_once() -> int:
            processed_data_manager = (
                Preprocess.pipeline(
                    data_manager=async_data_manager,
                    batch_size=batch_size,
                    queue_ab_size=queue_ab_size,
                    queue_bc_size=queue_bc_size,
                    keep_order=False,
                )
                .noise_reduction(method="ma_loop", window=10, numba_max_threads=20)
                .start()
            )
            try:
                return len(processed_data_manager.ms)
            finally:
                processed_data_manager.close()

        result_len = benchmark.pedantic(
            run_pipeline_once,
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        logger.info(
            f"async pipeline benchmark result spectra={result_len}, "
            f"batch_size={batch_size}, queue_ab_size={queue_ab_size}, queue_bc_size={queue_bc_size}"
        )
        assert result_len == len(async_data_manager.ms)
