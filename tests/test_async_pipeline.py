import time
import os
import subprocess
import pytest
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import Preprocessor
from massflow.tools.logger import get_logger

logger = get_logger("test_async_pipeline")
ROUND = 1


def drop_caches():
    if os.name == "nt":
        return
    subprocess.run(["sudo", "purge"], check=True)


def run_pipeline_once(batch_size, queue_ab_size, queue_bc_size):
    with MSDataManagerImzML(filepath="data/example.imzML") as dm:
        dm.load_head_data()
        processed_data_manager = (
            Preprocessor(
                data_manager=dm,
                batch_size=batch_size,
                queue_ab_size=queue_ab_size,
                queue_bc_size=queue_bc_size,
            )
            .noise_reduction(method="ma_numba", window=5)
            .start()
        )
        processed_data_manager.close()

class TestAsyncPipelineProfile:
    """Memory and runtime profile tests for async pipeline and each internal stage."""

    @pytest.mark.parametrize("batch_size", [1024])
    @pytest.mark.parametrize("queue_ab_size", [2])
    @pytest.mark.parametrize("queue_bc_size", [2])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_async_pipeline_benchmark(
        self,
        benchmark,
        batch_size: int,
        queue_ab_size: int,
        queue_bc_size: int,
    ) -> None:
        """Whole pipeline runtime benchmark."""
        benchmark.pedantic(
            run_pipeline_once,
            args=(batch_size, queue_ab_size, queue_bc_size),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        drop_caches()
