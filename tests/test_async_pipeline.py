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


def run_pipeline_once(
    data_manager,
    batch_size=256,
    queue_ab_size=2,
    queue_bc_size=2
):
    """Run the async pipeline once with specified parameters."""
    processed_data_manager = (
        Preprocessor(
            data_manager=data_manager,
            batch_size=batch_size,
            queue_ab_size=queue_ab_size,
            queue_bc_size=queue_bc_size,
        )
        .noise_reduction()
        .peak_pick()
        .peak_align()
        .start()
    )
    processed_data_manager.close()


class TestAsyncPipelineProfile:
    """Memory and runtime profile tests for async pipeline and each internal stage."""

    @pytest.fixture(scope="module", params=[r"C:\Users\dre\Desktop\data\280TopL, 250TopR, 270BottomL, 260BottomR-profile\280TopL, 250TopR, 270BottomL, 260BottomR-profile.imzML"])
    def ms_raw_data(self, request) -> MSDataManagerImzML:
        """Fixture providing MSDataManagerImzML instance with fully initialized spectra for noise reduction tests."""
        data_file_path = request.param
        dm = MSDataManagerImzML(filepath=data_file_path)
        dm.load_head_data()
        for _ in dm.batch_generator(batch_size=512):
            pass
        return dm

    @pytest.mark.parametrize("batch_size", [1024])
    @pytest.mark.parametrize("queue_ab_size", [2])
    @pytest.mark.parametrize("queue_bc_size", [2])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_async_pipeline_benchmark(
        self,
        ms_raw_data,
        benchmark,
        batch_size: int,
        queue_ab_size: int,
        queue_bc_size: int,
    ) -> None:
        """Whole pipeline runtime benchmark."""
        benchmark.pedantic(
            run_pipeline_once,
            args=(ms_raw_data, batch_size, queue_ab_size, queue_bc_size),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
        )
        drop_caches()

    @pytest.mark.parametrize(("profile", "centroid"), [(None, None), (True, True)])
    def test_pipeline_resolve_spectrum_type(
        self,
        ms_raw_data,
        profile,
        centroid,
    ):
        """Test that the pipeline correctly resolves spectrum type."""
        ms_raw_data.ms.meta.profile_spectrum = profile
        ms_raw_data.ms.meta.centroid_spectrum = centroid

        processed_data_manager = (
            Preprocessor(
                data_manager=ms_raw_data,
                batch_size=256,
                queue_ab_size=2,
                queue_bc_size=2,
            )
            .noise_reduction()
            .peak_pick()
            .peak_align()
            .start()
        )
        try:
            assert processed_data_manager.ms.meta.profile_spectrum is None
            assert processed_data_manager.ms.meta.centroid_spectrum is True
        finally:
            processed_data_manager.close()

    def test_pipeline_resolve_tasks(
        self,
        ms_raw_data,
    ):
        """Test that the pipeline correctly resolves tasks based on spectrum type."""
        ms_raw_data.ms.meta.profile_spectrum = None
        ms_raw_data.ms.meta.centroid_spectrum = True

        preprocessor = (
            Preprocessor(
                data_manager=ms_raw_data,
                batch_size=256,
                queue_ab_size=2,
                queue_bc_size=2,
            )
            .noise_reduction()
            .peak_pick()
            .peak_align()
        )
        ordered_task_names = [task.name for task in preprocessor._sorted_tasks("centroid")] # pylint: disable=protected-access
        assert ordered_task_names == ["peak_align"]

        unsupported_preprocessor = Preprocessor(data_manager=ms_raw_data).noise_reduction()
        with pytest.raises(ValueError, match="No supported preprocessing tasks remain"):
            unsupported_preprocessor.start()

        processed_data_manager = preprocessor.start()
        processed_data_manager.close()
