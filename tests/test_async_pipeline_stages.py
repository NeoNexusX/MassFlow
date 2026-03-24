import os
import tempfile
from queue import Queue
from threading import Event, Lock

import pytest

from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.preprocess.batch_pre_fun import BatchPreprocess
from massflow.preprocess.dm_pre_fun import Preprocess

pytestmark = pytest.mark.filterwarnings(
    r"ignore:This process .* is multi-threaded, use of fork\(\):DeprecationWarning"
)


def _make_temp_dir(prefix: str) -> str:
    os.makedirs("./temp_noise_async_data", exist_ok=True)
    return tempfile.mkdtemp(prefix=prefix, dir="./temp_noise_async_data")


@pytest.fixture(scope="module")
def async_data_manager() -> MSDataManagerImzML:
    mass_data = MassSpectrumSet()
    dm = MSDataManagerImzML(mass_data, filepath="data/example.imzML")
    dm.load_head_data()
    yield dm
    dm.close()


class TestAsyncPipelineStagesMemory:
    """Memory profile tests for each async pipeline stage."""

    def test_stage_a_memory_profile(self, async_data_manager: MSDataManagerImzML) -> None:
        """Stage A memory profile (run with: pytest --memray)."""
        pipeline = Preprocess.pipeline(
            data_manager=async_data_manager,
            batch_size=128,
            temp_dir=_make_temp_dir("stage_a_"),
            queue_ab_size=5,
            queue_bc_size=10,
            keep_order=False,
        )

        queue_ab = Queue(maxsize=5)
        stop_event = Event()
        error_holder: list[BaseException] = []
        error_lock = Lock()

        pipeline._stage_a_reader(
            queue_ab=queue_ab,
            stop_event=stop_event,
            error_holder=error_holder,
            error_lock=error_lock,
        )

        batch_count = 0
        while True:
            chunk = queue_ab.get()
            if chunk is None:
                break
            batch_count += 1
            async_data_manager.clear_batch_data_memory(chunk.batch)

        assert not error_holder
        assert batch_count > 0

    def test_stage_b_memory_profile(self, async_data_manager: MSDataManagerImzML) -> None:
        """Stage B memory profile (run with: pytest --memray)."""
        pipeline = (
            Preprocess.pipeline(
                data_manager=async_data_manager,
                batch_size=128,
                temp_dir=_make_temp_dir("stage_b_"),
                queue_ab_size=5,
                queue_bc_size=10,
                keep_order=False,
            )
            .noise_reduction(method="ma_loop", window=10, numba_max_threads=4)
        )

        ordered_tasks = pipeline._sorted_tasks()
        source_batch = next(async_data_manager.batch_generator(batch_size=128))

        queue_ab = Queue(maxsize=5)
        queue_bc = Queue(maxsize=10)
        stop_event = Event()
        error_holder: list[BaseException] = []
        error_lock = Lock()

        from massflow.preprocess.async_pipeline import BatchChunk

        queue_ab.put(BatchChunk(index=0, batch=source_batch))
        queue_ab.put(None)

        pipeline._stage_b_processor(
            ordered_tasks=ordered_tasks,
            queue_ab=queue_ab,
            queue_bc=queue_bc,
            stop_event=stop_event,
            error_holder=error_holder,
            error_lock=error_lock,
        )

        out_chunk = queue_bc.get()
        sentinel = queue_bc.get()

        assert sentinel is None
        assert out_chunk is not None
        assert len(out_chunk.batch) == len(source_batch)
        assert not error_holder

        async_data_manager.clear_batch_data_memory(source_batch)
        async_data_manager.clear_batch_data_memory(out_chunk.batch)

    def test_stage_c_memory_profile(self, async_data_manager: MSDataManagerImzML) -> None:
        """Stage C memory profile (run with: pytest --memray)."""
        pipeline = Preprocess.pipeline(
            data_manager=async_data_manager,
            batch_size=128,
            temp_dir=_make_temp_dir("stage_c_pipeline_"),
            queue_ab_size=5,
            queue_bc_size=10,
            keep_order=False,
        )

        input_batch = next(async_data_manager.batch_generator(batch_size=128))
        processed_batch = BatchPreprocess.noise_reduction_batch(
            batch_spectra=input_batch,
            method="ma_loop",
            window=10,
            numba_max_threads=4,
        )
        async_data_manager.clear_batch_data_memory(input_batch)

        processed_dm = MSDataManagerImzML(
            MassSpectrumSet(),
            temp_dir=_make_temp_dir("stage_c_output_"),
        )
        processed_dm.copy_meta(async_data_manager)
        writer = processed_dm.writer

        queue_bc = Queue(maxsize=10)
        stop_event = Event()
        error_holder: list[BaseException] = []
        error_lock = Lock()

        from massflow.preprocess.async_pipeline import BatchChunk

        queue_bc.put(BatchChunk(index=0, batch=list(processed_batch)))
        queue_bc.put(None)

        pipeline._stage_c_writer(
            queue_bc=queue_bc,
            stop_event=stop_event,
            error_holder=error_holder,
            error_lock=error_lock,
            processed_data_manager=processed_dm,
            writer=writer,
            total_batches=1,
        )

        processed_dm.close_writer()
        processed_dm.load_head_data()

        assert not error_holder
        assert len(processed_dm.ms) == len(processed_batch)

        processed_dm.close()
