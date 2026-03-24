from __future__ import annotations

from dataclasses import dataclass, field
from queue import Queue
from threading import Event, Lock, Thread
from typing import Any, Callable, Optional, Sequence

from massflow.data_manager.ms_data_manager import MSDataManager
from massflow.data_manager.ms_data_manager_imzml import MSDataManagerImzML
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.spectrum import Spectrum
from massflow.preprocess.batch_pre_fun import BatchPreprocess
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.async_pipeline")


@dataclass(slots=True)
class BatchChunk:
    """Transport object between pipeline stages."""

    index: int
    batch: list[Spectrum]


@dataclass(slots=True)
class PreprocessTask:
    """Lazy task registration model for batch preprocessing."""

    name: str
    apply_fn: Callable[..., Sequence[Spectrum]]
    kwargs: dict[str, Any] = field(default_factory=dict)
    sequence: int = 0


class AsyncPreprocessPipeline:
    """Three-stage asynchronous preprocessing pipeline with bounded queues."""

    _OPERATION_ORDER: dict[str, int] = {
        "baseline_correction": 10,
        "noise_reduction": 20,
        "normalization": 30,
        "peak_pick": 40,
        "peak_align": 50,
    }

    def __init__(
        self,
        data_manager: MSDataManager,
        *,
        batch_size: int = 256,
        temp_dir: str = "./temp_pipeline_data",
        queue_ab_size: int = 3,
        queue_bc_size: int = 3,
        keep_order: bool = False,
    ):
        if data_manager is None:
            raise ValueError("data_manager must be provided for async preprocess pipeline.")
        if batch_size <= 0:
            raise ValueError("batch_size must be a positive integer.")
        if queue_ab_size <= 0 or queue_bc_size <= 0:
            raise ValueError("queue_ab_size and queue_bc_size must be positive integers.")

        self.data_manager = data_manager
        self.batch_size = batch_size
        self.temp_dir = temp_dir
        self.queue_ab_size = queue_ab_size
        self.queue_bc_size = queue_bc_size
        self.keep_order = keep_order

        self._tasks: list[PreprocessTask] = []
        self._task_sequence = 0

    def _register_task(self, name: str, apply_fn: Callable[..., Sequence[Spectrum]], **kwargs) -> "AsyncPreprocessPipeline":
        self._task_sequence += 1
        self._tasks.append(
            PreprocessTask(
                name=name,
                apply_fn=apply_fn,
                kwargs=kwargs,
                sequence=self._task_sequence,
            )
        )
        return self

    def baseline_correction(self, **kwargs) -> "AsyncPreprocessPipeline":
        return self._register_task("baseline_correction", BatchPreprocess.baseline_correction_batch, **kwargs)

    def noise_reduction(self, **kwargs) -> "AsyncPreprocessPipeline":
        return self._register_task("noise_reduction", BatchPreprocess.noise_reduction_batch, **kwargs)

    def normalization(self, **kwargs) -> "AsyncPreprocessPipeline":
        return self._register_task("normalization", BatchPreprocess.normalization_batch, **kwargs)

    def peak_align(self, **kwargs) -> "AsyncPreprocessPipeline":
        return self._register_task("peak_align", BatchPreprocess.peak_align_batch, **kwargs)

    def peak_pick(self, **kwargs) -> "AsyncPreprocessPipeline":
        return self._register_task("peak_pick", BatchPreprocess.peak_pick_batch, **kwargs)

    def _sorted_tasks(self) -> list[PreprocessTask]:
        def sort_key(task: PreprocessTask) -> tuple[int, int]:
            return (self._OPERATION_ORDER.get(task.name, 10_000), task.sequence)

        ordered = sorted(self._tasks, key=sort_key)

        if ordered:
            task_names = [task.name for task in ordered]
            logger.info(f"async_pipeline_task_order: {task_names}")

        return ordered

    def _set_error(
        self,
        *,
        error_holder: list[BaseException],
        error_lock: Lock,
        stop_event: Event,
        exc: BaseException,
    ) -> None:
        with error_lock:
            if not error_holder:
                error_holder.append(exc)
        stop_event.set()

    def _stage_a_reader(
        self,
        *,
        queue_ab: Queue[Optional[BatchChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
        """Stage A: producer; reads batch and pushes to queue_ab."""
        try:
            for batch_idx, batch in enumerate(self.data_manager.batch_generator(batch_size=self.batch_size), start=0):
                if stop_event.is_set():
                    break
                queue_ab.put(BatchChunk(index=batch_idx, batch=batch))
        except BaseException as exc:  # pylint: disable=broad-exception-caught
            self._set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc)
        finally:
            queue_ab.put(None)

    def _stage_b_processor(
        self,
        *,
        ordered_tasks: list[PreprocessTask],
        queue_ab: Queue[Optional[BatchChunk]],
        queue_bc: Queue[Optional[BatchChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
        """Stage B: middleware; runs ordered preprocessing tasks."""
        try:
            while True:
                chunk = queue_ab.get()
                if chunk is None:
                    queue_bc.put(None)
                    break

                if stop_event.is_set():
                    continue

                current_batch = chunk.batch

                for task in ordered_tasks:
                    next_batch = list(task.apply_fn(batch_spectra=current_batch, **task.kwargs))
                    self.data_manager.clear_batch_data_memory(current_batch)
                    current_batch = next_batch

                chunk.batch = current_batch
                queue_bc.put(chunk)

        except BaseException as exc:  # pylint: disable=broad-exception-caught
            self._set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc,)
            queue_bc.put(None)

    def _flush_chunk(
        self,
        *,
        chunk: BatchChunk,
        processed_data_manager: MSDataManagerImzML,
        writer,
        written_batches: int,
        total_batches: int,
    ) -> int:
        processed_data_manager.swap_batch_data_out2disk(batch=chunk.batch, writer=writer)
        processed_data_manager.clear_batch_data_memory(batch=chunk.batch)
        written_batches += 1

        if total_batches <= 10 or written_batches % max(1, total_batches // 10) == 0 or written_batches == total_batches:
            progress = (written_batches / total_batches) * 100 if total_batches > 0 else 100.0
            logger.info(f"async_pipeline progress: {written_batches}/{total_batches} batches ({progress:.1f}%)")

        return written_batches

    def _flush_pending_in_order(
        self,
        *,
        pending: dict[int, BatchChunk],
        expected_index: int,
        processed_data_manager: MSDataManagerImzML,
        writer,
        written_batches: int,
        total_batches: int,
    ) -> tuple[int, int]:
        while expected_index in pending:
            ordered_chunk = pending.pop(expected_index)
            written_batches = self._flush_chunk(
                chunk=ordered_chunk,
                processed_data_manager=processed_data_manager,
                writer=writer,
                written_batches=written_batches,
                total_batches=total_batches,
            )
            expected_index += 1
        return expected_index, written_batches

    def _ensure_no_pending_gap(self, *, pending: dict[int, BatchChunk], expected_index: int) -> None:
        if not pending:
            return

        missing = [idx for idx in range(expected_index, max(pending.keys()) + 1) if idx not in pending]
        raise RuntimeError(
            f"Pipeline ended with missing chunks in keep_order mode. "
            f"expected_index={expected_index}, pending={sorted(pending.keys())}, missing={missing}"
        )

    def _stage_c_writer(
        self,
        *,
        queue_bc: Queue[Optional[BatchChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
        processed_data_manager: MSDataManagerImzML,
        writer,
        total_batches: int,
    ) -> None:
        """Stage C: consumer; persist processed chunks, optionally keeping order."""
        written_batches = 0
        expected_index = 0
        pending: dict[int, BatchChunk] = {}

        try:
            while True:
                chunk = queue_bc.get()
                if chunk is None:
                    if self.keep_order:
                        expected_index, written_batches = self._flush_pending_in_order(
                            pending=pending,
                            expected_index=expected_index,
                            processed_data_manager=processed_data_manager,
                            writer=writer,
                            written_batches=written_batches,
                            total_batches=total_batches,
                        )
                        self._ensure_no_pending_gap(pending=pending, expected_index=expected_index)
                    break

                if stop_event.is_set():
                    continue

                if not self.keep_order:
                    written_batches = self._flush_chunk(
                        chunk=chunk,
                        processed_data_manager=processed_data_manager,
                        writer=writer,
                        written_batches=written_batches,
                        total_batches=total_batches,
                    )
                    continue

                pending[chunk.index] = chunk
                expected_index, written_batches = self._flush_pending_in_order(
                    pending=pending,
                    expected_index=expected_index,
                    processed_data_manager=processed_data_manager,
                    writer=writer,
                    written_batches=written_batches,
                    total_batches=total_batches,
                )

        except BaseException as exc:  # pylint: disable=broad-exception-caught
            self._set_error(
                error_holder=error_holder,
                error_lock=error_lock,
                stop_event=stop_event,
                exc=exc,
            )

    def start(self) -> MSDataManagerImzML:
        if len(self._tasks) == 0:
            raise ValueError("No preprocessing task registered. Please register tasks before start().")

        ordered_tasks = self._sorted_tasks()
        total_batches = (len(self.data_manager.ms) + self.batch_size - 1) // self.batch_size

        processed_ms = MassSpectrumSet()
        processed_data_manager = MSDataManagerImzML(processed_ms, temp_dir=self.temp_dir)
        processed_data_manager.copy_meta(self.data_manager)
        writer = processed_data_manager.writer

        queue_ab: Queue[Optional[BatchChunk]] = Queue(maxsize=self.queue_ab_size)
        queue_bc: Queue[Optional[BatchChunk]] = Queue(maxsize=self.queue_bc_size)

        stop_event = Event()
        error_holder: list[BaseException] = []
        error_lock = Lock()

        t_a = Thread(
            target=self._stage_a_reader,
            kwargs={
                "queue_ab": queue_ab,
                "stop_event": stop_event,
                "error_holder": error_holder,
                "error_lock": error_lock,
            },
            name="massflow-stage-a-reader",
            daemon=True,
        )
        t_b = Thread(
            target=self._stage_b_processor,
            kwargs={
                "ordered_tasks": ordered_tasks,
                "queue_ab": queue_ab,
                "queue_bc": queue_bc,
                "stop_event": stop_event,
                "error_holder": error_holder,
                "error_lock": error_lock,
            },
            name="massflow-stage-b-processor",
            daemon=True,
        )
        t_c = Thread(
            target=self._stage_c_writer,
            kwargs={
                "queue_bc": queue_bc,
                "stop_event": stop_event,
                "error_holder": error_holder,
                "error_lock": error_lock,
                "processed_data_manager": processed_data_manager,
                "writer": writer,
                "total_batches": total_batches,
            },
            name="massflow-stage-c-writer",
            daemon=True,
        )

        t_a.start()
        t_b.start()
        t_c.start()

        t_a.join()
        t_b.join()
        t_c.join()

        if error_holder:
            processed_data_manager.close_writer()
            raise RuntimeError("Async preprocessing pipeline failed.") from error_holder[0]

        processed_data_manager.close_writer()
        processed_data_manager.load_head_data()
        return processed_data_manager
