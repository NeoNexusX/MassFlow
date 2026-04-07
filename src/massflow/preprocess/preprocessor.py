from __future__ import annotations

from dataclasses import dataclass, field
from queue import Queue
from threading import Event, Lock, Thread
from typing import Any, Callable, Optional, Sequence

from massflow.data_manager import MSDataManager, MSDataManagerImzML
from massflow.module import MassSpectrumSet, Spectrum
from massflow.preprocess.api import PreprocessorAPI, TaskScope
from massflow.preprocess.helper.peak_align_helper import compute_reference
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
    """Lazy task registration model for batch/dataset preprocessing."""

    name: str
    apply_fn: Callable[..., Sequence[Spectrum]] | Callable[..., MSDataManagerImzML]
    scope: TaskScope = "batch"
    kwargs: dict[str, Any] = field(default_factory=dict)
    sequence: int = 0

class Preprocessor(PreprocessorAPI):
    """
    Hybrid preprocessing executor.

    - batch tasks are executed by 3-stage async pipeline
    - dataset tasks are executed on whole data manager (for peak align and cardinal backends)
    """

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
        queue_ab_size: int = 1,
        queue_bc_size: int = 1,
        temp_dir: str | None = "./temp",
        keep_order: bool = False,
    ): # pylint: disable=super-init-not-called
        if data_manager is None or batch_size <= 0 or batch_size > 9056 or queue_ab_size <= 0 or queue_bc_size <= 0 :
            logger.error(f"Invalid parameter values. please check:"
                         f"data_manager={data_manager}, batch_size={batch_size},"
                         f"queuesize={queue_ab_size} & {queue_bc_size}? <0?")
            raise ValueError("Invalid parameter values.")

        self.data_manager = data_manager
        self.batch_size = batch_size
        self.temp_dir = temp_dir
        self.queue_ab_size = queue_ab_size
        self.queue_bc_size = queue_bc_size
        self.keep_order = keep_order

        self._tasks: list[PreprocessTask] = []
        self._task_sequence = 0

    def _register_task(
        self,
        name: str,
        *,
        scope: TaskScope = "batch",
        apply_fn: Callable[..., Sequence[Spectrum]] | Callable[..., MSDataManagerImzML],
        **kwargs: Any,
    ) -> "Preprocessor":
        if scope == "batch" and apply_fn is None:
            raise ValueError(f"Batch task {name} requires apply_fn.")

        self._task_sequence += 1
        self._tasks.append(
            PreprocessTask(
                name=name,
                scope=scope,
                apply_fn=apply_fn,
                kwargs=kwargs,
                sequence=self._task_sequence,
            )
        )
        return self

    def _sorted_tasks(self) -> list[PreprocessTask]:
        def sort_key(task: PreprocessTask) -> tuple[int, int]:
            return (self._OPERATION_ORDER.get(task.name, 10_000), task.sequence)

        ordered = sorted(self._tasks, key=sort_key)

        if ordered:
            task_names = [f"{task.name}:{task.scope}" for task in ordered]
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
        data_manager: MSDataManager,
        queue_ab: Queue[Optional[BatchChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
        """Stage A: producer; reads batch and pushes to queue_ab."""
        try:
            for batch_idx, batch in enumerate(data_manager.batch_generator(batch_size=self.batch_size), start=0):
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
        tasks: list[PreprocessTask],
        queue_ab: Queue[Optional[BatchChunk]],
        queue_bc: Queue[Optional[BatchChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
        """Stage B: middleware; runs ordered preprocessing tasks for each batch."""
        try:
            while True:
                chunk = queue_ab.get()
                if chunk is None:
                    queue_bc.put(None)
                    break

                if stop_event.is_set():
                    continue

                current_batch = chunk.batch

                for task in tasks:
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
            "Pipeline ended with missing chunks in keep_order mode. "
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

    def _run_batch_task(
        self,
        *,
        data_manager: MSDataManager,
        tasks: list[PreprocessTask],
    ) -> MSDataManagerImzML:
        total_batches = (len(data_manager.ms) + self.batch_size - 1) // self.batch_size

        processed_ms = MassSpectrumSet()
        processed_data_manager = MSDataManagerImzML(processed_ms, temp_dir=self.temp_dir)
        processed_data_manager.copy_meta(data_manager)
        writer = processed_data_manager.writer

        queue_ab: Queue[Optional[BatchChunk]] = Queue(maxsize=self.queue_ab_size)
        queue_bc: Queue[Optional[BatchChunk]] = Queue(maxsize=self.queue_bc_size)

        stop_event = Event()
        error_holder: list[BaseException] = []
        error_lock = Lock()

        t_a = Thread(
            target=self._stage_a_reader,
            kwargs={
                "data_manager": data_manager,
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
                "tasks": tasks,
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
            processed_data_manager.close()
            raise RuntimeError("Async preprocessing pipeline failed.") from error_holder[0]

        processed_data_manager.close_writer()
        processed_data_manager.load_head_data()
        return processed_data_manager

    def _prepare_peak_align_task(
        self,
        *,
        data_manager: MSDataManager,
        task: PreprocessTask,
    ) -> PreprocessTask:
        """Resolve python peak align parameters and build a batch task."""
        kwargs = dict(task.kwargs)

        units = kwargs.get("units", "ppm")
        reference = kwargs.get("reference")
        tolerance = kwargs.get("tolerance")

        if reference is None or tolerance is None:
            reference, tol_internal = compute_reference(
                data_manager=data_manager,
                reference=reference,
                tolerance=tolerance,
                units=units,
                binfun=kwargs.get("binfun", "median"),
                binratio=kwargs.get("binratio", 2.0),
                clear_memory=kwargs.get("clear_memory", False),
                batch_size=self.batch_size,
            )

            kwargs["reference"] = reference
            kwargs["tolerance"] = tol_internal * 1e6 if units == "ppm" else tol_internal

        kwargs.pop("binfun", None)
        kwargs.pop("binratio", None)
        kwargs.pop("clear_memory", None)

        return PreprocessTask(
            name=task.name,
            scope="batch",
            apply_fn=task.apply_fn,
            kwargs=kwargs,
            sequence=task.sequence,
        )

    def _run_dataset_task(
        self,
        *,
        data_manager: MSDataManager,
        task: PreprocessTask,
    ) -> MSDataManagerImzML:
        """Execute whole-dataset task."""
        if task.name == "peak_align" and task.apply_fn is BatchPreprocess.peak_align_batch:
            prepared_task = self._prepare_peak_align_task(
                data_manager=data_manager,
                task=task,
            )
            return self._run_batch_task(
                data_manager=data_manager,
                tasks=[prepared_task],
            )

        if task.name == "peak_align" or task.name == "peak_pick":
            return task.apply_fn(data_manager=data_manager, **task.kwargs)
        else:
            raise NotImplementedError(f"Dataset task {task.name} with backend Cardinal is not supported.")

    def start(self) -> MSDataManagerImzML:
        """Execute registered tasks and return new data manager with processed data. """
        if len(self._tasks) == 0:
            raise ValueError("No preprocessing task registered. Please register tasks before start().")

        ordered_tasks = self._sorted_tasks()
        current_data_manager: MSDataManager = self.data_manager
        task_cursor = 0

        while task_cursor < len(ordered_tasks):
            batch_tasks: list[PreprocessTask] = []

            while task_cursor < len(ordered_tasks) and ordered_tasks[task_cursor].scope == "batch":
                batch_tasks.append(ordered_tasks[task_cursor])
                task_cursor += 1

            if batch_tasks:
                next_data_manager = self._run_batch_task(
                    data_manager=current_data_manager,
                    tasks=batch_tasks,
                )

                if current_data_manager is not self.data_manager:
                    current_data_manager.close()

                current_data_manager = next_data_manager

            if task_cursor >= len(ordered_tasks):
                break

            dataset_task = ordered_tasks[task_cursor]
            task_cursor += 1

            next_data_manager = self._run_dataset_task(
                data_manager=current_data_manager,
                task=dataset_task,
            )

            if current_data_manager is not self.data_manager:
                current_data_manager.close()

            current_data_manager = next_data_manager

        return current_data_manager
