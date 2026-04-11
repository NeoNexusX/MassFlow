from __future__ import annotations

from dataclasses import dataclass, field
from queue import Queue
from threading import Event, Lock, Thread
from typing import Any, Callable, Optional, cast

import numpy as np

from massflow.data_manager import MSDataManager, MSDataManagerImzML
from massflow.module import MassSpectrumSet
from massflow.preprocess.api import PreprocessorAPI, TaskScope
from massflow.preprocess.helper.peak_align_helper import reference_computer
from massflow.preprocess.flat_pre_fun import FlatPreprocess, FlatBatchResult
from massflow.preprocess.numba.numba_runtime import apply_numba_runtime
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.async_pipeline")

@dataclass(slots=True)
class FlatChunk:
    """Transport object between pipeline stages."""

    index: int
    mz_data: Optional[np.ndarray]
    intensity: np.ndarray
    lengths: np.ndarray
    coordinates: np.ndarray


@dataclass(slots=True)
class PreprocessTask:
    """Lazy task registration model for batch/dataset preprocessing."""

    name: str
    apply_fn: Callable[..., FlatBatchResult] | Callable[..., MSDataManagerImzML]
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
        numba_max_threads: Optional[int] = None,
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
        self.numba_max_threads = numba_max_threads

        self._tasks: list[PreprocessTask] = []
        self._task_sequence = 0

    def _register_task(
        self,
        name: str,
        *,
        scope: TaskScope = "batch",
        apply_fn: Callable[..., FlatBatchResult] | Callable[..., MSDataManagerImzML],
        **kwargs: Any,
    ) -> "Preprocessor":
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
        queue_ab: Queue[Optional[FlatChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
        """Stage A: producer; reads flat and pushes to queue_ab."""
        try:
            for batch_idx, (mz_data, intensity, lengths, coords) in enumerate(
                data_manager.flat_generator(batch_size=self.batch_size, include_mz=True),
                start=0
            ):
                if stop_event.is_set():
                    break
                queue_ab.put(
                    FlatChunk(
                        index=batch_idx,
                        mz_data=np.asarray(mz_data),
                        intensity=np.asarray(intensity),
                        lengths=np.asarray(lengths, dtype=np.int32),
                        coordinates=np.asarray(coords),
                    )
                )
        except BaseException as exc:  # pylint: disable=broad-exception-caught
            self._set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc)
        finally:
            queue_ab.put(None)

    def _stage_b_processor(
        self,
        *,
        tasks: list[PreprocessTask],
        queue_ab: Queue[Optional[FlatChunk]],
        queue_bc: Queue[Optional[FlatChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
        """Stage B: middleware; runs ordered preprocessing tasks for each batch."""
        try:
            apply_numba_runtime(override_workers=self.numba_max_threads)
            while True:
                chunk = queue_ab.get()
                if chunk is None:
                    queue_bc.put(None)
                    break

                if stop_event.is_set():
                    continue

                current_chunk = chunk

                for task in tasks:
                    batch_apply_fn = cast(Callable[..., FlatBatchResult], task.apply_fn)
                    output = batch_apply_fn(
                        mz_data=current_chunk.mz_data,
                        intensity=current_chunk.intensity,
                        lengths=current_chunk.lengths,
                        **task.kwargs,
                    )
                    current_chunk = FlatChunk(
                        index=current_chunk.index,
                        mz_data=output.mz_data,
                        intensity=output.intensity,
                        lengths=output.lengths,
                        coordinates=current_chunk.coordinates,
                    )

                queue_bc.put(current_chunk)

        except BaseException as exc:  # pylint: disable=broad-exception-caught
            self._set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc,)
            queue_bc.put(None)

    def _flush_chunk(
        self,
        *,
        chunk: FlatChunk,
        processed_data_manager: MSDataManagerImzML,
        written_batches: int,
        total_batches: int,
    ) -> int:
        processed_data_manager.swap_flat_data_out2disk(
            mz_flat=chunk.mz_data,
            intensity_flat=chunk.intensity,
            lengths=chunk.lengths,
            coordinates=chunk.coordinates,
        )
        written_batches += 1

        if total_batches <= 10 or written_batches % max(1, total_batches // 10) == 0 or written_batches == total_batches:
            progress = (written_batches / total_batches) * 100 if total_batches > 0 else 100.0
            logger.info(f"async_pipeline progress: {written_batches}/{total_batches} batches ({progress:.1f}%)")

        return written_batches

    def _flush_pending_in_order(
        self,
        *,
        pending: dict[int, FlatChunk],
        expected_index: int,
        processed_data_manager: MSDataManagerImzML,
        written_batches: int,
        total_batches: int,
    ) -> tuple[int, int]:
        while expected_index in pending:
            ordered_chunk = pending.pop(expected_index)
            written_batches = self._flush_chunk(
                chunk=ordered_chunk,
                processed_data_manager=processed_data_manager,
                written_batches=written_batches,
                total_batches=total_batches,
            )
            expected_index += 1
        return expected_index, written_batches

    def _ensure_no_pending_gap(
        self,
        *,
        pending: dict[int, FlatChunk],
        expected_index: int,
    ) -> None:
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
        queue_bc: Queue[Optional[FlatChunk]],
        stop_event: Event,
        error_holder: list[BaseException],
        error_lock: Lock,
        processed_data_manager: MSDataManagerImzML,
        total_batches: int,
    ) -> None:
        """Stage C: consumer; persist processed chunks, optionally keeping order."""
        written_batches = 0
        expected_index = 0
        pending: dict[int, FlatChunk] = {}

        try:
            while True:
                chunk = queue_bc.get()
                if chunk is None:
                    if self.keep_order:
                        expected_index, written_batches = self._flush_pending_in_order(
                            pending=pending,
                            expected_index=expected_index,
                            processed_data_manager=processed_data_manager,
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
                        written_batches=written_batches,
                        total_batches=total_batches,
                    )
                    continue

                pending[chunk.index] = chunk
                expected_index, written_batches = self._flush_pending_in_order(
                    pending=pending,
                    expected_index=expected_index,
                    processed_data_manager=processed_data_manager,
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

    def _run_flat_batch_task(
        self,
        *,
        data_manager: MSDataManager,
        tasks: list[PreprocessTask],
    ) -> MSDataManagerImzML:
        total_batches = (len(data_manager.ms) + self.batch_size - 1) // self.batch_size

        processed_ms = MassSpectrumSet()
        processed_data_manager = MSDataManagerImzML(processed_ms, temp_dir=self.temp_dir)
        processed_data_manager.copy_meta(data_manager)

        queue_ab: Queue[Optional[FlatChunk]] = Queue(maxsize=self.queue_ab_size)
        queue_bc: Queue[Optional[FlatChunk]] = Queue(maxsize=self.queue_bc_size)

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
            reference, tolerance = reference_computer(
                data_manager=data_manager,
                reference=reference,
                tolerance=tolerance,
                units=units,
                binfun=kwargs.get("binfun", "median"),
                binratio=kwargs.get("binratio", 2.0),
                batch_size=self.batch_size,
            )
        else:
            tolerance = tolerance * 1e-6 if units == "ppm" else tolerance

        kwargs["reference"] = reference
        kwargs["tolerance"] = tolerance

        kwargs.pop("binfun", None)
        kwargs.pop("binratio", None)

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
        if task.name == "peak_align" and task.apply_fn is FlatPreprocess.peak_align_flat:
            prepared_task = self._prepare_peak_align_task(
                data_manager=data_manager,
                task=task,
            )
            return self._run_flat_batch_task(
                data_manager=data_manager,
                tasks=[prepared_task],
            )

        if task.name == "peak_align" or task.name == "peak_pick":
            dataset_apply_fn = cast(Callable[..., MSDataManagerImzML], task.apply_fn)
            return dataset_apply_fn(data_manager=data_manager, **task.kwargs)
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

            # prepare consecutive batch tasks to run in the pipeline together, then run dataset task separately for better performance
            while task_cursor < len(ordered_tasks) and ordered_tasks[task_cursor].scope == "batch":
                batch_tasks.append(ordered_tasks[task_cursor])
                task_cursor += 1

            if batch_tasks:
                next_data_manager = self._run_flat_batch_task(
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

        return cast(MSDataManagerImzML, current_data_manager)
