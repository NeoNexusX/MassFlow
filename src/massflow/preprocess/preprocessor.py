from __future__ import annotations
import tempfile
from dataclasses import dataclass, field
from queue import Queue
from threading import Event, Lock, Thread
from typing import Any, Callable, Literal, Optional, cast

import numpy as np

from massflow.data_manager import IonImageDataManagerBinary, IonImageDataManagerH5, MSDataManager, MSDataManagerImzML
from massflow.module import MassSpectrumSet
from massflow.preprocess.api import PreprocessorAPI, TaskScope
from massflow.preprocess.helper.peak_align_helper import reference_computer
from massflow.preprocess.flat_pre_fun import FlatPreprocess, FlatBatchResult
from massflow.preprocess.numba.numba_runtime import apply_numba_runtime
from massflow.tools.infer_spectrum_type import (
    SpectrumType,
    resolve_spectrum_type,
    set_spectrum_type_metadata,
)
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.async_pipeline")

@dataclass(slots=True)
class FlatChunk:
    """Transport object between pipeline stages."""

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
    requires_dataset_prepare: bool = False
    kwargs: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class StageCSink:
    """
    Unified Stage C sink wrapping the target data manager.
    `kind="imzml"` writes intermediate swap files; `kind="ion_image"` writes
    the final ion image output (binary or H5) and requires every chunk to
    carry mz_data.
    """

    data_manager: MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5
    kind: Literal["imzml", "ion_image"]

    def write_chunk(self, chunk: FlatChunk) -> None:
        if self.kind == "ion_image" and chunk.mz_data is None:
            raise ValueError("ion_image output requires mz_data from flat_generator.")
        self.data_manager.write_flat_batch(
            intensity_flat=chunk.intensity,
            lengths=chunk.lengths,
            coordinates=chunk.coordinates,
            mz_axis=chunk.mz_data,
        )

    def finalize(self) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5:
        if self.kind == "imzml":
            self.data_manager.close_writer()
            self.data_manager.load_head_data()
        else:
            self.data_manager.finalize()
        return self.data_manager

    def close(self) -> None:
        if self.kind == "imzml":
            self.data_manager.close_writer()
        self.data_manager.close()


class Preprocessor(PreprocessorAPI):
    """
    Hybrid preprocessing executor.

    - batch tasks are executed by 3-stage async pipeline
    - dataset tasks are executed on whole data manager (for peak align and cardinal backends)
    """

    _PROFILE_ONLY_TASKS = frozenset({"baseline_correction", "noise_reduction", "peak_pick"})
    _ION_IMAGE: dict[str, tuple[type, str]] = {
        "ion_image_binary": (IonImageDataManagerBinary, ".iib"),
        "ion_image_h5": (IonImageDataManagerH5, ".h5"),
    }

    def __init__(
        self,
        data_manager: MSDataManager,
        *,
        batch_size: int = 256,
        queue_ab_size: int = 1,
        queue_bc_size: int = 1,
        temp_dir: str | None = "./temp",
        numba_max_threads: Optional[int] = None,
        output_format: str = "imzml",
        output_path: str | None = None,
        output_max_workers: int = 2,
    ): # pylint: disable=super-init-not-called
        if data_manager is None or batch_size <= 0 or batch_size > 9056 or queue_ab_size <= 0 or queue_bc_size <= 0 :
            logger.error(f"Invalid parameter values. please check:"
                         f"data_manager={data_manager}, batch_size={batch_size},"
                         f"queuesize={queue_ab_size} & {queue_bc_size}? <0?")
            raise ValueError("Invalid parameter values.")

        if output_format not in {"imzml", "ion_image_binary", "ion_image_h5"}:
            raise ValueError("output_format must be 'imzml', 'ion_image_binary', or 'ion_image_h5'.")
        self.data_manager = data_manager
        self.batch_size = batch_size
        self.temp_dir = temp_dir
        self.queue_ab_size = queue_ab_size
        self.queue_bc_size = queue_bc_size
        self.numba_max_threads = numba_max_threads
        self.output_format = output_format
        self.output_path = output_path
        self.output_max_workers = output_max_workers

        self._tasks: list[PreprocessTask] = []

    def _register_task(
        self,
        name: str,
        *,
        scope: TaskScope = "batch",
        apply_fn: Callable[..., FlatBatchResult] | Callable[..., MSDataManagerImzML],
        requires_dataset_prepare: bool = False,
        **kwargs: Any,
    ) -> "Preprocessor":
        self._tasks.append(
            PreprocessTask(
                name=name,
                scope=scope,
                apply_fn=apply_fn,
                requires_dataset_prepare=requires_dataset_prepare,
                kwargs=kwargs,
            )
        )
        return self

    def _arrange_tasks(self, spectrum_type) -> list[PreprocessTask]:
        tasks = list(self._tasks)

        if spectrum_type == "centroid":
            dropped_task_names = [task.name for task in tasks if task.name in self._PROFILE_ONLY_TASKS]
            if dropped_task_names:
                logger.warning(f"Detected centroid spectra; dropping unsupported tasks: {dropped_task_names}")
                tasks = [task for task in tasks if task.name not in self._PROFILE_ONLY_TASKS]

        if tasks:
            task_names = [
                f"{task.name}:{task.scope}{':prepare' if task.requires_dataset_prepare else ''}"
                for task in tasks
            ]
            logger.info(f"async_pipeline_task_order: {task_names}")
        else:
            logger.warning("No supported preprocessing tasks remain after task arrangement.")

        return tasks

    def _resolve_input_spectrum_type(self) -> SpectrumType:
        """Resolve profile/centroid metadata before arranging tasks."""
        return resolve_spectrum_type(self.data_manager)

    @staticmethod
    def _resolve_output_spectrum_type(
        *,
        input_spectrum_type: SpectrumType,
        tasks: list[PreprocessTask],
    ) -> SpectrumType:
        """Infer output spectrum type after applying the provided tasks."""
        return "centroid" if any(task.name == "peak_pick" for task in tasks) else input_spectrum_type

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
            for mz_data, intensity, lengths, coords in data_manager.flat_generator(
                batch_size=self.batch_size, include_mz=True
            ):
                if stop_event.is_set():
                    break
                queue_ab.put(
                    FlatChunk(
                        mz_data=np.asarray(mz_data),
                        intensity=np.asarray(intensity),
                        lengths=np.asarray(lengths, dtype=np.uint32),
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
        sink: StageCSink,
        written_batches: int,
        total_batches: int,
    ) -> int:
        sink.write_chunk(chunk)
        written_batches += 1

        if total_batches <= 10 or written_batches % max(1, total_batches // 10) == 0 or written_batches == total_batches:
            progress = (written_batches / total_batches) * 100 if total_batches > 0 else 100.0
            logger.info(f"async_pipeline progress: {written_batches}/{total_batches} batches ({progress:.1f}%)")

        return written_batches

    def _stage_c_writer(
        self,
        *,
        queue_bc: Queue[Optional[FlatChunk]],
        stop_event: Event,
        sink: StageCSink,
        total_batches: int,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
        """Stage C: consumer; persist processed chunks."""
        written_batches = 0

        try:
            while True:
                chunk = queue_bc.get()
                if chunk is None:
                    break

                if stop_event.is_set():
                    continue

                written_batches = self._flush_chunk(
                    chunk=chunk,
                    sink=sink,
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
        output_spectrum_type: SpectrumType,
        final_output_format: str | None = None,
    ) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5:
        total_batches = (len(data_manager.ms) + self.batch_size - 1) // self.batch_size

        sink: StageCSink
        if final_output_format in self._ION_IMAGE:
            sink = StageCSink(self._make_ion_image_dm(data_manager, final_output_format), kind="ion_image")
        else:
            processed_ms = MassSpectrumSet()
            processed_data_manager = MSDataManagerImzML(processed_ms, temp_dir=self.temp_dir)
            processed_data_manager.copy_meta(data_manager)
            set_spectrum_type_metadata(processed_data_manager.ms.meta, output_spectrum_type)
            sink = StageCSink(processed_data_manager, kind="imzml")

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
                "sink": sink,
                "total_batches": total_batches,
                "error_holder": error_holder,
                "error_lock": error_lock,
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
            sink.close()
            raise RuntimeError("Async preprocessing pipeline failed.") from error_holder[0]

        return sink.finalize()

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
                min_count=kwargs.get("min_count"),
                min_frequency=kwargs.get("min_frequency"),
                batch_size=self.batch_size,
            )
        else:
            tolerance = tolerance * 1e-6 if units == "ppm" else tolerance

        kwargs["reference"] = reference
        kwargs["tolerance"] = tolerance

        kwargs.pop("binfun", None)
        kwargs.pop("binratio", None)
        kwargs.pop("min_count", None)
        kwargs.pop("min_frequency", None)

        return PreprocessTask(
            name=task.name,
            scope="batch",
            apply_fn=task.apply_fn,
            requires_dataset_prepare=False,
            kwargs=kwargs,
        )

    def _run_dataset_task(
        self,
        *,
        data_manager: MSDataManager,
        task: PreprocessTask,
        output_spectrum_type: SpectrumType,
        final_output_format: str | None = None,
    ) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5:
        """Execute whole-dataset task."""
        if task.name == "peak_align" and task.apply_fn is FlatPreprocess.peak_align_flat:
            prepared_task = self._prepare_peak_align_task(
                data_manager=data_manager,
                task=task,
            )
            return self._run_flat_batch_task(
                data_manager=data_manager,
                tasks=[prepared_task],
                output_spectrum_type=output_spectrum_type,
                final_output_format=final_output_format,
            )

        if task.name == "peak_align" or task.name == "peak_pick":
            dataset_apply_fn = cast(Callable[..., MSDataManagerImzML], task.apply_fn)
            result = dataset_apply_fn(data_manager=data_manager, **task.kwargs)
            set_spectrum_type_metadata(result.ms.meta, output_spectrum_type)
            return result
        else:
            raise NotImplementedError(f"Dataset task {task.name} with backend Cardinal is not supported.")

    def _default_ion_image_output_path(self, suffix: str = ".iib") -> str:
        if self.output_path is not None:
            return self.output_path
        base_dir = self.temp_dir if self.temp_dir is not None else tempfile.gettempdir()
        return f"{base_dir}/ion_image_output{suffix}"

    def _make_ion_image_dm(
        self, data_manager: MSDataManager, fmt: str
    ) -> IonImageDataManagerBinary | IonImageDataManagerH5:
        cls, suffix = self._ION_IMAGE[fmt]
        return cls.from_data_manager(
            data_manager, self._default_ion_image_output_path(suffix),
            temp_dir=self.temp_dir,
            temp_workers=max(1, self.output_max_workers),
            merge_workers=max(1, self.output_max_workers),
        )

    def _export_ion_image(
        self, data_manager: MSDataManager
    ) -> IonImageDataManagerBinary | IonImageDataManagerH5:
        return self._write_through_generator(
            data_manager,
            self._make_ion_image_dm(data_manager, self.output_format)
        )

    def _write_through_generator(
        self,
        data_manager: MSDataManager,
        output_dm: IonImageDataManagerBinary | IonImageDataManagerH5,
    ) -> IonImageDataManagerBinary | IonImageDataManagerH5:
        """通过 flat_generator 将数据写入 output_dm。"""
        try:
            for mz_axis, intensity_flat, lengths, coordinates in data_manager.flat_generator(
                batch_size=self.batch_size,
                include_mz=True,
                max_threads=max(2, self.output_max_workers),
            ):
                output_dm.write_flat_batch(
                    intensity_flat=intensity_flat,
                    lengths=lengths,
                    coordinates=coordinates,
                    mz_axis=mz_axis,
                )
            output_dm.finalize()
            return output_dm
        except BaseException:
            output_dm.close()
            raise

    def _maybe_close(self, dm: MSDataManager) -> None:
        if dm is not self.data_manager:
            dm.close()

    def start(self) -> Any:
        """Execute registered tasks and return new data manager with processed data."""
        current_data_manager: MSDataManager = self.data_manager
        current_spectrum_type = self._resolve_input_spectrum_type()
        ordered_tasks = self._arrange_tasks(current_spectrum_type)
        task_cursor = 0
        _ion_out = self.output_format in self._ION_IMAGE

        while task_cursor < len(ordered_tasks):
            batch_tasks: list[PreprocessTask] = []

            while (
                task_cursor < len(ordered_tasks)
                and ordered_tasks[task_cursor].scope == "batch"
                and not ordered_tasks[task_cursor].requires_dataset_prepare
            ):
                batch_tasks.append(ordered_tasks[task_cursor])
                task_cursor += 1

            if batch_tasks:
                fmt = self.output_format if (_ion_out and task_cursor >= len(ordered_tasks)) else None
                output_spectrum_type = self._resolve_output_spectrum_type(
                    input_spectrum_type=current_spectrum_type,
                    tasks=batch_tasks,
                )
                next_data_manager = self._run_flat_batch_task(
                    data_manager=current_data_manager,
                    tasks=batch_tasks,
                    output_spectrum_type=output_spectrum_type,
                    final_output_format=fmt,
                )
                self._maybe_close(current_data_manager)
                if fmt is not None:
                    return next_data_manager
                current_data_manager = cast(MSDataManager, next_data_manager)
                current_spectrum_type = output_spectrum_type

            if task_cursor >= len(ordered_tasks):
                break

            prepared_batch_task = ordered_tasks[task_cursor]
            if prepared_batch_task.scope == "batch" and prepared_batch_task.requires_dataset_prepare:
                task_cursor += 1
                fmt = self.output_format if (_ion_out and task_cursor >= len(ordered_tasks)) else None
                output_spectrum_type = self._resolve_output_spectrum_type(
                    input_spectrum_type=current_spectrum_type,
                    tasks=[prepared_batch_task],
                )
                next_data_manager = self._run_flat_batch_task(
                    data_manager=current_data_manager,
                    tasks=[
                        self._prepare_peak_align_task(
                            data_manager=current_data_manager,
                            task=prepared_batch_task,
                        )
                    ],
                    output_spectrum_type=output_spectrum_type,
                    final_output_format=fmt,
                )
                self._maybe_close(current_data_manager)
                if fmt is not None:
                    return next_data_manager
                current_data_manager = cast(MSDataManager, next_data_manager)
                current_spectrum_type = output_spectrum_type
                continue

            dataset_task = ordered_tasks[task_cursor]
            task_cursor += 1
            fmt = self.output_format if (_ion_out and task_cursor >= len(ordered_tasks)) else None
            output_spectrum_type = self._resolve_output_spectrum_type(
                input_spectrum_type=current_spectrum_type,
                tasks=[dataset_task],
            )
            next_data_manager = self._run_dataset_task(
                data_manager=current_data_manager,
                task=dataset_task,
                output_spectrum_type=output_spectrum_type,
                final_output_format=fmt,
            )
            self._maybe_close(current_data_manager)
            if fmt is not None and not isinstance(next_data_manager, MSDataManagerImzML):
                return next_data_manager
            current_data_manager = cast(MSDataManager, next_data_manager)
            current_spectrum_type = output_spectrum_type

        if _ion_out:
            output_dm = self._export_ion_image(current_data_manager)
            self._maybe_close(current_data_manager)
            return output_dm

        return cast(MSDataManagerImzML, current_data_manager)
