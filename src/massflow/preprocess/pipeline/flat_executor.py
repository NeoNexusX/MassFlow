from __future__ import annotations

from queue import Queue
from threading import Event, Lock, Thread
from typing import Callable, Optional, cast
from dataclasses import dataclass

import numpy as np

from massflow.data_manager import MSDataManager, MSDataManagerImzML, IonImageDataManagerBinary,IonImageDataManagerH5
from massflow.preprocess.flat_pre_fun import FlatBatchResult
from massflow.preprocess.pipeline.types import (
    ION_IMAGE_OUTPUT_FORMATS,
    OUTPUT_FORMAT_ION_IMAGE_BINARY,
    InputMetaData,
    OutputMetaData,
    PreprocessTask,
)
from massflow.preprocess.numba.numba_runtime import apply_numba_runtime
from massflow.tools.infer_spectrum_type import set_spectrum_type_metadata
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.pipeline.flat_executor")


@dataclass(slots=True)
class FlatChunk:
    """Transport object between pipeline stages."""

    mz_data: Optional[np.ndarray]
    intensity: np.ndarray
    lengths: np.ndarray
    coordinates: np.ndarray

@dataclass(slots=True)
class FlatOutputSink:
    """Output sink for flat segments; ion-image writers are created lazily."""

    source_data_manager: MSDataManager
    input_meta: InputMetaData
    output_meta: OutputMetaData
    data_manager: MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5 | None = None

    @classmethod
    def create(
        cls,
        *,
        source_data_manager: MSDataManager,
        input_meta: InputMetaData,
        output_meta: OutputMetaData,
    ) -> "FlatOutputSink":
        sink = cls(
            source_data_manager=source_data_manager,
            input_meta=input_meta,
            output_meta=output_meta,
        )
        if output_meta.output_format not in ION_IMAGE_OUTPUT_FORMATS:
            processed_data = MSDataManagerImzML(temp_dir=input_meta.temp_dir)
            processed_data.copy_meta(source_data_manager)
            set_spectrum_type_metadata(processed_data.ms.meta, output_meta.spectrum_type)
            sink.data_manager = processed_data
        return sink

    def _make_ion_image_data_manager(self, chunk: FlatChunk) -> IonImageDataManagerBinary | IonImageDataManagerH5:
        if chunk.mz_data is None:
            raise ValueError("ion image output requires mz_data from processed flat chunks.")
        if self.output_meta.file_path is None:
            raise ValueError("ion image output requires a file_path.")

        ion_image_cls = (
            IonImageDataManagerBinary
            if self.output_meta.output_format == OUTPUT_FORMAT_ION_IMAGE_BINARY
            else IonImageDataManagerH5
        )
        return ion_image_cls.from_flat_output_contract(
            self.source_data_manager,
            self.output_meta.file_path,
            shared_mz_axis=chunk.mz_data,
            output_spectrum_type=self.output_meta.spectrum_type,
            temp_dir=self.input_meta.temp_dir,
            temp_workers=self.input_meta.output_max_workers,
            merge_workers=self.input_meta.output_max_workers,
        )

    def write_chunk(self, chunk: FlatChunk) -> None:
        if self.data_manager is None:
            self.data_manager = self._make_ion_image_data_manager(chunk)
        self.data_manager.write_flat_batch(
            intensity_flat=chunk.intensity,
            lengths=chunk.lengths,
            coordinates=chunk.coordinates,
            mz_axis=chunk.mz_data, # type: ignore
        )

    def finalize(self) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5:
        if self.data_manager is None:
            raise ValueError("No batches were written.")
        if isinstance(self.data_manager, MSDataManagerImzML):
            self.data_manager.close_writer()
            self.data_manager.load_head_data()
        else:
            self.data_manager.finalize()
        return self.data_manager

    def close(self) -> None:
        if isinstance(self.data_manager, MSDataManagerImzML):
            self.data_manager.close_writer()
        if self.data_manager is not None:
            self.data_manager.close()

def _set_error(
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
    input_meta: InputMetaData,
    data_manager: MSDataManager,
    queue_ab: Queue[Optional[FlatChunk]],
    stop_event: Event,
    error_holder: list[BaseException],
    error_lock: Lock,
) -> None:
    """Stage A: producer; reads flat and pushes to queue_ab."""
    try:
        for mz_data, intensity, lengths, coords in data_manager.flat_generator(
            batch_size=input_meta.batch_size, include_mz=True # include_mz must be True
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
        _set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc)
    finally:
        queue_ab.put(None)


def _stage_b_processor(
    input_meta: InputMetaData,
    tasks: list[PreprocessTask],
    queue_ab: Queue[Optional[FlatChunk]],
    queue_bc: Queue[Optional[FlatChunk]],
    stop_event: Event,
    error_holder: list[BaseException],
    error_lock: Lock,
) -> None:
    """Stage B: middleware; runs ordered preprocessing tasks for each batch."""
    try:
        apply_numba_runtime(override_workers=input_meta.numba_max_threads)
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
        _set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc)
        queue_bc.put(None)

def _flush_chunk(
    data_sink: FlatOutputSink,
    chunk: FlatChunk,
    written_batches: int,
    total_batches: int,
) -> int:
    data_sink.write_chunk(chunk)
    written_batches += 1

    if total_batches <= 10 or written_batches % max(1, total_batches // 10) == 0 or written_batches == total_batches:
        progress = (written_batches / total_batches) * 100 if total_batches > 0 else 100.0
        logger.info(f"async_pipeline progress: {written_batches}/{total_batches} batches ({progress:.1f}%)")

    return written_batches

def _stage_c_writer(
    data_sink: FlatOutputSink,
    queue_bc: Queue[Optional[FlatChunk]],
    stop_event: Event,
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

            written_batches = _flush_chunk(
                data_sink=data_sink,
                chunk=chunk,
                written_batches=written_batches,
                total_batches=total_batches,
            )

    except BaseException as exc:  # pylint: disable=broad-exception-caught
        _set_error(
            error_holder=error_holder,
            error_lock=error_lock,
            stop_event=stop_event,
            exc=exc,
        )

def run_flat_task(
    data_manager: MSDataManager,
    tasks: list[PreprocessTask],
    input_meta: InputMetaData,
    output_meta: OutputMetaData,
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5:
    total_batches = (len(data_manager.ms) + input_meta.batch_size - 1) // input_meta.batch_size

    data_sink = FlatOutputSink.create(
        source_data_manager=data_manager,
        input_meta=input_meta,
        output_meta=output_meta,
    )

    queue_ab: Queue[Optional[FlatChunk]] = Queue(maxsize=input_meta.queue_ab_size)
    queue_bc: Queue[Optional[FlatChunk]] = Queue(maxsize=input_meta.queue_bc_size)

    stop_event = Event()
    error_holder: list[BaseException] = []
    error_lock = Lock()

    t_a = Thread(
        target=_stage_a_reader,
        kwargs={
            "input_meta": input_meta,
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
        target=_stage_b_processor,
        kwargs={
            "tasks": tasks,
            "queue_ab": queue_ab,
            "queue_bc": queue_bc,
            "stop_event": stop_event,
            "input_meta": input_meta,
            "error_holder": error_holder,
            "error_lock": error_lock,
        },
        name="massflow-stage-b-processor",
        daemon=True,
    )
    t_c = Thread(
        target=_stage_c_writer,
        kwargs={
            "data_sink": data_sink,
            "queue_bc": queue_bc,
            "stop_event": stop_event,
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
        data_sink.close()
        raise RuntimeError("Async preprocessing pipeline failed.") from error_holder[0]

    return data_sink.finalize()
