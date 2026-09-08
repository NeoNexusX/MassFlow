from __future__ import annotations

from queue import Empty, Full, Queue
from threading import Event, Lock, Thread
from time import monotonic
from typing import Callable, Optional, cast
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from massflow.data_manager import (
    MSDataManager,
    MSDataManagerImzML,
    MSDataManagerZarr,
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
)
from massflow.preprocess.flat_pre_fun import FlatBatchResult
from massflow.preprocess.pipeline.types import (
    OUTPUT_FORMAT_IMZML,
    OUTPUT_FORMAT_ION_IMAGE_BINARY,
    OUTPUT_FORMAT_ION_IMAGE_H5,
    OUTPUT_FORMAT_ION_IMAGE_ZARR,
    OUTPUT_FORMAT_SPECTRA_ZARR,
    InputMetaData,
    OutputMetaData,
    PreprocessTask,
)
from massflow.preprocess.numba.numba_runtime import apply_numba_runtime
from massflow.tools.infer_spectrum_type import set_spectrum_type_metadata
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.pipeline.flat_executor")

_QUEUE_WAIT_SECONDS = 0.1
_THREAD_SHUTDOWN_SECONDS = 1.0


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
    data_manager: MSDataManagerImzML | MSDataManagerZarr | IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr | None = None

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
        if output_meta.output_format in {
            OUTPUT_FORMAT_IMZML,
            OUTPUT_FORMAT_SPECTRA_ZARR,
        }:
            if output_meta.file_path is not None:
                Path(output_meta.file_path).parent.mkdir(parents=True, exist_ok=True)
            if output_meta.output_format == OUTPUT_FORMAT_SPECTRA_ZARR:
                processed_data: MSDataManagerImzML | MSDataManagerZarr = (
                    MSDataManagerZarr(
                        filepath=output_meta.file_path,
                        mode="w",
                        temp_dir=input_meta.temp_dir,
                        chunk_size=input_meta.zarr_chunk_size,
                        compute_tic=input_meta.compute_tic,
                        compute_rms=input_meta.compute_rms,
                    )
                )
            else:
                processed_data = MSDataManagerImzML(
                    filepath=output_meta.file_path,
                    temp_dir=input_meta.temp_dir,
                )
            processed_data.copy_meta(source_data_manager)
            set_spectrum_type_metadata(processed_data.ms.meta, output_meta.spectrum_type)
            sink.data_manager = processed_data
        return sink

    def _make_ion_image_data_manager(self, chunk: FlatChunk) -> IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr:
        if chunk.mz_data is None:
            raise ValueError("ion image output requires mz_data from processed flat chunks.")
        if self.output_meta.file_path is None:
            raise ValueError("ion image output requires a file_path.")

        if self.output_meta.output_format == OUTPUT_FORMAT_ION_IMAGE_BINARY:
            ion_image_cls = IonImageDataManagerBinary
        elif self.output_meta.output_format == OUTPUT_FORMAT_ION_IMAGE_H5:
            ion_image_cls = IonImageDataManagerH5
        elif self.output_meta.output_format == OUTPUT_FORMAT_ION_IMAGE_ZARR:
            ion_image_cls = IonImageDataManagerZarr
        else:
            raise ValueError(f"Unsupported ion image output format: {self.output_meta.output_format}")

        zarr_options = (
            {
                "chunk_size": self.input_meta.zarr_chunk_size,
                "compute_rms": self.input_meta.compute_rms,
            }
            if self.output_meta.output_format == OUTPUT_FORMAT_ION_IMAGE_ZARR
            else {}
        )
        return ion_image_cls.from_flat_output_contract(
            self.source_data_manager,
            self.output_meta.file_path,
            shared_mz_axis=chunk.mz_data,
            output_spectrum_type=self.output_meta.spectrum_type,
            temp_dir=self.input_meta.temp_dir,
            temp_workers=self.input_meta.output_max_workers,
            merge_workers=self.input_meta.output_max_workers,
            ion_block=self.input_meta.ion_block,
            compute_mean_spectrum=self.input_meta.compute_mean_spectrum,
            **zarr_options, # type: ignore
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

    def finalize(self) -> MSDataManagerImzML | MSDataManagerZarr | IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr:
        if self.data_manager is None:
            raise ValueError("No batches were written.")
        if isinstance(self.data_manager, MSDataManagerImzML):
            self.data_manager.close_writer()
            self.data_manager.load_head_data()
        elif isinstance(self.data_manager, MSDataManagerZarr):
            self.data_manager.finalize()
            self.data_manager.close()
            result = MSDataManagerZarr(
                filepath=self.output_meta.file_path,
                mode="r",
            )
            result.load_head_data()
            self.data_manager = result
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

def _mark_progress(last_progress: list[float], progress_lock: Lock) -> None:
    with progress_lock:
        last_progress[0] = monotonic()

def _put_interruptibly(
    queue: Queue[Optional[FlatChunk]],
    item: Optional[FlatChunk],
    stop_event: Event,
    last_progress: list[float],
    progress_lock: Lock,
) -> bool:
    while not stop_event.is_set():
        try:
            queue.put(item, timeout=_QUEUE_WAIT_SECONDS)
        except Full:
            continue
        _mark_progress(last_progress, progress_lock)
        return True
    return False

def _stage_a_reader(
    input_meta: InputMetaData,
    data_manager: MSDataManager,
    queue_ab: Queue[Optional[FlatChunk]],
    stop_event: Event,
    error_holder: list[BaseException],
    error_lock: Lock,
    last_progress: list[float],
    progress_lock: Lock,
) -> None:
    """Stage A: producer; reads flat and pushes to queue_ab."""
    try:
        for mz_data, intensity, lengths, coords in data_manager.flat_generator(
            batch_size=input_meta.batch_size, include_mz=True # include_mz must be True
        ):
            if stop_event.is_set():
                break
            if not _put_interruptibly(
                queue_ab,
                FlatChunk(
                    mz_data=np.asarray(mz_data),
                    intensity=np.asarray(intensity),
                    lengths=np.asarray(lengths, dtype=np.uint32),
                    coordinates=np.asarray(coords),
                ),
                stop_event,
                last_progress,
                progress_lock,
            ):
                break
    except BaseException as exc:  # pylint: disable=broad-exception-caught
        _set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc)
    finally:
        _put_interruptibly(
            queue_ab,
            None,
            stop_event,
            last_progress,
            progress_lock,
        )


def _stage_b_processor(
    input_meta: InputMetaData,
    tasks: list[PreprocessTask],
    queue_ab: Queue[Optional[FlatChunk]],
    queue_bc: Queue[Optional[FlatChunk]],
    stop_event: Event,
    error_holder: list[BaseException],
    error_lock: Lock,
    last_progress: list[float],
    progress_lock: Lock,
) -> None:
    """Stage B: middleware; runs ordered preprocessing tasks for each batch."""
    try:
        apply_numba_runtime(override_workers=input_meta.numba_max_threads)
        while True:
            if stop_event.is_set():
                break
            try:
                chunk = queue_ab.get(timeout=_QUEUE_WAIT_SECONDS)
            except Empty:
                continue
            _mark_progress(last_progress, progress_lock)
            if chunk is None:
                break

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

            if not _put_interruptibly(
                queue_bc,
                current_chunk,
                stop_event,
                last_progress,
                progress_lock,
            ):
                break

    except BaseException as exc:  # pylint: disable=broad-exception-caught
        _set_error(error_holder=error_holder,error_lock=error_lock,stop_event=stop_event,exc=exc)
    finally:
        _put_interruptibly(
            queue_bc,
            None,
            stop_event,
            last_progress,
            progress_lock,
        )

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
    last_progress: list[float],
    progress_lock: Lock,
) -> None:
    """Stage C: consumer; persist processed chunks."""
    written_batches = 0

    try:
        while True:
            if stop_event.is_set():
                break
            try:
                chunk = queue_bc.get(timeout=_QUEUE_WAIT_SECONDS)
            except Empty:
                continue
            _mark_progress(last_progress, progress_lock)
            if chunk is None:
                break

            written_batches = _flush_chunk(
                data_sink=data_sink,
                chunk=chunk,
                written_batches=written_batches,
                total_batches=total_batches,
            )
            _mark_progress(last_progress, progress_lock)

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
) -> MSDataManagerImzML | MSDataManagerZarr | IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr:
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
    last_progress = [monotonic()]
    progress_lock = Lock()

    t_a = Thread(
        target=_stage_a_reader,
        kwargs={
            "input_meta": input_meta,
            "data_manager": data_manager,
            "queue_ab": queue_ab,
            "stop_event": stop_event,
            "error_holder": error_holder,
            "error_lock": error_lock,
            "last_progress": last_progress,
            "progress_lock": progress_lock,
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
            "last_progress": last_progress,
            "progress_lock": progress_lock,
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
            "last_progress": last_progress,
            "progress_lock": progress_lock,
        },
        name="massflow-stage-c-writer",
        daemon=True,
    )

    t_a.start()
    t_b.start()
    t_c.start()

    threads = (t_a, t_b, t_c)
    while any(thread.is_alive() for thread in threads) and not error_holder:
        stall_timeout = input_meta.pipeline_stall_timeout
        if stall_timeout is not None:
            with progress_lock:
                stalled_for = monotonic() - last_progress[0]
            if stalled_for >= stall_timeout:
                stalled_threads = ", ".join(
                    thread.name for thread in threads if thread.is_alive()
                )
                _set_error(
                    error_holder=error_holder,
                    error_lock=error_lock,
                    stop_event=stop_event,
                    exc=TimeoutError(
                        "Async preprocessing pipeline made no progress for "
                        f"{stall_timeout:g} seconds; active stages: {stalled_threads}."
                    ),
                )
                break
        for thread in threads:
            thread.join(timeout=_QUEUE_WAIT_SECONDS)

    if error_holder:
        stop_event.set()
        for thread in threads:
            thread.join(timeout=_THREAD_SHUTDOWN_SECONDS)
        if not t_c.is_alive():
            data_sink.close()
        raise RuntimeError("Async preprocessing pipeline failed.") from error_holder[0]

    return data_sink.finalize()
