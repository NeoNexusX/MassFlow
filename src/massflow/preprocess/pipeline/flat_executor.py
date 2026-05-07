from __future__ import annotations

from queue import Queue
from threading import Event, Lock, Thread
from typing import Callable, Optional, cast

import numpy as np

from massflow.data_manager import MSDataManager
from massflow.preprocess.flat_pre_fun import FlatBatchResult
from massflow.preprocess.numba.numba_runtime import apply_numba_runtime
from massflow.tools.logger import get_logger

from massflow.preprocess.pipeline.models import FlatChunk, PipelineConfig, PipelineResult, PreprocessTask
from massflow.preprocess.pipeline.sinks import OutputSink

logger = get_logger("massflow.preprocess.pipeline.flat_executor")


class FlatPipelineExecutor:
    """Run the three-stage flat batch pipeline against a supplied output sink."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    @staticmethod
    def _set_error(
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
        try:
            for mz_data, intensity, lengths, coords in data_manager.flat_generator(
                batch_size=self.config.batch_size,
                include_mz=True,
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
        except BaseException as exc:
            self._set_error(
                error_holder=error_holder,
                error_lock=error_lock,
                stop_event=stop_event,
                exc=exc,
            )
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
        try:
            apply_numba_runtime(override_workers=self.config.numba_max_threads)
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

        except BaseException as exc:
            self._set_error(
                error_holder=error_holder,
                error_lock=error_lock,
                stop_event=stop_event,
                exc=exc,
            )
            queue_bc.put(None)

    @staticmethod
    def _flush_chunk(
        *,
        chunk: FlatChunk,
        sink: OutputSink,
        written_batches: int,
        total_batches: int,
    ) -> int:
        sink.write_chunk(chunk)
        written_batches += 1

        if total_batches <= 10 or written_batches % max(1, total_batches // 10) == 0 or written_batches == total_batches:
            progress = (written_batches / total_batches) * 100 if total_batches > 0 else 100.0
            logger.info(f"flat_pipeline progress: {written_batches}/{total_batches} batches ({progress:.1f}%)")

        return written_batches

    def _stage_c_writer(
        self,
        *,
        queue_bc: Queue[Optional[FlatChunk]],
        stop_event: Event,
        sink: OutputSink,
        total_batches: int,
        error_holder: list[BaseException],
        error_lock: Lock,
    ) -> None:
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

        except BaseException as exc:
            self._set_error(
                error_holder=error_holder,
                error_lock=error_lock,
                stop_event=stop_event,
                exc=exc,
            )

    def run(
        self,
        *,
        data_manager: MSDataManager,
        tasks: list[PreprocessTask],
        sink: OutputSink,
    ) -> PipelineResult:
        total_batches = (len(data_manager.ms) + self.config.batch_size - 1) // self.config.batch_size
        queue_ab: Queue[Optional[FlatChunk]] = Queue(maxsize=self.config.queue_ab_size)
        queue_bc: Queue[Optional[FlatChunk]] = Queue(maxsize=self.config.queue_bc_size)

        stop_event = Event()
        error_holder: list[BaseException] = []
        error_lock = Lock()

        threads = [
            Thread(
                target=self._stage_a_reader,
                kwargs={
                    "data_manager": data_manager,
                    "queue_ab": queue_ab,
                    "stop_event": stop_event,
                    "error_holder": error_holder,
                    "error_lock": error_lock,
                },
                name="massflow-flat-reader",
                daemon=True,
            ),
            Thread(
                target=self._stage_b_processor,
                kwargs={
                    "tasks": tasks,
                    "queue_ab": queue_ab,
                    "queue_bc": queue_bc,
                    "stop_event": stop_event,
                    "error_holder": error_holder,
                    "error_lock": error_lock,
                },
                name="massflow-flat-processor",
                daemon=True,
            ),
            Thread(
                target=self._stage_c_writer,
                kwargs={
                    "queue_bc": queue_bc,
                    "stop_event": stop_event,
                    "sink": sink,
                    "total_batches": total_batches,
                    "error_holder": error_holder,
                    "error_lock": error_lock,
                },
                name="massflow-flat-writer",
                daemon=True,
            ),
        ]

        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        if error_holder:
            sink.close()
            raise RuntimeError("Flat preprocessing pipeline failed.") from error_holder[0]

        return sink.finalize()
