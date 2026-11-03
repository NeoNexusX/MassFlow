"""Temporary pixel-major to ion-major transpose storage."""

from __future__ import annotations

import os
import shutil
import tempfile
import threading
import uuid
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from massflow.concurrency import resolve_workers
from massflow.concurrency.runtime import get_io_executor

_BINARY_OPEN_FLAG = getattr(os, "O_BINARY", 0)
_OS_PREAD = getattr(os, "pread", None)
_OS_PWRITE = getattr(os, "pwrite", None)
_POSITIONED_IO_LOCK = threading.Lock()
_TEMP_WRITE_ION_BLOCK = 2048
_BATCH_FILE_NAME = "batches.bin"


def positioned_read(fd: int, nbytes: int, offset: int) -> bytes:
    """Read exactly from an explicit file offset."""
    if _OS_PREAD is not None:
        return _OS_PREAD(fd, nbytes, offset)
    with _POSITIONED_IO_LOCK:
        os.lseek(fd, offset, os.SEEK_SET)
        return os.read(fd, nbytes)


def positioned_write(fd: int, data: Any, offset: int) -> None:
    """Write the complete buffer at an explicit file offset."""
    view = memoryview(data).cast("B")
    if _OS_PWRITE is not None:
        while view:
            written = _OS_PWRITE(fd, view, offset)
            if written == 0:
                raise OSError("Positioned write returned zero bytes.")
            view = view[written:]
            offset += written
        return
    with _POSITIONED_IO_LOCK:
        os.lseek(fd, offset, os.SEEK_SET)
        while view:
            written = os.write(fd, view)
            if written == 0:
                raise OSError("Write returned zero bytes.")
            view = view[written:]


@dataclass(frozen=True, slots=True)
class BatchInfo:
    """One transposed temporary batch."""

    path: Path
    n_pixels: int
    global_pixel_start: int


class IonMajorTempStore:
    """Own temporary transposed batches used by ion-major writers.

    All batches live in one shared file: batch ``k`` occupies the contiguous
    ``[n_ions, n_pixels]`` segment starting at
    ``global_pixel_start * n_ions * itemsize``. Segments never overlap, so
    concurrent batch writes need no coordination, and the store holds exactly
    one open descriptor regardless of batch count (BUG-05).
    """

    def __init__(
        self,
        *,
        temp_dir: str | os.PathLike[str] | None,
        dtype: np.dtype | type,
        workers: int,
    ) -> None:
        base = Path(temp_dir) if temp_dir is not None else Path(tempfile.gettempdir())
        base.mkdir(parents=True, exist_ok=True)
        self._root = base / f"massflow_ion_major_{uuid.uuid4().hex}"
        self._root.mkdir(parents=True, exist_ok=False)
        self._batch_path = self._root / _BATCH_FILE_NAME
        self.dtype = np.dtype(dtype)
        self.workers = resolve_workers("io_write", workers, default=4)
        self._executor: ThreadPoolExecutor | None = None
        self._pending: set[Future[tuple[int, BatchInfo]]] = set()
        self._batches: list[tuple[int, BatchInfo]] = []
        self._batch_index = 0
        self._pixel_count = 0
        self._n_ions: int | None = None
        self._fd: int | None = None
        self._finished = False
        self._closed = False

    def _ensure_fd(self) -> int:
        if self._fd is None:
            self._fd = os.open(
                self._batch_path,
                os.O_RDWR | os.O_CREAT | _BINARY_OPEN_FLAG,
                0o600,
            )
        return self._fd

    def append_batch(
        self,
        intensity_flat: np.ndarray,
        *,
        n_pixels: int,
        n_ions: int,
    ) -> None:
        """Transpose and queue one owned pixel-major batch."""
        if self._closed or self._finished:
            raise ValueError("Cannot append after finish or close.")
        n_pixels = int(n_pixels)
        n_ions = int(n_ions)
        if self._n_ions is None:
            self._n_ions = n_ions
        elif n_ions != self._n_ions:
            raise ValueError(
                f"n_ions must stay constant across batches: "
                f"{self._n_ions} != {n_ions}."
            )
        values = np.array(
            intensity_flat,
            dtype=self.dtype,
            order="C",
            copy=True,
        )
        if values.ndim != 1 or values.size != n_pixels * n_ions:
            raise ValueError("intensity_flat size must equal n_pixels * n_ions.")

        batch_index = self._batch_index
        self._batch_index += 1
        info = BatchInfo(self._batch_path, n_pixels, self._pixel_count)
        self._pixel_count += n_pixels
        fd = self._ensure_fd()
        segment_offset = info.global_pixel_start * n_ions * self.dtype.itemsize
        executor = get_io_executor()
        if executor is None:
            if self._executor is None:
                self._executor = ThreadPoolExecutor(max_workers=self.workers)
            executor = self._executor
        future = executor.submit(
            self._write_batch,
            batch_index,
            info,
            fd,
            segment_offset,
            values,
            n_ions,
        )
        self._pending.add(future)
        workers = resolve_workers("io_write", self.workers, default=4)
        while len(self._pending) >= workers:
            self._collect(return_when=FIRST_COMPLETED)

    def finish(self) -> list[BatchInfo]:
        """Wait for queued writes and return batches in input order."""
        if self._closed:
            raise ValueError("Temporary store is closed.")
        while self._pending:
            self._collect(return_when=FIRST_COMPLETED)
        self._batches.sort(key=lambda item: item[0])
        self._finished = True
        return [info for _, info in self._batches]

    def read_ion_block(
        self,
        batch: BatchInfo,
        ion_start: int,
        ion_end: int,
    ) -> np.ndarray:
        """Read ``[ion_start:ion_end, :]`` from one temporary batch."""
        if not self._finished:
            raise ValueError("finish() must be called before reading.")
        ion_start = int(ion_start)
        ion_end = int(ion_end)
        n_ions = ion_end - ion_start
        if ion_start < 0 or n_ions < 0:
            raise ValueError("Invalid ion block.")
        if self._fd is None or self._n_ions is None:
            raise ValueError("No batches were appended.")
        if ion_end > self._n_ions:
            raise ValueError(f"ion_end must not exceed n_ions ({self._n_ions}).")
        itemsize = self.dtype.itemsize
        nbytes = n_ions * batch.n_pixels * itemsize
        offset = (
            batch.global_pixel_start * self._n_ions * itemsize
            + ion_start * batch.n_pixels * itemsize
        )
        raw = positioned_read(self._fd, nbytes, offset)
        if len(raw) != nbytes:
            raise OSError(
                f"Short read from {self._batch_path}: "
                f"expected {nbytes}, got {len(raw)}"
            )
        return np.frombuffer(raw, dtype=self.dtype).reshape(n_ions, batch.n_pixels)

    def close(self) -> None:
        """Release workers, descriptors and temporary files."""
        if self._closed:
            return
        self._closed = True
        wait(self._pending)
        if self._executor is not None:
            self._executor.shutdown(wait=True)
        if self._fd is not None:
            os.close(self._fd)
            self._fd = None
        shutil.rmtree(self._root, ignore_errors=True)

    @staticmethod
    def _write_batch(
        batch_index: int,
        info: BatchInfo,
        fd: int,
        segment_offset: int,
        values: np.ndarray,
        n_ions: int,
    ) -> tuple[int, BatchInfo]:
        itemsize = values.dtype.itemsize
        matrix = values.reshape(info.n_pixels, n_ions)
        for ion_start in range(0, n_ions, _TEMP_WRITE_ION_BLOCK):
            ion_end = min(n_ions, ion_start + _TEMP_WRITE_ION_BLOCK)
            block = np.ascontiguousarray(matrix[:, ion_start:ion_end].T)
            positioned_write(
                fd,
                block,
                segment_offset + ion_start * info.n_pixels * itemsize,
            )
        return batch_index, info

    def _collect(self, *, return_when: str) -> None:
        done, pending = wait(self._pending, return_when=return_when)
        self._pending = set(pending)
        for future in done:
            self._batches.append(future.result())

    def __enter__(self) -> "IonMajorTempStore":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()
