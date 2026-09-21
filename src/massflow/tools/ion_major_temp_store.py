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
_POSITIONED_IO_LOCK = threading.Lock()
_TEMP_WRITE_ION_BLOCK = 2048


def positioned_read(fd: int, nbytes: int, offset: int) -> bytes:
    """Read exactly from an explicit file offset."""
    if _OS_PREAD is not None:
        return _OS_PREAD(fd, nbytes, offset)
    with _POSITIONED_IO_LOCK:
        os.lseek(fd, offset, os.SEEK_SET)
        return os.read(fd, nbytes)


@dataclass(frozen=True, slots=True)
class BatchInfo:
    """One transposed temporary batch."""

    path: Path
    n_pixels: int
    global_pixel_start: int


class IonMajorTempStore:
    """Own temporary transposed batches used by ion-major writers."""

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
        self.dtype = np.dtype(dtype)
        self.workers = resolve_workers("io_write", workers, default=4)
        self._executor: ThreadPoolExecutor | None = None
        self._pending: set[Future[tuple[int, BatchInfo]]] = set()
        self._batches: list[tuple[int, BatchInfo]] = []
        self._batch_index = 0
        self._pixel_count = 0
        self._read_fds: dict[Path, int] = {}
        self._fd_lock = threading.Lock()
        self._finished = False
        self._closed = False

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
        batch_path = self._root / f"batch_{batch_index:08d}.bin"
        info = BatchInfo(batch_path, n_pixels, self._pixel_count)
        self._pixel_count += n_pixels
        executor = get_io_executor()
        if executor is None:
            if self._executor is None:
                self._executor = ThreadPoolExecutor(max_workers=self.workers)
            executor = self._executor
        future = executor.submit(
            self._write_batch,
            batch_index,
            info,
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
        with self._fd_lock:
            fd = self._read_fds.get(batch.path)
            if fd is None:
                fd = os.open(batch.path, os.O_RDONLY | _BINARY_OPEN_FLAG)
                self._read_fds[batch.path] = fd
        nbytes = n_ions * batch.n_pixels * self.dtype.itemsize
        offset = ion_start * batch.n_pixels * self.dtype.itemsize
        raw = positioned_read(fd, nbytes, offset)
        if len(raw) != nbytes:
            raise OSError(
                f"Short read from {batch.path}: expected {nbytes}, got {len(raw)}"
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
        for fd in self._read_fds.values():
            os.close(fd)
        self._read_fds.clear()
        shutil.rmtree(self._root, ignore_errors=True)

    @staticmethod
    def _write_batch(
        batch_index: int,
        info: BatchInfo,
        values: np.ndarray,
        n_ions: int,
    ) -> tuple[int, BatchInfo]:
        matrix = values.reshape(info.n_pixels, n_ions)
        with info.path.open("wb") as output:
            for ion_start in range(0, n_ions, _TEMP_WRITE_ION_BLOCK):
                ion_end = min(n_ions, ion_start + _TEMP_WRITE_ION_BLOCK)
                block = np.ascontiguousarray(matrix[:, ion_start:ion_end].T)
                output.write(memoryview(block))
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
