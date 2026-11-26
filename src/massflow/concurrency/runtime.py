"""Scoped native-library and Zarr executor limits for admitted phases."""

from __future__ import annotations

import asyncio
from concurrent.futures import Executor, Future, ThreadPoolExecutor, wait
from contextlib import ExitStack, contextmanager
from contextvars import ContextVar
from threading import Lock
from typing import Generator

from threadpoolctl import threadpool_limits
from zarr import config as zarr_config
from zarr.core.sync import sync

from massflow.concurrency.numba_bridge import apply_mask, restore_mask

_IO_EXECUTOR: ContextVar[ThreadPoolExecutor | None] = ContextVar("massflow_io_executor", default=None)


def get_io_executor() -> ThreadPoolExecutor | None:
    """Return the calling phase's shared IO pool without transferring ownership."""
    return _IO_EXECUTOR.get()


class _BorrowedExecutor(Executor):
    """Track this borrower's leaf tasks so resources outlive pending IO."""

    def __init__(self, executor: ThreadPoolExecutor) -> None:
        self._executor = executor
        self._pending: set[Future] = set()
        self._lock = Lock()

    def submit(self, fn, /, *args, **kwargs):
        future = self._executor.submit(fn, *args, **kwargs)
        with self._lock:
            self._pending.add(future)
        future.add_done_callback(self._discard)
        return future

    def _discard(self, future: Future) -> None:
        with self._lock:
            self._pending.discard(future)

    def drain(self) -> None:
        with self._lock:
            pending = tuple(self._pending)
        wait(pending)


@contextmanager
def io_executor(*, max_workers: int, thread_name_prefix: str = "") -> Generator[Executor, None, None]:
    """Borrow the phase pool for leaf IO, waiting for this scope's tasks on exit.

    Coordinators that synchronously wait for Zarr or other futures must stay
    outside this pool. Uncontrolled callers retain their own executor.
    """
    shared = get_io_executor()
    if shared is None:
        with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix=thread_name_prefix) as executor:
            yield executor
    else:
        borrowed = _BorrowedExecutor(shared)
        try:
            yield borrowed
        finally:
            borrowed.drain()


async def _swap_zarr_executor(executor: ThreadPoolExecutor | None) -> ThreadPoolExecutor | None:
    loop = asyncio.get_running_loop()
    # asyncio has no getter/restorer for an initially absent default executor.
    # Keep the private compatibility access confined to this backend adapter.
    previous = loop._default_executor  # type: ignore[attr-defined]  # pylint: disable=protected-access
    if executor is None:
        loop._default_executor = None  # type: ignore[attr-defined]  # pylint: disable=protected-access
    else:
        loop.set_default_executor(executor)
    return previous


@contextmanager
def runtime_limits(
    allocations: dict[str, int],
    *,
    bounded: bool,
    sequential: bool = False,
) -> Generator[None, None, None]:
    """Restore masks, native-library limits and the existing Zarr executor."""
    with ExitStack() as stack:
        compute = allocations.get("compute", 1 if bounded else None)
        if compute is not None:
            previous_mask = apply_mask(compute)
            stack.callback(restore_mask, previous_mask)
        if bounded or "openmp" in allocations:
            stack.enter_context(
                threadpool_limits(limits={"blas": 1, "openmp": allocations.get("openmp", 1)})
            )
        if bounded or "io_read" in allocations or "io_write" in allocations:
            io_allocations = [allocations[pool] for pool in ("io_read", "io_write") if pool in allocations]
            io_workers = (max(io_allocations, default=1) if sequential else sum(io_allocations) or 1)
            executor = stack.enter_context(
                ThreadPoolExecutor(max_workers=io_workers, thread_name_prefix="massflow-io")
            )
            token = _IO_EXECUTOR.set(executor)
            stack.callback(_IO_EXECUTOR.reset, token)
            stack.enter_context(_zarr_limits(executor, io_workers))
        yield


@contextmanager
def _zarr_limits(executor: ThreadPoolExecutor, workers: int) -> Generator[None, None, None]:
    """Share the phase's IO allocation with readers, writers and prefetch.

    Zarr's threading.max_workers only affects initial executor creation, so
    replace the event loop executor even when Zarr has already performed IO.
    The admission lock remains held until all pipeline work has finished.
    """
    with zarr_config.set({"async.concurrency": workers, "threading.max_workers": None}):
        previous_executor = sync(_swap_zarr_executor(executor))
        try:
            yield
        finally:
            sync(_swap_zarr_executor(previous_executor))
