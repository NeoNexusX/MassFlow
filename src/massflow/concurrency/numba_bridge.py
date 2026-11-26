"""Numba thread-mask bridge for the concurrency controller.

Numba's thread pool has a fixed size; the only lever is the calling
thread's mask ( numba.set_num_threads ), so every worker thread that runs
numba parallel kernels must apply its mask itself via apply_compute_mask().
"""

from __future__ import annotations

from numba import config, get_num_threads, set_num_threads

from massflow.concurrency.capacity import (
    detect_performance_core_workers,
    get_logical_cpu_count,
)
from massflow.concurrency.controller import resolve_workers
from massflow.tools.logger import get_logger

logger = get_logger("massflow.concurrency.numba_bridge")


def resolve_compute_workers(target: int) -> int:
    """Clamp a compute worker count to the available CPU and Numba pool limit."""
    if not isinstance(target, int):
        raise TypeError("compute workers must be an integer")
    if target < 1:
        raise ValueError("compute workers must be >= 1")

    limit = min(get_logical_cpu_count(), config.NUMBA_NUM_THREADS)  # type: ignore[attr-defined]  # pylint: disable=no-member
    if target > limit:
        logger.warning(
            "Requested compute workers=%d exceeds the available limit=%d; using %d workers.",
            target,
            limit,
            limit,
        )
        return limit
    return target


def apply_mask(target: int) -> int:
    """Set the calling thread's numba mask to *target*; return the previous mask."""
    previous = get_num_threads()
    resolved = resolve_compute_workers(target)
    if resolved != previous:
        set_num_threads(resolved)
    return previous


def restore_mask(previous: int) -> None:
    """Restore the calling thread's numba mask to a previously saved value."""
    set_num_threads(previous)


def apply_compute_mask() -> int:
    """Apply the controller's compute allocation to the calling thread."""
    target = resolve_workers("compute", None, default=detect_performance_core_workers())
    target = resolve_compute_workers(target)
    if get_num_threads() != target:
        set_num_threads(target)
    return target
