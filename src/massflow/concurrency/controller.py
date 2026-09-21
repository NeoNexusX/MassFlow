"""Process-wide admission and worker budgets for MassFlow execution phases."""

from __future__ import annotations

import os
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from enum import Enum
from threading import RLock, get_ident
from typing import Generator, Mapping, TypedDict, overload

from massflow.concurrency.capacity import get_logical_cpu_count
from massflow.tools.logger import get_logger

logger = get_logger("massflow.concurrency.controller")
POOL_KEYS = ("compute", "io_read", "io_write", "openmp")
_PHASE_RATIOS: dict[str, dict[str, float]] = {
    "preprocess_pipeline": {"io_read": 0.2, "compute": 0.6, "io_write": 0.2},
    "dataset_task": {"compute": 0.8, "io_read": 0.2},
    "io_transfer": {"io_read": 0.2, "io_write": 0.8},
    "io_read": {"io_read": 1.0},
    "compute": {"compute": 1.0},
    "umap": {"compute": 1.0},
    "kmeans": {"openmp": 1.0},
}
_BUDGET_ENV_VAR = "MASSFLOW_BUDGET"


class _Unset(Enum):
    VALUE = "unset"


class Status(TypedDict):
    """Configuration and the calling context's effective phase allocation."""

    budget: int | None
    pool_overrides: dict[str, int]
    phase: str | None
    allocations: dict[str, int]
    effective_budget: int | None
    sequential: bool


@dataclass(frozen=True)
class _Phase:
    name: str
    allocations: dict[str, int]
    bounded: bool
    owner: int
    budget: int | None
    sequential: bool


@dataclass
class _Configuration:
    budget: int | None = None
    overrides: dict[str, int] = field(default_factory=dict)
    running: int = 0


_STATE_LOCK = RLock()
_RUN_LOCK = RLock()
_ACTIVE: ContextVar[_Phase | None] = ContextVar("massflow_phase", default=None)


def _read_initial_budget() -> int:
    """Use an environment override or fall back to detected CPU capacity."""
    raw = os.environ.get(_BUDGET_ENV_VAR)
    if not raw:
        return get_logical_cpu_count()
    try:
        value = int(raw)
        if value < 1:
            raise ValueError
    except ValueError:
        logger.warning("Ignoring invalid %s=%r; must be a positive integer.", _BUDGET_ENV_VAR, raw)
        return get_logical_cpu_count()
    return value


_CONFIG = _Configuration(budget=_read_initial_budget())


# ---------------------------------------------------------------------------
# Budget configuration and allocation
# ---------------------------------------------------------------------------


def _validate_workers(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("workers must be an integer")
    if value < 1:
        raise ValueError("workers must be >= 1")


def _compute_allocations(phase_name: str, budget: int) -> dict[str, int]:
    ratios = _PHASE_RATIOS[phase_name]
    if budget < len(ratios):
        # Sequential stages reuse the budget instead of reserving a worker
        # per stage. Single-worker writes finish before the next batch starts.
        return {pool: budget if pool == "compute" else 1 for pool in ratios}
    raw = {pool: budget * ratio for pool, ratio in ratios.items()}
    allocations = {pool: max(1, int(value)) for pool, value in raw.items()}
    while sum(allocations.values()) > budget:
        largest = max(allocations, key=allocations.__getitem__)
        allocations[largest] -= 1
    for pool in sorted(raw, key=lambda key: raw[key] - int(raw[key]), reverse=True):
        if sum(allocations.values()) == budget:
            break
        allocations[pool] += 1
    return allocations


def set_budget(
    budget: int | None | _Unset = _Unset.VALUE,
    *,
    pools: Mapping[str, int | None] | None = None,
) -> Status:
    """Configure idle execution. Omitted budget preserves it; None disables it.

    Pool and per-call overrides are caps on the phase allocation. Configuration
    cannot change during execution: existing pools and masks cannot be resized
    atomically. Apply changes between runs instead.
    """
    if budget is not None and budget is not _Unset.VALUE:
        _validate_workers(budget)
    updates = {} if pools is None else pools
    for pool, value in updates.items():
        if pool not in POOL_KEYS:
            raise ValueError(f"Unknown pool key: {pool!r}")
        if value is not None:
            _validate_workers(value)
    with _STATE_LOCK:
        if _CONFIG.running:
            raise RuntimeError("Cannot change concurrency configuration during an active phase.")
        if budget is not _Unset.VALUE:
            _CONFIG.budget = budget
            if budget is None:
                _CONFIG.overrides.clear()
        for pool, value in updates.items():
            if value is None:
                _CONFIG.overrides.pop(pool, None)
            else:
                _CONFIG.overrides[pool] = value
    return get_status()


def get_status() -> Status:
    """Return a copy of the configuration and this context's active phase."""
    active = _ACTIVE.get()
    with _STATE_LOCK:
        return {
            "budget": _CONFIG.budget,
            "pool_overrides": dict(_CONFIG.overrides),
            "phase": active.name if active else None,
            "allocations": dict(active.allocations) if active else {},
            "effective_budget": get_effective_budget(),
            "sequential": runs_sequentially(),
        }


def get_effective_budget() -> int | None:
    """Return this phase's worker budget after CPU-capacity resolution."""
    active = _ACTIVE.get()
    if active is not None:
        return active.budget
    with _STATE_LOCK:
        budget = _CONFIG.budget
    return None if budget is None else min(budget, get_logical_cpu_count())


def runs_sequentially() -> bool:
    """Whether this phase reuses workers across its stages instead of overlapping them."""
    active = _ACTIVE.get()
    return active.sequential if active is not None else False


@overload
def get_allocation(pool: str, default: int) -> int: ...


@overload
def get_allocation(pool: str, default: None = None) -> int | None: ...


def get_allocation(pool: str, default: int | None = None) -> int | None:
    """Resolve a pool cap; inactive pools in a phase run sequentially."""
    if pool not in POOL_KEYS:
        raise ValueError(f"Unknown pool key: {pool!r}")
    active = _ACTIVE.get()
    if active is not None:
        # UMAP reads before computing, using one IO worker; other incidental
        # library pools must not introduce another independently sized pool.
        return active.allocations.get(pool, 1 if active.bounded else default)
    with _STATE_LOCK:
        budget = _CONFIG.budget
        override = _CONFIG.overrides.get(pool)
    if budget is not None:
        limit = min(budget, get_logical_cpu_count())
        return min(override, limit) if override is not None else limit
    return override if override is not None else default


def resolve_workers(pool: str, requested: int | None, *, default: int) -> int:
    """Resolve an API parameter without allowing it to exceed a pool cap."""
    if requested is not None:
        _validate_workers(requested)
    cap = get_allocation(pool)
    if cap is None:
        return default if requested is None else requested
    return cap if requested is None else min(requested, cap)


# ---------------------------------------------------------------------------
# Phase admission and runtime configuration
# ---------------------------------------------------------------------------


@contextmanager
def phase(
    name: str,
    *,
    pools: Mapping[str, int | None] | None = None,
) -> Generator[dict[str, int], None, None]:
    """Serialize controlled top-level phases; preserve nested caller contexts.

    Internal worker threads must receive a copy_context() snapshot. They read
    their phase allocation without acquiring the admission lock. Nested phases
    on the owner thread replace allocations only while the parent is suspended.
    """
    if name not in _PHASE_RATIOS:
        raise ValueError(f"Unknown phase: {name!r}")
    overrides = {} if pools is None else pools
    for pool, value in overrides.items():
        if pool not in POOL_KEYS:
            raise ValueError(f"Unknown pool key: {pool!r}")
        if value is not None:
            _validate_workers(value)

    active = _ACTIVE.get()
    if active is not None and active.owner != get_ident():
        raise RuntimeError(
            "Worker threads must use their inherited phase; only its owner can enter a nested phase."
        )

    _RUN_LOCK.acquire()
    with _STATE_LOCK:
        budget = _CONFIG.budget
        configured = dict(_CONFIG.overrides)
        _CONFIG.running += 1
    controlled = (
        budget is not None
        or bool(configured)
        or any(value is not None for value in overrides.values())
    )
    if not controlled:
        _RUN_LOCK.release()
    try:
        if not controlled:
            yield {}
            return
        effective_budget = None if budget is None else min(budget, get_logical_cpu_count())
        allocations = _resolve_phase_allocations(name, effective_budget, configured, overrides)
        token = _ACTIVE.set(
            _Phase(
                name,
                allocations,
                bounded=budget is not None,
                owner=get_ident(),
                budget=effective_budget,
                sequential=effective_budget is not None
                and effective_budget < len(_PHASE_RATIOS[name]),
            )
        )
        try:
            # Imported lazily: the runtime bridge queries this controller.
            # pylint: disable-next=import-outside-toplevel
            from massflow.concurrency.runtime import runtime_limits

            with runtime_limits(
                allocations,
                bounded=budget is not None,
                sequential=effective_budget is not None
                and effective_budget < len(_PHASE_RATIOS[name]),
            ):
                yield dict(allocations)
        finally:
            _ACTIVE.reset(token)
    finally:
        with _STATE_LOCK:
            _CONFIG.running -= 1
        if controlled:
            _RUN_LOCK.release()


def _resolve_phase_allocations(
    name: str,
    budget: int | None,
    configured: dict[str, int],
    overrides: Mapping[str, int | None],
) -> dict[str, int]:
    requested = configured | {pool: value for pool, value in overrides.items() if value is not None}
    if budget is None:
        return requested
    allocations = _compute_allocations(name, budget)
    for pool, value in requested.items():
        if pool in allocations:
            allocations[pool] = min(allocations[pool], value)
    return allocations
