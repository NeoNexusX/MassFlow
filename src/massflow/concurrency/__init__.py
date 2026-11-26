"""Top-level concurrency control for MassFlow.

Public API:
- set_budget / get_status: global worker budget and per-pool overrides
- phase: context manager activating a processing phase's pool allocations
- get_allocation: pool lookup for integration points ( default-aware )
- apply_compute_mask: apply the compute allocation to the calling thread's
  numba mask
- get_logical_cpu_count / detect_performance_core_workers: CPU capacity
"""

from massflow.concurrency.capacity import (
    detect_performance_core_workers,
    get_logical_cpu_count,
)
from massflow.concurrency.controller import (
    get_allocation,
    get_effective_budget,
    runs_sequentially,
    get_status,
    phase,
    set_budget,
    resolve_workers,
)
from massflow.concurrency.numba_bridge import (
    apply_compute_mask,
    apply_mask,
    resolve_compute_workers,
    restore_mask,
)

__all__ = [
    "apply_compute_mask",
    "apply_mask",
    "detect_performance_core_workers",
    "get_allocation",
    "get_effective_budget",
    "runs_sequentially",
    "get_logical_cpu_count",
    "get_status",
    "phase",
    "resolve_compute_workers",
    "resolve_workers",
    "restore_mask",
    "set_budget",
]
