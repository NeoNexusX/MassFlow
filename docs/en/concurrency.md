# Concurrency Control

This page explains how to set a concurrency budget through `massflow.concurrency` and how MassFlow distributes it across data reading, preprocessing, writing, and UMAP/KMeans computation.

Authors: Leqiang Dong, Jingwei Dong

**Feedback:** If you find any problems while using MassFlow or reading the documentation, you can submit feedback through [GitHub Issues](https://github.com/NeoNexusX/MassFlow/issues), or contact us by email at 707244007@qq.com.

## Overview

MassFlow uses reading threads, Numba computation threads, and writing threads to process data. A concurrency budget coordinates this work so that each thread pool does not independently use the machine's full core count.

Concurrency control is enabled by default, with an initial budget based on the detected available CPU capacity. You can usually run preprocessing or segmentation directly. To leave computing resources for other programs, adjust the budget with `set_budget` before running a task.

- **Base budget**: the worker count used to scale each pool independently in the current Python process. The automatic startup configuration uses the detected CPU capacity as the base.
- **Processing phase**: a stage such as batch preprocessing, a whole-dataset task, or UMAP. Each phase uses its own allocation ratios.
- **Resource pool**: an allocation category for reading, computation, writing, or OpenMP. These categories do not necessarily have separate thread pools.

The budget limits concurrency for managed work. It is not the total process thread count shown by a system monitor: pipeline coordination threads, event loops, and idle library threads may still exist.

## Usage Example

Set the budget before preprocessing, then run the workflow as follows:

```python
import massflow.concurrency as mc
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor

mc.set_budget(8)
print(mc.get_status())

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()

normalized = (
    Preprocessor(dm, output_format="spectra_zarr", output_path="normalized.zarr")
    .normalization(method="tic_numba")
    .start()
)

normalized.close()
dm.close()
```

Here, `8` is the configured value. If the process can use only 4 logical cores, the effective budget is 4. Subsequent UMAP or KMeans calls use the same global configuration.

You can also set the initial budget before starting Python:

```bash
export MASSFLOW_BUDGET=8
```

`MASSFLOW_BUDGET` is read when the concurrency module is first imported and must be a positive integer. Invalid values are ignored, falling back to automatic detection. After import, use `set_budget` to change the budget.

## Function Description

### set_budget (Set the Base Budget and Pool Limits)

```python
mc.set_budget(8, pools={"compute": 4})
```

- Function: sets the base budget and optional pool limits for subsequent processing phases.
- Parameters:
  - `budget`: a positive integer or `None`. Omitting it preserves the current budget; passing `None` disables the base budget and clears existing pool settings.
  - `pools`: a dictionary with keys `compute`, `io_read`, `io_write`, or `openmp`. A positive integer sets the pool limit; `None` removes that pool's setting. Pools not listed keep their existing settings.
- Output: the same status dictionary returned by `get_status()`.

Common configuration options:

```python
mc.set_budget(8)                          # Set the base budget; keep existing pool settings
mc.set_budget(pools={"compute": 4})       # Change only the compute pool limit
mc.set_budget(pools={"compute": None})    # Remove the compute pool limit
mc.set_budget(None)                       # Disable the base budget and clear all pool settings
mc.set_budget(mc.get_logical_cpu_count()) # Enable using the current available CPU capacity
```

With a base budget, pool settings can only reduce a phase's allocation. They cannot increase it, and unused capacity is not redistributed. Calling `set_budget()` without arguments preserves the current configuration; it does not restore the defaults.

If you pass `budget=None` and `pools` in the same call, the old configuration is cleared before the new pool settings are applied. For example, `mc.set_budget(None, pools={"compute": 4})` still limits the compute pool but no longer sets a base budget. To disable control completely, call `mc.set_budget(None)` without `pools` and do not specify pool limits for subsequent `phase` calls.

Change the configuration between tasks. Calling `set_budget` while any phase is running raises `RuntimeError`. Unknown pool names, non-integers, and non-positive values also raise errors without partially changing the configuration.

### get_status (Inspect the Configuration and Current Phase)

```python
status = mc.get_status()
print(status)
```

| Field | Meaning |
|-------|---------|
| `budget` | Configured base budget; `None` means no base budget is set |
| `pool_overrides` | Current global pool settings |
| `phase` | Phase name in the calling context; `None` outside a phase |
| `allocations` | Current phase's pool allocations; an empty dictionary outside a phase |
| `effective_budget` | Effective budget after applying the available CPU capacity limit; `None` when no base budget is set |
| `sequential` | Whether the current phase uses sequential execution; `False` outside a phase |

This function returns the calling context's status, not a live view of all worker threads. To inspect a phase's allocations, call it inside that phase.

### phase (Choose a Phase for a Custom Workflow)

The main preprocessing and UMAP/KMeans entry points enter the appropriate phases automatically. When reading directly through a data manager, you can specify a phase explicitly so that reading uses the shared IO executor:

```python
import massflow.concurrency as mc
from massflow.data_manager import MSDataManagerImzML

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()

with mc.phase("io_read") as allocations:
    print(allocations)
    total_points = 0
    for mz, intensity, lengths, coordinates in dm.flat_generator(batch_size=256):
        total_points += intensity.size

print(total_points)
dm.close()
```

- `name`: one of the phase names in the table below. Unknown names raise `ValueError`.
- `pools`: optional limits for this phase, such as `mc.phase("compute", pools={"compute": 2})`. A non-`None` value replaces the global setting for the same pool and is then capped by the phase's allocation. `None` does not clear the global setting. Without a base budget, pool settings become the allocations directly.
- Output: `with ... as allocations` provides the phase's allocation dictionary. It is empty when control is completely disabled.

The generator must be **iterated** inside the `with` block. Creating it inside the phase and iterating after leaving the phase does not preserve the runtime limits. You can also use `phase` as a function decorator.

## How the Budget Is Allocated

### Resource Pools and Processing Phases

| Resource pool | Main use |
|---------------|----------|
| `compute` | Numba parallel computation and UMAP computation threads |
| `io_read` | The reading share of IO work, such as fragmented imzML reads |
| `io_write` | The writing share of IO work, such as temporary writes, transposition, and merging |
| `openmp` | Threads for OpenMP computation, such as KMeans |

| Phase | Used by | Pool multipliers |
|-------|---------|-------------------|
| `preprocess_pipeline` | Batch processing segments in `Preprocessor.start()` | Read 10% / compute 100% / write 10% |
| `dataset_task` | Whole-dataset tasks such as peak alignment and statistical filtering | Compute 80% / read 20% |
| `io_transfer` | Converting intermediate imzML results to ion-image storage during preprocessing | Read 20% / write 80% |
| `io_read` | Preprocessing input validation or an explicitly scoped reading workflow | Read 100% |
| `compute` | Computation entry points such as `flat_speed_process` | Compute 100% |
| `umap` | UMAP dimensionality reduction | Compute 100% |
| `kmeans` | KMeans clustering | OpenMP 100% |

`run_umap_kmeans` uses the `umap` and `kmeans` phases in sequence; `plot_umap_image` uses only `umap`. `io_write` is a resource pool name, not a phase name accepted by `phase`.

Each pool's allocation is the effective budget multiplied by its multiplier. Multipliers above 100% use 100% instead. The result is rounded down with a minimum of one worker per pool, then any configured pool limit is applied. Each pool is calculated independently, so the combined allocation can exceed the effective budget. With `set_budget(10)` and at least 10 available logical CPUs, batch preprocessing receives the following allocation by default.

```text
io_read = 1, compute = 10, io_write = 1
```

### Shared IO Capacity for Reading, Writing, and Zarr

Within a controlled phase, fragmented imzML reads, temporary writes, binary merges, H5/Zarr background prefetch, and Zarr local IO and codec work share one phase IO executor.

- When reading and writing can overlap, shared IO capacity is `io_read + io_write`. If there is only one IO pool, its allocation is used.
- In sequential mode, capacity is the larger of the read and write allocations.
- If a base budget exists but the phase has no IO allocation, 1 auxiliary IO worker is provided.

The budget-10 example above allocates 10 computation threads and 2 shared IO workers. Zarr, other reading and writing, and prefetch work draw on the same IO executor.

Zarr's asynchronous concurrency matches the shared IO capacity, and its event loop uses the same executor. On exit, the phase restores the previous Zarr settings and executor, then shuts down the executor it created. Borrowers do not shut down shared resources early.

### Sequential Execution with Small Budgets

A phase is marked as sequential when its effective budget is smaller than its number of resource pools. For example, batch preprocessing with a budget of 1–2 switches to “read one batch → compute → write”, finishing that batch before reading the next.

The compute pool then receives the full effective budget, while reading and writing each keep 1 worker and run at separate times. Adding the displayed pool allocations does not give the number of threads working simultaneously. This path still processes data in batches and does not load the whole dataset because the budget is smaller. At an effective budget of 1, H5/Zarr final writes disable background prefetch.

### CPU Capacity and Computation Libraries

The effective budget is the smaller of the configured value and available CPU capacity. Detection considers logical cores, process CPU affinity where supported (the set of CPUs the process may use), and Linux cgroup v1/v2 quotas. The default base budget uses available logical cores; it is not automatically reduced to the performance-core count.

Numba limits parallel computation through the calling thread's thread mask. The effective value is also capped by Numba's fixed thread-pool size. UMAP's `n_jobs` uses the compute allocation. When `random_state` is set, the current UMAP implementation changes `n_jobs` to 1, so increasing the budget may not improve fitting speed.

With a base budget, BLAS is limited to 1 thread. OpenMP uses its allocation in the `kmeans` phase and is limited to 1 in other phases. If the base budget is disabled but an `openmp` allocation is set, the corresponding OpenMP limit and single-thread BLAS limit still apply.

## Scope and Usage Notes

- **Multiple tasks in one process**: controlled top-level phases run one after another, each using its phase budget. Serialization applies to phases, not the entire `Preprocessor.start()` call.
- **Direct reading or writing**: some interfaces query global allocations, but calling a data manager directly does not create a phase. For the full phase limits, wrap the actual work in `phase`.
- **Restoring runtime settings**: normal and exceptional exits restore the Numba mask, native-library limits, and Zarr settings from before the phase. Native-library and Zarr settings affect the whole process. Phase admission does not manage separate processes, arbitrary external library calls, or R/Cardinal internal threads.
- **Custom threads and nested phases**: internal threads must inherit context through `contextvars.copy_context()`. Threads running Numba kernels must also apply the compute mask and finish their work before the phase exits. Only the thread that entered the phase may enter a nested phase; wait for outer background work to finish first. Worker threads with inherited context cannot enter another phase.
- **Async callers**: if an asyncio task holds a controlled phase across an `await`, other tasks on the same thread must wait for that phase to exit before entering their own.
- **Stalls and shutdown**: the pipeline can detect stalls but cannot forcibly interrupt native computation or blocking IO that is already running. Shutdown still waits for worker threads to finish.

Budget, batch size, and storage layout affect concurrency, memory use, and IO behavior separately. A larger budget does not guarantee better speed. When comparing performance, keep the data, batch size, output format, and cache conditions the same.

## Extension Interfaces

These interfaces are mainly for integrating custom execution workflows. Standard preprocessing and segmentation usually need only `set_budget` and `get_status`.

| Interface | Purpose |
|-----------|---------|
| `get_effective_budget()` | Returns the current effective base budget; `None` when no base budget is set |
| `runs_sequentially()` | Checks whether the current phase uses sequential execution |
| `get_allocation(pool, default=None)` | Looks up a pool allocation without creating a phase; unallocated auxiliary pools return 1 inside a phase with a base budget |
| `resolve_workers(pool, requested, default=...)` | Caps a request at the pool allocation; uses the request or default when no allocation exists |
| `apply_compute_mask()` | Applies the compute allocation to the calling thread and returns the effective value |
| `apply_mask(target)` / `restore_mask(previous)` | Sets the Numba mask and returns its previous value / restores the previous mask |
| `resolve_compute_workers(target)` | Caps computation threads at the available CPU capacity and Numba pool limit |
| `get_logical_cpu_count()` | Returns the logical CPU capacity available to the process |
| `detect_performance_core_workers()` | Returns a computation thread count that prefers performance cores, capped by available capacity and Numba's pool size |

`apply_compute_mask()` falls back to performance-core detection only when no compute allocation is available. Calling mask interfaces directly does not create a phase or restore the mask automatically; custom code must manage their scope.

## References

- [Normalization](./normalization.md): preprocessing examples and parameters
- [UMAP Dimensionality Reduction and Visualization](./umap.md), [KMeans Tissue Region Segmentation](./segmentation.md): segmentation workflows and result interpretation
- [`concurrency/controller.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/controller.py): budget configuration and phase allocations
- [`concurrency/runtime.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/runtime.py): shared IO executor and runtime limits
- [`concurrency/numba_bridge.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/numba_bridge.py), [`concurrency/capacity.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/capacity.py): Numba masks and CPU capacity detection
