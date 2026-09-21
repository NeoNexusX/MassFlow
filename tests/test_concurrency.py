from contextlib import nullcontext
from contextvars import copy_context
from pathlib import Path
from threading import Event, Lock, Thread
from concurrent.futures import ThreadPoolExecutor
import asyncio

import numpy as np
import zarr
from zarr.core.sync import sync

import pytest
from numba import get_num_threads, set_num_threads

from massflow.concurrency import (
    apply_compute_mask,
    get_allocation,
    get_effective_budget,
    runs_sequentially,
    get_logical_cpu_count,
    get_status,
    phase,
    resolve_compute_workers,
    resolve_workers,
    set_budget,
)
from massflow.concurrency import capacity, controller, numba_bridge, runtime
from massflow.concurrency.capacity import detect_performance_core_workers


@pytest.fixture(autouse=True)
def _restore_concurrency_state():
    prev_threads = get_num_threads()
    prev_status = get_status()
    set_budget(None)
    yield
    set_budget(None)
    if prev_status["budget"] is not None or prev_status["pool_overrides"]:
        set_budget(prev_status["budget"], pools=prev_status["pool_overrides"] or None)
    set_num_threads(prev_threads)


def _patch_runtime(monkeypatch: pytest.MonkeyPatch) -> None:
    """Neutralize numba mask application for pure allocation-math tests."""
    monkeypatch.setattr(runtime, "runtime_limits", lambda allocations, **_kwargs: nullcontext())
    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: 64)


# ---------------------------------------------------------------------------
# Controller: budget and allocation math
# ---------------------------------------------------------------------------


def test_set_budget_roundtrip() -> None:
    set_budget(4)

    status = get_status()
    assert status["budget"] == 4
    assert status["phase"] is None
    assert status["allocations"] == {}


def test_set_budget_validates_type() -> None:
    with pytest.raises(TypeError):
        set_budget(1.5)  # type: ignore[arg-type]


def test_set_budget_validates_value() -> None:
    with pytest.raises(ValueError):
        set_budget(0)


def test_set_budget_none_resets_overrides() -> None:
    set_budget(4, pools={"compute": 2})

    set_budget(None)

    status = get_status()
    assert status["budget"] is None
    assert status["pool_overrides"] == {}


def test_pool_override_validates_key_and_value() -> None:
    with pytest.raises(ValueError):
        set_budget(4, pools={"not_a_pool": 2})

    with pytest.raises(ValueError):
        set_budget(4, pools={"compute": 0})


def test_initial_budget_uses_environment_or_detection(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: 7)
    monkeypatch.setenv("MASSFLOW_BUDGET", "6")
    assert controller._read_initial_budget() == 6

    monkeypatch.setenv("MASSFLOW_BUDGET", "not-a-number")
    assert controller._read_initial_budget() == 7

    monkeypatch.setenv("MASSFLOW_BUDGET", "0")
    assert controller._read_initial_budget() == 7

    monkeypatch.delenv("MASSFLOW_BUDGET")
    assert controller._read_initial_budget() == 7


def test_get_allocation_returns_default_when_disabled() -> None:
    assert get_allocation("compute", default=-1) == -1
    assert get_allocation("io_read", default=0) == 0
    assert get_allocation("openmp") is None


def test_get_allocation_rejects_unknown_pool() -> None:
    with pytest.raises(ValueError):
        get_allocation("not_a_pool")


def test_phase_allocation_math(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_runtime(monkeypatch)
    set_budget(10)

    with phase("preprocess_pipeline") as allocations:
        assert allocations == {"io_read": 2, "compute": 6, "io_write": 2}
        assert get_allocation("io_read") == 2

    with phase("dataset_task") as allocations:
        assert allocations == {"compute": 8, "io_read": 2}

    with phase("umap") as allocations:
        assert allocations == {"compute": 10}

    with phase("kmeans") as allocations:
        assert allocations == {"openmp": 10}

    assert get_status()["phase"] is None
    assert get_status()["allocations"] == {}


def test_phase_allocation_uses_full_budget(monkeypatch: pytest.MonkeyPatch) -> None:
    """Remainder workers are handed out; the budget is not truncated away."""
    _patch_runtime(monkeypatch)
    set_budget(8)

    with phase("preprocess_pipeline") as allocations:
        assert sum(allocations.values()) == 8
        assert allocations["compute"] == 5

    with phase("dataset_task") as allocations:
        assert sum(allocations.values()) == 8


@pytest.mark.parametrize("budget", [1, 2, 3, 8])
@pytest.mark.parametrize("cpu_count", [1, 2, 8])
@pytest.mark.parametrize(
    "phase_name", list(controller._PHASE_RATIOS)
)  # pylint: disable=protected-access
def test_all_phases_adapt_to_available_budget(
    monkeypatch: pytest.MonkeyPatch,
    budget: int,
    cpu_count: int,
    phase_name: str,
) -> None:
    _patch_runtime(monkeypatch)
    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: cpu_count)
    set_budget(budget)
    effective = min(budget, cpu_count)
    with phase(phase_name) as allocations:
        assert get_effective_budget() == effective
        assert all(1 <= value <= effective for value in allocations.values())
        assert runs_sequentially() == (effective < len(allocations))
        if not runs_sequentially():
            assert sum(allocations.values()) <= effective
    assert get_status()["phase"] is None
    assert not runs_sequentially()


def test_phase_requires_known_name(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_runtime(monkeypatch)
    set_budget(4)

    with pytest.raises(ValueError):
        with phase("not_a_phase"):
            pass


def test_phase_is_noop_when_disabled() -> None:
    prev_threads = get_num_threads()

    with phase("umap") as allocations:
        assert allocations == {}
        assert get_allocation("compute", default=-1) == -1

    assert get_num_threads() == prev_threads


def test_pool_override_beats_phase_allocation(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_runtime(monkeypatch)
    set_budget(10, pools={"compute": 3})

    with phase("umap"):
        assert get_allocation("compute") == 3

    set_budget(10, pools={"compute": None})

    with phase("umap"):
        assert get_allocation("compute") == 10


def test_nested_phase_restores_previous_allocations(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_runtime(monkeypatch)
    set_budget(10)

    with phase("umap"):
        with phase("kmeans"):
            assert get_allocation("openmp") == 10
            assert get_allocation("compute") == 1
        assert get_status()["phase"] == "umap"
        assert get_allocation("compute") == 10
    assert get_status()["phase"] is None


def test_set_budget_rejects_active_reconfiguration(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_runtime(monkeypatch)
    set_budget(10)

    with phase("umap"):
        with pytest.raises(RuntimeError, match="active phase"):
            set_budget(5)
        assert get_allocation("compute") == 10


# ---------------------------------------------------------------------------
# Numba bridge: mask application
# ---------------------------------------------------------------------------


def test_resolve_compute_workers_clamps_with_warning(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # NUMBA_NUM_THREADS is populated dynamically by Numba.
    # pylint: disable-next=no-member
    limit = min(get_logical_cpu_count(), capacity.config.NUMBA_NUM_THREADS)  # type: ignore[attr-defined]
    warnings: list[str] = []
    monkeypatch.setattr(
        numba_bridge.logger,
        "warning",
        lambda msg, *args: warnings.append(msg % args),
    )

    assert resolve_compute_workers(10**9) == limit
    assert warnings and "exceeds the available limit" in warnings[0]


def test_resolve_compute_workers_validates() -> None:
    with pytest.raises(TypeError):
        resolve_compute_workers(1.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        resolve_compute_workers(0)


def test_apply_compute_mask_uses_phase_allocation() -> None:
    set_budget(2)
    expected = min(2, get_logical_cpu_count(), capacity.config.NUMBA_NUM_THREADS)  # type: ignore[attr-defined]  # pylint: disable=no-member

    with phase("umap"):
        applied = apply_compute_mask()
        assert applied == expected
        assert get_num_threads() == expected


def test_apply_compute_mask_falls_back_to_detection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(numba_bridge, "detect_performance_core_workers", lambda: 1)

    assert apply_compute_mask() == 1
    assert get_num_threads() == 1


# ---------------------------------------------------------------------------
# Capacity detection ( migrated from test_numba_runtime )
# ---------------------------------------------------------------------------


def test_get_logical_cpu_count_returns_positive_int() -> None:
    count = get_logical_cpu_count()
    print(f"Logical CPU count: {count}")
    assert isinstance(count, int)
    assert count >= 1


def test_get_logical_cpu_count_respects_affinity_and_cgroup(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(capacity.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(capacity.os, "sched_getaffinity", lambda _pid: set(range(8)), raising=False)
    monkeypatch.setattr(capacity.platform, "system", lambda: "Linux")
    monkeypatch.setattr(capacity, "_detect_linux_cgroup_cpu_count", lambda: 2)

    assert get_logical_cpu_count() == 2


def test_detect_linux_cgroup_v2_cpu_count(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    cgroup_root = tmp_path / "cgroup"
    current_cgroup = cgroup_root / "user.slice" / "app.scope"
    current_cgroup.mkdir(parents=True)
    (cgroup_root / "cpu.max").write_text("max 100000\n", encoding="utf-8")
    (cgroup_root / "user.slice" / "cpu.max").write_text("100000 100000\n", encoding="utf-8")
    (current_cgroup / "cpu.max").write_text("200000 100000\n", encoding="utf-8")

    proc_cgroup = tmp_path / "self.cgroup"
    proc_cgroup.write_text("0::/user.slice/app.scope\n", encoding="utf-8")
    proc_mountinfo = tmp_path / "mountinfo"
    proc_mountinfo.write_text(
        f"36 25 0:32 / {cgroup_root} rw - cgroup2 cgroup rw\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(capacity, "_PROC_SELF_CGROUP", proc_cgroup)
    monkeypatch.setattr(capacity, "_PROC_SELF_MOUNTINFO", proc_mountinfo)

    assert capacity._detect_linux_cgroup_cpu_count() == 1


def test_detect_linux_cgroup_v2_namespace_root(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    cgroup_root = tmp_path / "cgroup"
    current_cgroup = cgroup_root / "child.scope"
    current_cgroup.mkdir(parents=True)
    (cgroup_root / "cpu.max").write_text("max 100000\n", encoding="utf-8")
    (current_cgroup / "cpu.max").write_text("150000 100000\n", encoding="utf-8")

    proc_cgroup = tmp_path / "self.cgroup"
    proc_cgroup.write_text("0::/child.scope\n", encoding="utf-8")
    proc_mountinfo = tmp_path / "mountinfo"
    proc_mountinfo.write_text(
        f"36 25 0:32 / {cgroup_root} rw - cgroup2 cgroup rw\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(capacity, "_PROC_SELF_CGROUP", proc_cgroup)
    monkeypatch.setattr(capacity, "_PROC_SELF_MOUNTINFO", proc_mountinfo)

    assert capacity._detect_linux_cgroup_cpu_count() == 2


def test_detect_linux_cgroup_v1_cpu_count(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    cgroup_root = tmp_path / "cpu,cpuacct"
    current_cgroup = cgroup_root / "system.slice" / "massflow.service"
    current_cgroup.mkdir(parents=True)
    for cgroup_dir, quota in (
        (cgroup_root, "-1"),
        (cgroup_root / "system.slice", "200000"),
        (current_cgroup, "300000"),
    ):
        (cgroup_dir / "cpu.cfs_quota_us").write_text(f"{quota}\n", encoding="utf-8")
        (cgroup_dir / "cpu.cfs_period_us").write_text("100000\n", encoding="utf-8")

    proc_cgroup = tmp_path / "self.cgroup"
    proc_cgroup.write_text("2:cpu,cpuacct:/system.slice/massflow.service\n", encoding="utf-8")
    proc_mountinfo = tmp_path / "mountinfo"
    proc_mountinfo.write_text(
        f"31 25 0:28 / {cgroup_root} rw - cgroup cgroup rw,cpu,cpuacct\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(capacity, "_PROC_SELF_CGROUP", proc_cgroup)
    monkeypatch.setattr(capacity, "_PROC_SELF_MOUNTINFO", proc_mountinfo)

    assert capacity._detect_linux_cgroup_cpu_count() == 2


def test_detect_linux_cgroup_cpu_count_uses_fixed_root_fallback(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    cgroup_root = tmp_path / "cgroup"
    cgroup_root.mkdir()
    (cgroup_root / "cpu.max").write_text("150000 100000\n", encoding="utf-8")

    monkeypatch.setattr(capacity, "_PROC_SELF_CGROUP", tmp_path / "missing")
    monkeypatch.setattr(capacity, "_PROC_SELF_MOUNTINFO", tmp_path / "missing")
    monkeypatch.setattr(capacity, "_CGROUP_V2_FALLBACK_ROOT", cgroup_root)
    monkeypatch.setattr(capacity, "_CGROUP_V1_FALLBACK_ROOTS", ())

    assert capacity._detect_linux_cgroup_cpu_count() == 2


def test_detect_performance_core_workers_returns_positive_int() -> None:
    count = detect_performance_core_workers()
    assert isinstance(count, int)
    assert count >= 1


def test_detect_performance_core_workers_windows_branch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(capacity.platform, "system", lambda: "Windows")
    monkeypatch.setattr(capacity, "_detect_windows_intel_performance_workers", lambda: 6)
    monkeypatch.setattr(capacity, "get_logical_cpu_count", lambda: 8)
    monkeypatch.setattr(capacity.config, "NUMBA_NUM_THREADS", 8)

    assert detect_performance_core_workers() == 6


def test_detect_performance_core_workers_fallback_minimum_one(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(capacity.platform, "system", lambda: "Unknown")
    monkeypatch.setattr(capacity, "get_logical_cpu_count", lambda: 1)

    assert detect_performance_core_workers() == 1


def test_detect_performance_core_workers_fallback_uses_available_workers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(capacity.platform, "system", lambda: "Unknown")
    monkeypatch.setattr(capacity, "get_logical_cpu_count", lambda: 8)
    monkeypatch.setattr(capacity.config, "NUMBA_NUM_THREADS", 8)

    assert detect_performance_core_workers() == 8


def test_detect_performance_core_workers_caps_numba_thread_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(capacity.platform, "system", lambda: "Unknown")
    monkeypatch.setattr(capacity, "get_logical_cpu_count", lambda: 64)
    monkeypatch.setattr(capacity.config, "NUMBA_NUM_THREADS", 2)

    assert detect_performance_core_workers() == 2


# ---------------------------------------------------------------------------
# Admission, worker propagation and backend integration
# ---------------------------------------------------------------------------


def test_pool_update_preserves_budget_and_invalid_update_is_atomic() -> None:
    set_budget(8)
    set_budget(pools={"io_read": 2})
    assert get_status()["budget"] == 8
    previous = get_status()
    with pytest.raises(ValueError):
        set_budget(pools={"io_write": 3, "unknown": 1})
    assert get_status() == previous
    set_budget(None)
    assert get_status()["pool_overrides"] == {}


def test_budget_applies_outside_phase() -> None:
    set_budget(1)
    assert resolve_workers("compute", 8, default=8) == 1
    assert resolve_workers("io_read", None, default=8) == 1


def test_phase_override_controls_dataset_mask() -> None:
    previous = get_num_threads()
    set_budget(8, pools={"compute": 1})
    with phase("dataset_task"):
        assert get_num_threads() == get_allocation("compute") == 1
    assert get_num_threads() == previous


def test_worker_context_sees_owner_allocation(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_runtime(monkeypatch)
    set_budget(8)
    with phase("preprocess_pipeline"):
        with ThreadPoolExecutor(max_workers=1) as executor:
            result = executor.submit(copy_context().run, get_status).result(timeout=5)
        assert result["allocations"] == get_status()["allocations"]
        assert result["allocations"]["io_write"] == 1


def test_parallel_phases_wait_for_admission(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_runtime(monkeypatch)
    set_budget(8)
    waiting = Event()
    entered = Event()
    statuses: list[dict] = []

    def run_second() -> None:
        waiting.set()
        with phase("kmeans"):
            statuses.append(dict(get_status()))
            entered.set()

    with phase("preprocess_pipeline"):
        worker = Thread(target=run_second)
        worker.start()
        assert waiting.wait(timeout=5)
        assert not entered.wait(timeout=0.05)
        assert get_status()["phase"] == "preprocess_pipeline"
    worker.join(timeout=5)
    assert not worker.is_alive()
    assert entered.is_set()
    assert statuses[0]["phase"] == "kmeans"
    assert get_status()["phase"] is None


def test_exception_restores_phase_and_mask() -> None:
    previous = get_num_threads()
    set_budget(1)
    with pytest.raises(RuntimeError, match="task failed"):
        with phase("compute"):
            assert get_num_threads() == 1
            raise RuntimeError("task failed")
    assert get_status()["phase"] is None
    assert get_num_threads() == previous
    set_budget(None)


def test_flat_speed_process_uses_budget() -> None:
    from massflow.tools.dm_process import flat_speed_process

    seen: list[int] = []
    set_budget(1)
    flat_speed_process([None], lambda **_kwargs: seen.append(get_num_threads()), {})
    assert seen == [1]


def test_read_pool_tracks_phase_allocation(tmp_path: Path) -> None:
    from massflow.data_manager import MSDataManagerImzML

    manager = MSDataManagerImzML(temp_dir=str(tmp_path))
    try:
        assert manager.read_workers == 8
        set_budget(1)
        with phase("io_read"):
            assert manager.read_workers == 1
            manager.threads_executor.submit(lambda: None).result(timeout=5)
        set_budget(None)
        assert manager.read_workers == 8
    finally:
        manager.close()


@pytest.mark.parametrize("budget", [1, 8])
@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize("fail", [False, True])
def test_zarr_restores_existing_executor_and_config(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    budget: int,
    nested: bool,
    fail: bool,
) -> None:
    async def get_executor():
        loop = asyncio.get_running_loop()
        await asyncio.to_thread(lambda: None)
        return loop._default_executor  # type: ignore[attr-defined]  # pylint: disable=protected-access

    original_executor = sync(get_executor())
    original_concurrency = zarr.config.get("async.concurrency")
    original_workers = zarr.config.get("threading.max_workers")
    original_mask = get_num_threads()
    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: 8)
    set_budget(budget)
    with pytest.raises(RuntimeError, match="task failed") if fail else nullcontext():
        with phase("io_read"):
            outer_executor = runtime.get_io_executor()
            assert outer_executor is not None
            assert outer_executor is not original_executor
            assert zarr.config.get("async.concurrency") == budget
            assert sync(get_executor()) is outer_executor
            array = zarr.create_array(
                str(tmp_path / "data.zarr"), shape=(16,), chunks=(2,), dtype="int32"
            )
            array[:] = np.arange(16)
            np.testing.assert_array_equal(array[:], np.arange(16))
            if nested:
                with phase("compute"):
                    inner_executor = runtime.get_io_executor()
                    assert inner_executor is not None
                    assert inner_executor is not outer_executor
                    assert sync(get_executor()) is inner_executor
                    assert zarr.config.get("async.concurrency") == 1
                    np.testing.assert_array_equal(array[:], np.arange(16))
                assert runtime.get_io_executor() is outer_executor
                assert sync(get_executor()) is outer_executor
                assert zarr.config.get("async.concurrency") == budget
            if fail:
                raise RuntimeError("task failed")
    assert sync(get_executor()) is original_executor
    assert zarr.config.get("async.concurrency") == original_concurrency
    assert zarr.config.get("threading.max_workers") == original_workers
    assert get_num_threads() == original_mask
    assert runtime.get_io_executor() is None


def test_inherited_worker_cannot_deadlock_on_nested_admission(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_runtime(monkeypatch)
    set_budget(8)

    def enter_phase() -> None:
        with phase("compute"):
            pytest.fail("Worker acquired owner-only admission")

    with phase("preprocess_pipeline"):
        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(copy_context().run, enter_phase)
            with pytest.raises(RuntimeError, match="only its owner"):
                future.result(timeout=5)


@pytest.mark.parametrize("cpu_count,budget", [(1, 1), (2, 8), (8, 8)])
@pytest.mark.parametrize(
    "output_format", ["imzml", "spectra_zarr", "ion_image_binary", "ion_image_h5", "ion_image_zarr"]
)
def test_pipeline_adapts_without_changing_results(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    cpu_count: int,
    budget: int,
    output_format: str,
) -> None:
    from typing import cast

    from massflow.data_manager import MSDataManagerImzML, MSDataManagerZarr
    from massflow.preprocess import Preprocessor
    from massflow.preprocess.flat_pre_fun import FlatBatchResult
    from massflow.preprocess.pipeline import flat_executor
    from massflow.preprocess.pipeline.types import OutputFormat
    from massflow.tools.stream_imzml_writer import ImzMLWriter

    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: cpu_count)
    mz_axis = np.arange(100, 132, dtype=np.float64)
    source_matrix = np.arange(128, dtype=np.float32).reshape(4, 32) + 1
    source_path = tmp_path / "source.imzML"
    with ImzMLWriter(str(source_path), mode="continuous") as writer:
        for index, intensities in enumerate(source_matrix):
            writer.add_spectrum(mz_axis, intensities, (index + 1, 1, 1))

    source = MSDataManagerImzML(filepath=str(source_path))
    source.load_head_data()
    started: list[str] = []

    def make_thread(*args, **kwargs):
        started.append(kwargs["name"])
        return Thread(*args, **kwargs)

    def scale(mz_data, intensity, lengths) -> FlatBatchResult:
        return FlatBatchResult(mz_data, intensity * 2, lengths)

    def offset(mz_data, intensity, lengths) -> FlatBatchResult:
        return FlatBatchResult(mz_data, intensity + 3, lengths)

    monkeypatch.setattr(flat_executor, "Thread", make_thread)
    if min(cpu_count, budget) == 1:
        from massflow.data_manager import ion_image_data_manager_h5
        from massflow.msi_zarr import writer as zarr_writer

        def reject_prefetch(*_args, **_kwargs):
            pytest.fail("A single-worker output must not start a background prefetch pool")

        monkeypatch.setattr(ion_image_data_manager_h5, "io_executor", reject_prefetch)
        monkeypatch.setattr(zarr_writer, "io_executor", reject_prefetch)
    set_budget(budget)
    pipeline = Preprocessor(
        source,
        batch_size=2,
        temp_dir=str(tmp_path),
        output_format=cast(OutputFormat, output_format),
        output_path=str(tmp_path / "result"),
    )
    # Two ordered tasks exercise the shared batch processor in both modes.
    pipeline._register_task("scale", apply_fn=scale)  # pylint: disable=protected-access
    pipeline._register_task("offset", apply_fn=offset)  # pylint: disable=protected-access
    result = None
    try:
        result = pipeline.start()
        expected = source_matrix * 2 + 3
        assert len(started) == (0 if min(cpu_count, budget) < 3 else 3)
        if isinstance(result, (MSDataManagerImzML, MSDataManagerZarr)):
            actual = np.concatenate(
                [intensity for _, intensity, _, _ in result.flat_generator(batch_size=2)]
            ).reshape(4, 32)
        else:
            actual = np.column_stack([result.get_ion_vector(index) for index in range(32)])
        np.testing.assert_array_equal(actual, expected)
    finally:
        if result is not None:
            result.close()
        source.close()


@pytest.mark.parametrize("budget", [1, 2])
def test_low_budget_write_finishes_before_returning(tmp_path: Path, budget: int) -> None:
    from massflow.tools.ion_major_temp_store import IonMajorTempStore

    # A store can outlive a phase; the current allocation must govern backpressure.
    with IonMajorTempStore(temp_dir=tmp_path, dtype=np.float32, workers=8) as store:
        set_budget(budget)
        with phase("preprocess_pipeline"):
            store.append_batch(np.arange(6, dtype=np.float32), n_pixels=2, n_ions=3)
            assert store.workers == 8
            assert not store._pending  # pylint: disable=protected-access
            batches = store.finish()
            np.testing.assert_array_equal(
                store.read_ion_block(batches[0], 0, 3),
                np.array([[0, 3], [1, 4], [2, 5]], dtype=np.float32),
            )


@pytest.mark.parametrize("budget", [1, 2])
def test_sequential_failure_closes_output_and_restores_context(budget: int) -> None:
    from unittest.mock import create_autospec

    from massflow.data_manager import MSDataManager
    from massflow.preprocess.pipeline import flat_executor
    from massflow.preprocess.pipeline.types import InputMetaData, PreprocessTask

    source = create_autospec(MSDataManager, instance=True)
    sink = create_autospec(flat_executor.FlatOutputSink, instance=True)
    metadata = InputMetaData(
        batch_size=1,
        temp_dir=None,
        queue_ab_size=1,
        queue_bc_size=1,
        result_path=None,
        result_format="imzml",
    )

    def batches():
        yield np.array([100.0]), np.array([2.0]), np.array([1]), np.array([[1, 1]])
        pytest.fail("The next batch was read after a processing failure")

    def fail(**_kwargs):
        raise ValueError("processing failed")

    source.flat_generator.return_value = batches()
    previous_mask = get_num_threads()
    set_budget(budget)
    with pytest.raises(ValueError, match="processing failed"):
        with phase("preprocess_pipeline"):
            flat_executor._run_sequential_flat_task(  # pylint: disable=protected-access
                source,
                [PreprocessTask(name="fail", apply_fn=fail)],
                metadata,
                sink,
                2,
            )
    sink.close.assert_called_once()
    sink.write_chunk.assert_not_called()
    sink.finalize.assert_not_called()
    assert get_status()["phase"] is None
    assert get_num_threads() == previous_mask


# ---------------------------------------------------------------------------
# Shared IO executor ownership and scheduling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "phase_name,budget,expected_workers",
    [
        ("preprocess_pipeline", 1, 1),
        ("preprocess_pipeline", 2, 1),
        ("preprocess_pipeline", 8, 3),
        ("io_transfer", 1, 1),
        ("io_transfer", 2, 2),
        ("io_transfer", 8, 8),
    ],
)
def test_shared_io_executor_identity_and_peak(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    phase_name: str,
    budget: int,
    expected_workers: int,
) -> None:
    from massflow.data_manager import MSDataManagerImzML

    async def zarr_executor():
        return asyncio.get_running_loop()._default_executor  # type: ignore[attr-defined]  # pylint: disable=protected-access

    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: 8)
    set_budget(budget)
    release = Event()
    saturated = Event()
    lock = Lock()
    active = 0
    peak = 0

    def work() -> None:
        nonlocal active, peak
        with lock:
            active += 1
            peak = max(peak, active)
            if active == expected_workers:
                saturated.set()
        try:
            assert release.wait(timeout=5)
        finally:
            with lock:
                active -= 1

    manager = MSDataManagerImzML(temp_dir=str(tmp_path))
    try:
        with phase(phase_name):
            shared = runtime.get_io_executor()
            assert shared is not None
            assert manager.threads_executor is shared
            assert sync(zarr_executor()) is shared
            assert shared._max_workers == expected_workers  # pylint: disable=protected-access
            with ThreadPoolExecutor(max_workers=1) as caller:
                inherited = caller.submit(copy_context().run, runtime.get_io_executor)
                assert inherited.result(timeout=5) is shared
            with runtime.io_executor(max_workers=99) as borrowed:
                try:
                    futures = [borrowed.submit(work) for _ in range(expected_workers * 2)]
                    assert saturated.wait(timeout=5)
                    with lock:
                        assert active == expected_workers
                finally:
                    release.set()
                for future in futures:
                    future.result(timeout=5)
            assert peak == expected_workers
            assert shared.submit(lambda: "still open").result(timeout=5) == "still open"
        assert runtime.get_io_executor() is None
    finally:
        release.set()
        manager.close()


def test_borrowed_io_scope_drains_on_exception_without_closing_shared_pool() -> None:
    started = Event()
    raising = Event()
    release = Event()
    finished = Event()
    exited = Event()
    errors: list[BaseException] = []

    def work() -> None:
        started.set()
        assert release.wait(timeout=5)
        finished.set()

    def caller() -> None:
        try:
            with runtime.io_executor(max_workers=1) as executor:
                executor.submit(work)
                raising.set()
                raise RuntimeError("scope failed")
        except BaseException as exc:
            errors.append(exc)
        finally:
            exited.set()

    set_budget(1)
    with phase("io_read"):
        shared = runtime.get_io_executor()
        assert shared is not None
        thread = Thread(target=copy_context().run, args=(caller,))
        thread.start()
        try:
            assert started.wait(timeout=5)
            assert raising.wait(timeout=5)
            assert not exited.wait(timeout=0.05)
        finally:
            release.set()
            thread.join(timeout=5)
        assert not thread.is_alive()
        assert finished.is_set()
        assert len(errors) == 1
        assert isinstance(errors[0], RuntimeError)
        assert str(errors[0]) == "scope failed"
        assert shared.submit(lambda: 42).result(timeout=5) == 42


def test_uncontrolled_io_scope_owns_its_executor() -> None:
    assert runtime.get_io_executor() is None
    with runtime.io_executor(max_workers=1) as executor:
        assert executor.submit(lambda: 42).result(timeout=5) == 42
    with pytest.raises(RuntimeError, match="shutdown"):
        executor.submit(lambda: None)


@pytest.mark.parametrize("budget", [1, 2, 8])
def test_closing_io_borrowers_keeps_shared_executor_usable(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, budget: int
) -> None:
    from massflow.data_manager import MSDataManagerImzML
    from massflow.tools.ion_major_temp_store import IonMajorTempStore

    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: 8)
    manager = MSDataManagerImzML(temp_dir=str(tmp_path))
    set_budget(budget)
    try:
        with phase("io_transfer"):
            shared = runtime.get_io_executor()
            assert shared is not None
            assert manager.threads_executor is shared
            submitted: list[str] = []
            original_submit = shared.submit

            def record_submit(fn, *args, **kwargs):
                submitted.append(getattr(fn, "__name__", ""))
                return original_submit(fn, *args, **kwargs)

            monkeypatch.setattr(shared, "submit", record_submit)
            with IonMajorTempStore(temp_dir=tmp_path, dtype=np.float32, workers=8) as store:
                store.append_batch(np.arange(6, dtype=np.float32), n_pixels=2, n_ions=3)
                batches = store.finish()
                np.testing.assert_array_equal(
                    store.read_ion_block(batches[0], 0, 3),
                    np.array([[0, 3], [1, 4], [2, 5]], dtype=np.float32),
                )
            assert "_write_batch" in submitted
            assert not batches[0].path.exists()
            manager.close()
            assert shared.submit(lambda: 42).result(timeout=5) == 42
            array = zarr.create_array(
                str(tmp_path / "after_close.zarr"), shape=(4,), chunks=(2,), dtype="int32"
            )
            array[:] = np.arange(4)
            np.testing.assert_array_equal(array[:], np.arange(4))
    finally:
        manager.close()


def test_temp_store_can_append_in_successive_phases(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from massflow.tools.ion_major_temp_store import IonMajorTempStore

    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: 8)
    executors: list[ThreadPoolExecutor] = []
    with IonMajorTempStore(temp_dir=tmp_path, dtype=np.float32, workers=8) as store:
        for index, budget in enumerate((1, 2, 8)):
            set_budget(budget)
            with phase("io_transfer"):
                shared = runtime.get_io_executor()
                assert shared is not None
                assert shared not in executors
                executors.append(shared)
                submitted: list[str] = []
                original_submit = shared.submit

                def record_submit(fn, *args, **kwargs):
                    submitted.append(getattr(fn, "__name__", ""))
                    return original_submit(fn, *args, **kwargs)

                with monkeypatch.context() as patch:
                    patch.setattr(shared, "submit", record_submit)
                    store.append_batch(
                        np.arange(6, dtype=np.float32) + index * 6,
                        n_pixels=2,
                        n_ions=3,
                    )
                assert "_write_batch" in submitted
        batches = store.finish()
        assert len(batches) == 3
        for index, batch in enumerate(batches):
            np.testing.assert_array_equal(
                store.read_ion_block(batch, 0, 3),
                (np.arange(6, dtype=np.float32) + index * 6).reshape(2, 3).T,
            )


@pytest.mark.parametrize("fail", [False, True])
def test_zarr_restores_initially_absent_executor(tmp_path: Path, fail: bool) -> None:
    async def current_executor():
        return asyncio.get_running_loop()._default_executor  # type: ignore[attr-defined]  # pylint: disable=protected-access

    with zarr.config.set({"threading.max_workers": None}):
        original = sync(runtime._swap_zarr_executor(None))  # pylint: disable=protected-access
        try:
            set_budget(1)
            with pytest.raises(RuntimeError, match="task failed") if fail else nullcontext():
                with phase("io_read"):
                    assert sync(current_executor()) is runtime.get_io_executor()
                    array = zarr.create_array(
                        str(tmp_path / "initially_absent.zarr"),
                        shape=(4,),
                        chunks=(2,),
                        dtype="int32",
                    )
                    array[:] = np.arange(4)
                    np.testing.assert_array_equal(array[:], np.arange(4))
                    if fail:
                        raise RuntimeError("task failed")
            assert sync(current_executor()) is None
        finally:
            sync(runtime._swap_zarr_executor(original))  # pylint: disable=protected-access


@pytest.mark.parametrize("budget", [1, 2, 8])
@pytest.mark.parametrize("backend", ["h5", "zarr"])
def test_multiband_prefetch_uses_shared_executor_and_preserves_values(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, budget: int, backend: str
) -> None:
    from massflow.data_manager import IonImageDataManagerH5
    from massflow.msi_zarr import MSIZarrReader, MSIZarrWriter, ZARR_MODE_CONTINUOUS

    monkeypatch.setattr(controller, "get_logical_cpu_count", lambda: 8)
    matrix = np.arange(4 * 48, dtype=np.float32).reshape(4, 48)
    mz_axis = np.arange(100, 148, dtype=np.float64)
    coordinates = np.column_stack((np.arange(1, 5), np.ones(4, dtype=np.int64)))
    path = tmp_path / f"multiband.{backend}"
    set_budget(budget)
    with phase("io_transfer"):
        shared = runtime.get_io_executor()
        assert shared is not None
        submitted: list[str] = []
        original_submit = shared.submit

        def record_submit(fn, *args, **kwargs):
            submitted.append(getattr(fn, "__name__", ""))
            return original_submit(fn, *args, **kwargs)

        monkeypatch.setattr(shared, "submit", record_submit)
        writer: IonImageDataManagerH5 | MSIZarrWriter
        if backend == "h5":
            writer = IonImageDataManagerH5(
                path, mode="w", width=4, height=1, temp_dir=tmp_path, ion_block=16
            )
            prefetch_task = "_read_scatter"
        else:
            writer = MSIZarrWriter(
                path,
                row_axis="ion",
                encoding=ZARR_MODE_CONTINUOUS,
                metadata=None,
                width=4,
                height=1,
                temp_dir=tmp_path,
                ion_block=16,
                target_chunk_bytes=16 * 4 * matrix.dtype.itemsize,
            )
            prefetch_task = "assemble_band"
        try:
            for start in (0, 2):
                writer.write_flat_batch(
                    intensity_flat=matrix[start:start + 2].ravel(),
                    lengths=np.full(2, 48, dtype=np.uint32),
                    coordinates=coordinates[start:start + 2],
                    mz_axis=mz_axis,
                )
            writer.finalize()
            assert submitted.count(prefetch_task) == (0 if budget == 1 else 2)
            if isinstance(writer, IonImageDataManagerH5):
                actual = writer.get_ion_block(0, 48)
            else:
                reader = MSIZarrReader(path)
                try:
                    actual = reader.read_rows(0, 48)
                finally:
                    reader.close()
            np.testing.assert_array_equal(actual, matrix.T)
            assert shared.submit(lambda: 42).result(timeout=5) == 42
        finally:
            writer.close()
