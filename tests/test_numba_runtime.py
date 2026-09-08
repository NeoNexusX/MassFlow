from pathlib import Path

import pytest
from numba import get_num_threads, set_num_threads
from massflow.preprocess.numba import numba_runtime

from massflow.preprocess.numba.numba_runtime import (
    apply_numba_runtime,
    detect_performance_core_workers,
    get_global_numba_runtime,
    get_logical_cpu_count,
    set_global_numba_runtime,
)


@pytest.fixture(autouse=True)
def _restore_numba_runtime_state():
    prev_runtime_threads = get_num_threads()
    prev_global_runtime = get_global_numba_runtime()
    yield
    set_global_numba_runtime(prev_global_runtime["max_workers"])
    set_num_threads(prev_runtime_threads)


def test_set_global_numba_runtime_roundtrip() -> None:
    target = 1

    configured = set_global_numba_runtime(target)

    assert configured["max_workers"] == target
    assert get_global_numba_runtime()["max_workers"] == target
    assert get_num_threads() == target


def test_apply_numba_runtime_uses_global_when_not_provided() -> None:
    target = 1

    set_global_numba_runtime(target)
    applied = apply_numba_runtime()

    assert applied == target
    assert get_num_threads() == target


def test_apply_numba_runtime_prefers_explicit_override() -> None:
    prev_threads = get_num_threads()
    supports_two_threads = True
    try:
        set_num_threads(2)
    except ValueError:
        supports_two_threads = False
    finally:
        set_num_threads(prev_threads)

    if not supports_two_threads:
        pytest.skip("Need at least 2 threads to verify explicit override precedence")

    set_global_numba_runtime(1)
    applied = apply_numba_runtime(2)

    assert applied == 2
    assert get_num_threads() == 2


def test_set_global_numba_runtime_validates_range() -> None:
    with pytest.raises(ValueError):
        set_global_numba_runtime(0)

    with pytest.raises(ValueError):
        set_global_numba_runtime(10**9)


def test_set_global_numba_runtime_validates_type() -> None:
    with pytest.raises(TypeError):
        set_global_numba_runtime(1.5)  # type: ignore[arg-type]


def test_get_logical_cpu_count_returns_positive_int() -> None:
    count = get_logical_cpu_count()
    print(f"Logical CPU count: {count}")
    assert isinstance(count, int)
    assert count >= 1


def test_get_logical_cpu_count_respects_affinity_and_cgroup(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(numba_runtime.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(
        numba_runtime.os, "sched_getaffinity", lambda _pid: set(range(8)), raising=False
    )
    monkeypatch.setattr(numba_runtime.platform, "system", lambda: "Linux")
    monkeypatch.setattr(numba_runtime, "_detect_linux_cgroup_cpu_count", lambda: 2)

    assert get_logical_cpu_count() == 2


def test_detect_linux_cgroup_v2_cpu_count(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    cgroup_root = tmp_path / "cgroup"
    current_cgroup = cgroup_root / "user.slice" / "app.scope"
    current_cgroup.mkdir(parents=True)
    (cgroup_root / "cpu.max").write_text("max 100000\n", encoding="utf-8")
    (cgroup_root / "user.slice" / "cpu.max").write_text(
        "100000 100000\n", encoding="utf-8"
    )
    (current_cgroup / "cpu.max").write_text("200000 100000\n", encoding="utf-8")

    proc_cgroup = tmp_path / "self.cgroup"
    proc_cgroup.write_text("0::/user.slice/app.scope\n", encoding="utf-8")
    proc_mountinfo = tmp_path / "mountinfo"
    proc_mountinfo.write_text(
        f"36 25 0:32 / {cgroup_root} rw - cgroup2 cgroup rw\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(numba_runtime, "_PROC_SELF_CGROUP", proc_cgroup)
    monkeypatch.setattr(numba_runtime, "_PROC_SELF_MOUNTINFO", proc_mountinfo)

    assert numba_runtime._detect_linux_cgroup_cpu_count() == 1


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
    monkeypatch.setattr(numba_runtime, "_PROC_SELF_CGROUP", proc_cgroup)
    monkeypatch.setattr(numba_runtime, "_PROC_SELF_MOUNTINFO", proc_mountinfo)

    assert numba_runtime._detect_linux_cgroup_cpu_count() == 2


def test_detect_linux_cgroup_v1_cpu_count(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
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
    proc_cgroup.write_text(
        "2:cpu,cpuacct:/system.slice/massflow.service\n", encoding="utf-8"
    )
    proc_mountinfo = tmp_path / "mountinfo"
    proc_mountinfo.write_text(
        f"31 25 0:28 / {cgroup_root} rw - cgroup cgroup rw,cpu,cpuacct\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(numba_runtime, "_PROC_SELF_CGROUP", proc_cgroup)
    monkeypatch.setattr(numba_runtime, "_PROC_SELF_MOUNTINFO", proc_mountinfo)

    assert numba_runtime._detect_linux_cgroup_cpu_count() == 2


def test_detect_linux_cgroup_cpu_count_uses_fixed_root_fallback(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    cgroup_root = tmp_path / "cgroup"
    cgroup_root.mkdir()
    (cgroup_root / "cpu.max").write_text("150000 100000\n", encoding="utf-8")

    monkeypatch.setattr(numba_runtime, "_PROC_SELF_CGROUP", tmp_path / "missing")
    monkeypatch.setattr(numba_runtime, "_PROC_SELF_MOUNTINFO", tmp_path / "missing")
    monkeypatch.setattr(numba_runtime, "_CGROUP_V2_FALLBACK_ROOT", cgroup_root)
    monkeypatch.setattr(numba_runtime, "_CGROUP_V1_FALLBACK_ROOTS", ())

    assert numba_runtime._detect_linux_cgroup_cpu_count() == 2


def test_detect_performance_core_workers_returns_positive_int() -> None:
    count = detect_performance_core_workers()
    assert isinstance(count, int)
    assert count >= 1


def test_detect_performance_core_workers_windows_branch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(numba_runtime.platform, "system", lambda: "Windows")
    monkeypatch.setattr(
        numba_runtime, "_detect_windows_intel_performance_workers", lambda: 6
    )
    monkeypatch.setattr(numba_runtime, "get_logical_cpu_count", lambda: 8)
    monkeypatch.setattr(numba_runtime.config, "NUMBA_NUM_THREADS", 8)

    assert detect_performance_core_workers() == 6


def test_detect_performance_core_workers_fallback_minimum_one(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(numba_runtime.platform, "system", lambda: "Unknown")
    monkeypatch.setattr(numba_runtime, "get_logical_cpu_count", lambda: 1)

    assert detect_performance_core_workers() == 1


def test_detect_performance_core_workers_fallback_uses_available_workers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(numba_runtime.platform, "system", lambda: "Unknown")
    monkeypatch.setattr(numba_runtime, "get_logical_cpu_count", lambda: 8)
    monkeypatch.setattr(numba_runtime.config, "NUMBA_NUM_THREADS", 8)

    assert detect_performance_core_workers() == 8


def test_detect_performance_core_workers_caps_numba_thread_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(numba_runtime.platform, "system", lambda: "Unknown")
    monkeypatch.setattr(numba_runtime, "get_logical_cpu_count", lambda: 64)
    monkeypatch.setattr(numba_runtime.config, "NUMBA_NUM_THREADS", 2)

    assert detect_performance_core_workers() == 2
