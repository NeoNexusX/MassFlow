from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path, PurePosixPath
from typing import Iterator, Optional

from numba import config, get_num_threads, set_num_threads


_NUMBA_RUNTIME_STATE: dict[str, Optional[int]] = {"max_workers": None}
_PROC_SELF_CGROUP = Path("/proc/self/cgroup")
_PROC_SELF_MOUNTINFO = Path("/proc/self/mountinfo")
_CGROUP_V2_FALLBACK_ROOT = Path("/sys/fs/cgroup")
_CGROUP_V1_FALLBACK_ROOTS = (
    Path("/sys/fs/cgroup/cpu"),
    Path("/sys/fs/cgroup/cpu,cpuacct"),
)


# Core detection logic part
def get_logical_cpu_count() -> int:
    """Return the logical CPU capacity available to the current process."""
    limits: list[int] = []

    if (system_count := os.cpu_count()) is not None and system_count > 0:
        limits.append(system_count)

    try:
        affinity_count = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        pass
    else:
        if affinity_count > 0:
            limits.append(affinity_count)

    if platform.system() == "Linux":
        if (cgroup_count := _detect_linux_cgroup_cpu_count()) is not None:
            limits.append(cgroup_count)

    return min(limits) if limits else 1


def _to_positive_int(raw: str) -> Optional[int]:
    try:
        return val if (val := int(raw.strip())) > 0 else None
    except (TypeError, ValueError):
        return None


def _read_sysctl_int(key: str) -> Optional[int]:
    try:
        output = subprocess.check_output(
            ["sysctl", "-n", key], text=True, stderr=subprocess.DEVNULL
        )
        return _to_positive_int(output)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _read_int_file(file_path: Path) -> Optional[int]:
    try:
        return _to_positive_int(file_path.read_text(encoding="utf-8", errors="ignore"))
    except OSError:
        return None


def _read_text_file(file_path: Path) -> Optional[str]:
    try:
        return file_path.read_text(encoding="utf-8", errors="ignore").strip()
    except OSError:
        return None


def _cpu_count_from_quota(quota_raw: str, period_raw: str) -> Optional[int]:
    if quota_raw == "max":
        return None

    quota = _to_positive_int(quota_raw)
    period = _to_positive_int(period_raw)
    if quota is None or period is None:
        return None

    return (quota + period - 1) // period


def _unescape_mountinfo_path(raw: str) -> str:
    for escaped, value in (
        (r"\040", " "),
        (r"\011", "\t"),
        (r"\012", "\n"),
        (r"\134", "\\"),
    ):
        raw = raw.replace(escaped, value)
    return raw


def _read_current_cgroup_paths() -> tuple[Optional[str], Optional[str]]:
    content = _read_text_file(_PROC_SELF_CGROUP)
    if content is None:
        return None, None

    v2_path = None
    v1_cpu_path = None
    for line in content.splitlines():
        fields = line.split(":", 2)
        if len(fields) != 3:
            continue

        hierarchy_id, controllers, cgroup_path = fields
        if hierarchy_id == "0" and not controllers and cgroup_path.startswith("/"):
            v2_path = cgroup_path
        elif "cpu" in controllers.split(",") and cgroup_path.startswith("/"):
            v1_cpu_path = cgroup_path

    return v2_path, v1_cpu_path


def _read_cpu_cgroup_mounts() -> list[tuple[str, Path, str]]:
    content = _read_text_file(_PROC_SELF_MOUNTINFO)
    if content is None:
        return []

    mounts: list[tuple[str, Path, str]] = []
    for line in content.splitlines():
        fields = line.split()
        try:
            separator = fields.index("-")
            filesystem = fields[separator + 1]
            super_options = fields[separator + 3].split(",")
            mount_root = _unescape_mountinfo_path(fields[3])
            mount_point = Path(_unescape_mountinfo_path(fields[4]))
        except (IndexError, ValueError):
            continue

        if filesystem == "cgroup2":
            mounts.append((filesystem, mount_point, mount_root))
        elif filesystem == "cgroup" and "cpu" in super_options:
            mounts.append((filesystem, mount_point, mount_root))

    return mounts


def _resolve_current_cgroup_dir(
    mount_point: Path, mount_root_raw: str, cgroup_path_raw: str
) -> Optional[Path]:
    mount_root = PurePosixPath(mount_root_raw)
    cgroup_path = PurePosixPath(cgroup_path_raw)

    try:
        relative_path = cgroup_path.relative_to(mount_root)
    except ValueError:
        return None

    if relative_path == PurePosixPath("."):
        return mount_point
    return mount_point.joinpath(*relative_path.parts)


def _iter_cgroup_ancestors(current: Path, mount_point: Path) -> Iterator[Path]:
    if current != mount_point and mount_point not in current.parents:
        return

    cgroup_dir = current
    while True:
        yield cgroup_dir
        if cgroup_dir == mount_point:
            return
        cgroup_dir = cgroup_dir.parent


def _read_v2_cpu_count(cgroup_dir: Path) -> Optional[int]:
    cpu_max = _read_text_file(cgroup_dir / "cpu.max")
    if cpu_max is None:
        return None

    fields = cpu_max.split()
    if len(fields) != 2:
        return None
    return _cpu_count_from_quota(fields[0], fields[1])


def _read_v1_cpu_count(cgroup_dir: Path) -> Optional[int]:
    quota = _read_text_file(cgroup_dir / "cpu.cfs_quota_us")
    period = _read_text_file(cgroup_dir / "cpu.cfs_period_us")
    if quota is None or period is None:
        return None
    return _cpu_count_from_quota(quota, period)


def _detect_linux_cgroup_cpu_count() -> Optional[int]:
    v2_path, v1_cpu_path = _read_current_cgroup_paths()
    locations: list[tuple[str, Path, Path]] = []
    for filesystem, mount_point, mount_root in _read_cpu_cgroup_mounts():
        cgroup_path = v2_path if filesystem == "cgroup2" else v1_cpu_path
        if cgroup_path is None:
            continue
        current = _resolve_current_cgroup_dir(mount_point, mount_root, cgroup_path)
        if current is not None:
            locations.append((filesystem, current, mount_point))

    if locations:
        limits = [
            count
            for filesystem, current, mount_point in locations
            for cgroup_dir in _iter_cgroup_ancestors(current, mount_point)
            if (
                count := (
                    _read_v2_cpu_count(cgroup_dir)
                    if filesystem == "cgroup2"
                    else _read_v1_cpu_count(cgroup_dir)
                )
            )
            is not None
        ]
        return min(limits) if limits else None

    fallback_limits = [
        count
        for count in (
            _read_v2_cpu_count(_CGROUP_V2_FALLBACK_ROOT),
            *(_read_v1_cpu_count(root) for root in _CGROUP_V1_FALLBACK_ROOTS),
        )
        if count is not None
    ]

    return min(fallback_limits) if fallback_limits else None


def _detect_apple_performance_workers() -> Optional[int]:
    for key in ("hw.perflevel0.logicalcpu", "hw.perflevel0.physicalcpu"):
        if (count := _read_sysctl_int(key)) is not None:
            return count
    return None


def _detect_linux_intel_performance_workers() -> Optional[int]:
    try:
        if "GenuineIntel" not in Path("/proc/cpuinfo").read_text(
            encoding="utf-8", errors="ignore"
        ):
            return None

    except OSError:
        return None

    core_types = [
        val
        for p in Path("/sys/devices/system/cpu").glob("cpu[0-9]*/topology/core_type")
        if (val := _read_int_file(p)) is not None
    ]

    if not core_types:
        return None

    perf_workers = core_types.count(max(core_types))
    return perf_workers if perf_workers > 0 else None


def _run_windows_shell_text(script: str) -> Optional[str]:
    commands = (
        [
            "powershell",
            "-NoProfile",
            "-NonInteractive",
            "-ExecutionPolicy",
            "Bypass",
            "-Command",
            script,
        ],
        ["pwsh", "-NoProfile", "-NonInteractive", "-Command", script],
    )

    for cmd in commands:
        try:
            output = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL)
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue

        text = output.strip()
        if text:
            return text

    return None


def _detect_windows_intel_performance_workers() -> Optional[int]:
    proc_id = os.environ.get("PROCESSOR_IDENTIFIER", "")
    if "intel" not in proc_id.lower():
        manufacturer = _run_windows_shell_text(
            "(Get-CimInstance Win32_Processor -ErrorAction SilentlyContinue | "
            "Select-Object -First 1 -ExpandProperty Manufacturer)"
        )
        if manufacturer is None or "intel" not in manufacturer.lower():
            return None

    detected = _run_windows_shell_text(
        "$sum = (Get-CimInstance Win32_Processor -ErrorAction SilentlyContinue | "
        "Measure-Object -Property NumberOfPerformanceCores -Sum).Sum; "
        "if ($null -ne $sum -and [int]$sum -gt 0) { [int]$sum }"
    )
    if detected is None:
        return None

    return _to_positive_int(detected)


def detect_performance_core_workers() -> int:
    """Detect preferred worker count from performance-core topology when available."""
    available_workers = get_logical_cpu_count()
    sys_name = platform.system()
    detected = None

    if sys_name == "Darwin":
        detected = _detect_apple_performance_workers()
    elif sys_name == "Linux":
        detected = _detect_linux_intel_performance_workers()
    elif sys_name == "Windows":
        detected = _detect_windows_intel_performance_workers()

    preferred_workers = detected if detected is not None else available_workers
    return max(
        1,
        min(preferred_workers, available_workers, config.NUMBA_NUM_THREADS),
    )


# Thread control part


def _set_threads_if_needed(target: int) -> None:
    if not isinstance(target, int):
        raise TypeError("max_workers must be an integer or None")
    if target < 1:
        raise ValueError("max_workers must be >= 1")
    if get_num_threads() != target:
        set_num_threads(target)


def set_global_numba_runtime(
    max_workers: Optional[int] = None,
) -> dict[str, Optional[int]]:
    """
    Set process-global Numba runtime options.
    Current scope manages only worker count. Passing None clears the override.
    """
    if max_workers is not None:
        _set_threads_if_needed(max_workers)

    _NUMBA_RUNTIME_STATE["max_workers"] = max_workers
    return get_global_numba_runtime()


def get_global_numba_runtime() -> dict[str, Optional[int]]:
    """Get process-global Numba runtime options."""
    return dict(_NUMBA_RUNTIME_STATE)


def apply_numba_runtime(override_workers: Optional[int] = None) -> int:
    """Apply per-call worker override, otherwise use global runtime setting.
    If the global setting is not initialized, it will optimize it once.
    """
    if override_workers is not None:
        target = override_workers
    else:
        target = _NUMBA_RUNTIME_STATE["max_workers"]
        if target is None:
            target = detect_performance_core_workers()

    _set_threads_if_needed(target)
    return target
