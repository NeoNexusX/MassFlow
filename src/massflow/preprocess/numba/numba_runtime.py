from __future__ import annotations

import ctypes
import os
import platform
import struct
import subprocess
from pathlib import Path, PurePosixPath
from threading import Lock
from typing import Iterator, Optional

from numba import config, get_num_threads, set_num_threads

from massflow.tools.logger import get_logger


logger = get_logger(__name__)
_NUMBA_RUNTIME_STATE: dict[str, Optional[int]] = {"max_workers": None}
_AUTO_NUMBA_WORKERS: Optional[int] = None
_NUMBA_RUNTIME_LOCK = Lock()
_PROC_SELF_CGROUP = Path("/proc/self/cgroup")
_PROC_SELF_MOUNTINFO = Path("/proc/self/mountinfo")
_CGROUP_V2_FALLBACK_ROOT = Path("/sys/fs/cgroup")
_CGROUP_V1_FALLBACK_ROOTS = (
    Path("/sys/fs/cgroup/cpu"),
    Path("/sys/fs/cgroup/cpu,cpuacct"),
)


# ---------------------------------------------------------------------------
# Available CPU capacity
# ---------------------------------------------------------------------------


def get_logical_cpu_count() -> int:
    """Return the logical CPU capacity available to the current process."""
    limits: list[int] = []

    if (system_count := os.cpu_count()) is not None and system_count > 0:
        limits.append(system_count)

    try:
        # Some platform stubs omit this optional API; unsupported calls are caught.
        affinity_count = len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
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


def _read_text_file(file_path: Path) -> Optional[str]:
    try:
        return file_path.read_text(encoding="utf-8", errors="ignore").strip()
    except OSError:
        return None


# --- Linux cgroup quotas ----------------------------------------------------


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
        limits: list[int] = []
        for filesystem, current, mount_point in locations:
            read_count = (
                _read_v2_cpu_count if filesystem == "cgroup2" else _read_v1_cpu_count
            )
            for cgroup_dir in _iter_cgroup_ancestors(current, mount_point):
                count = read_count(cgroup_dir)
                if count is not None:
                    limits.append(count)
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


# ---------------------------------------------------------------------------
# Performance-core topology
# ---------------------------------------------------------------------------

# --- Apple ------------------------------------------------------------------


def _read_sysctl_int(key: str) -> Optional[int]:
    try:
        output = subprocess.check_output(
            ["sysctl", "-n", key], text=True, stderr=subprocess.DEVNULL
        )
        return _to_positive_int(output)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _detect_apple_performance_workers() -> Optional[int]:
    for key in ("hw.perflevel0.logicalcpu", "hw.perflevel0.physicalcpu"):
        if (count := _read_sysctl_int(key)) is not None:
            return count
    return None


# --- Linux Intel ------------------------------------------------------------


def _detect_linux_intel_performance_workers() -> Optional[int]:
    cpu_list = _read_text_file(Path("/sys/devices/cpu_core/cpus"))
    if not cpu_list:
        return None

    performance_cpus: set[int] = set()
    try:
        for entry in cpu_list.split(","):
            first, separator, last = entry.partition("-")
            start = int(first)
            end = int(last) if separator else start
            if start < 0 or end < start:
                return None
            performance_cpus.update(range(start, end + 1))
    except ValueError:
        return None

    try:
        performance_cpus.intersection_update(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except (AttributeError, OSError):
        pass

    return len(performance_cpus) or None


# --- Windows ---------------------------------------------------------------


def _detect_windows_intel_performance_workers() -> Optional[int]:
    # RelationProcessorCore records contain exactly one GROUP_AFFINITY.
    # Layout: Relationship, Size, Flags, EfficiencyClass, Reserved[20],
    # GroupCount, GroupMask.Mask (pointer-sized), GroupMask.Group, Reserved[3].
    record = struct.Struct("@IIBB20xHPH6x")
    # This backend runs only on Windows; POSIX stubs omit its ctypes APIs.
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)  # type: ignore[attr-defined]
    query = kernel32.GetLogicalProcessorInformationEx
    query.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint32)]
    query.restype = ctypes.c_int
    relation_processor_core = 0
    error_insufficient_buffer = 122

    size = ctypes.c_uint32(0)
    query(relation_processor_core, None, ctypes.byref(size))
    if ctypes.get_last_error() != error_insufficient_buffer:  # type: ignore[attr-defined]
        return None

    buffer = ctypes.create_string_buffer(size.value)
    if not query(relation_processor_core, buffer, ctypes.byref(size)):
        return None

    workers_by_class: dict[int, int] = {}
    offset = 0
    while offset < size.value:
        if size.value - offset < record.size:
            return None
        relationship, length, _, efficiency, group_count, mask, _ = record.unpack_from(
            buffer, offset
        )
        if (
            relationship != relation_processor_core
            or group_count != 1
            or length < record.size
            or offset + length > size.value
        ):
            return None
        workers_by_class[efficiency] = (
            workers_by_class.get(efficiency, 0) + mask.bit_count()
        )
        offset += length

    # Higher EfficiencyClass values identify higher-performance cores.
    return workers_by_class[max(workers_by_class)] if workers_by_class else None


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
    # Numba injects NUMBA_NUM_THREADS into its config module at runtime.
    return max(
        1,
        min(preferred_workers, available_workers, config.NUMBA_NUM_THREADS),  # type: ignore[attr-defined]  # pylint: disable=no-member
    )


# ---------------------------------------------------------------------------
# Runtime configuration and thread control
# ---------------------------------------------------------------------------


def _resolve_workers(target: int) -> int:
    if not isinstance(target, int):
        raise TypeError("max_workers must be an integer or None")
    if target < 1:
        raise ValueError("max_workers must be >= 1")

    limit = min(get_logical_cpu_count(), config.NUMBA_NUM_THREADS)  # type: ignore[attr-defined]  # pylint: disable=no-member
    if target > limit:
        logger.warning(
            "Requested max_workers=%d exceeds the available limit=%d; using %d workers.",
            target,
            limit,
            limit,
        )
        return limit
    return target


def _set_threads_if_needed(target: int) -> None:
    if get_num_threads() != target:
        set_num_threads(target)


def set_global_numba_runtime(
    max_workers: Optional[int] = None,
) -> dict[str, Optional[int]]:
    """Set shared worker options and apply an explicit value to this thread.

    Other Python threads must call apply_numba_runtime() to use these options.
    Passing None clears the override and cached automatic count; this thread's
    mask is unchanged until its next apply_numba_runtime() call. Positive values
    above the available CPU or Numba pool limit are capped with a warning;
    the effective value is stored in the shared options.
    """
    global _AUTO_NUMBA_WORKERS  # pylint: disable=global-statement

    if max_workers is not None:
        max_workers = _resolve_workers(max_workers)
        _set_threads_if_needed(max_workers)

    with _NUMBA_RUNTIME_LOCK:
        _NUMBA_RUNTIME_STATE["max_workers"] = max_workers
        if max_workers is None:
            _AUTO_NUMBA_WORKERS = None
        return dict(_NUMBA_RUNTIME_STATE)


def get_global_numba_runtime() -> dict[str, Optional[int]]:
    """Get shared worker options, not the calling thread's active thread count."""
    with _NUMBA_RUNTIME_LOCK:
        return dict(_NUMBA_RUNTIME_STATE)


def apply_numba_runtime(override_workers: Optional[int] = None) -> int:
    """Apply a worker count to the calling thread and return the effective count.

    Precedence is per-call override, shared override, then cached automatic
    detection. Positive counts are capped by available CPUs and Numba's pool
    limit, with a warning when capped. Automatic detection is cached until
    set_global_numba_runtime(None).
    The resulting thread mask remains active after this call.
    """
    global _AUTO_NUMBA_WORKERS  # pylint: disable=global-statement

    if override_workers is not None:
        target = override_workers
    else:
        with _NUMBA_RUNTIME_LOCK:
            target = _NUMBA_RUNTIME_STATE["max_workers"]
            if target is None:
                if _AUTO_NUMBA_WORKERS is None:
                    _AUTO_NUMBA_WORKERS = detect_performance_core_workers()
                target = _AUTO_NUMBA_WORKERS

    target = _resolve_workers(target)
    _set_threads_if_needed(target)
    return target
