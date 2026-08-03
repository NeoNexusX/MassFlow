from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path
from typing import Optional

from numba import get_num_threads, set_num_threads


_NUMBA_RUNTIME_STATE: dict[str, Optional[int]] = {"max_workers": None}

# Core detection logic part
def get_logical_cpu_count() -> int:
    """Return the logical CPU count visible to the current process."""
    return int(count) if (count := os.cpu_count()) and count > 0 else 1


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
    sys_name = platform.system()
    detected = None

    if sys_name == "Darwin":
        detected = _detect_apple_performance_workers()
    elif sys_name == "Linux":
        detected = _detect_linux_intel_performance_workers()
    elif sys_name == "Windows":
        detected = _detect_windows_intel_performance_workers()

    return detected if detected is not None else max(2, get_logical_cpu_count() // 2)


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
