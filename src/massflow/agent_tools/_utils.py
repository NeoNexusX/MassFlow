from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence
import json
import math
import shutil
import tempfile

import numpy as np

from .schemas import ToolResult


def make_json_safe(value: Any) -> Any:
    """Convert numpy-heavy nested objects into JSON-safe Python values."""
    if isinstance(value, dict):
        return {str(key): make_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [make_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        return None
    return value


def write_result_json(
    result: ToolResult,
    output_dir: str | None,
    file_name: str,
) -> str | None:
    """Persist a tool result as JSON and return the artifact path."""
    if output_dir is None:
        return None

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    summary_path = output_path / file_name
    summary_path.write_text(
        json.dumps(make_json_safe(result.to_dict()), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return str(summary_path)


def ensure_directory(output_dir: str | None, *, temp_prefix: str = "massflow_agent") -> Path:
    """Return a writable directory, creating a temp directory when absent."""
    if output_dir is None:
        return Path(tempfile.mkdtemp(prefix=f"{temp_prefix}_"))

    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def copy_artifact(source_path: str | None, destination_path: Path) -> str | None:
    """Copy a file artifact when it exists and return the copied path."""
    if source_path is None:
        return None

    source = Path(source_path)
    if not source.exists():
        return None

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination_path)
    return str(destination_path)


def build_even_sample_indices(total: int, sample_size: int) -> list[int]:
    """Pick evenly distributed sample indices without duplicates."""
    if total <= 0 or sample_size <= 0:
        return []

    if sample_size >= total:
        return list(range(total))

    raw = np.linspace(0, total - 1, num=sample_size, dtype=int)
    unique = sorted(set(int(index) for index in raw.tolist()))

    if len(unique) == sample_size:
        return unique

    seen = set(unique)
    for candidate in range(total):
        if candidate not in seen:
            unique.append(candidate)
            seen.add(candidate)
        if len(unique) == sample_size:
            break

    return sorted(unique)


def summarize_numeric(values: Sequence[float]) -> dict[str, float | int | None]:
    """Return a compact summary for a numeric sequence."""
    if not values:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "std": None,
        }

    array = np.asarray(values, dtype=np.float64)
    return {
        "count": int(array.size),
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "mean": float(np.mean(array)),
        "median": float(np.median(array)),
        "std": float(np.std(array)),
    }


def summarize_coordinates(
    coordinates: Iterable[tuple[int, int, int]],
) -> dict[str, int | float | None]:
    """Summarize a coordinate stream for agent consumption."""
    coords = list(coordinates)
    if not coords:
        return {
            "spectra_count": 0,
            "min_x": None,
            "max_x": None,
            "min_y": None,
            "max_y": None,
            "min_z": None,
            "max_z": None,
            "unique_z_count": 0,
        }

    xs = np.asarray([coord[0] for coord in coords], dtype=np.int32)
    ys = np.asarray([coord[1] for coord in coords], dtype=np.int32)
    zs = np.asarray([coord[2] for coord in coords], dtype=np.int32)
    return {
        "spectra_count": int(len(coords)),
        "min_x": int(np.min(xs)),
        "max_x": int(np.max(xs)),
        "min_y": int(np.min(ys)),
        "max_y": int(np.max(ys)),
        "min_z": int(np.min(zs)),
        "max_z": int(np.max(zs)),
        "unique_z_count": int(len(set(int(z) for z in zs.tolist()))),
    }
