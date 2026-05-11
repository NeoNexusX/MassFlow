from typing import Optional
from dataclasses import dataclass
from pathlib import Path
import tempfile
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from massflow.preprocess.flat_pre_fun import FlatBatchResult
from massflow.preprocess.api import TaskScope
from massflow.tools.infer_spectrum_type import SpectrumType
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.pipeline.planner")


_OPERATION_ORDER: dict[str, int] = {
    "baseline_correction": 10,
    "noise_reduction": 20,
    "normalization": 30,
    "peak_pick": 40,
    "reference_computer": 50,
    "peak_align": 60,
    }
_PROFILE_ONLY_TASKS = frozenset({"baseline_correction", "noise_reduction", "peak_pick"})

@dataclass(slots=True)
class InputMetaData:
    """Structured metadata for output datasets."""
    batch_size: int
    temp_dir: str
    numba_max_threads: Optional[int]
    output_max_workers: int
    queue_ab_size: int
    queue_bc_size: int
    result_path: str
    result_format: str

@dataclass(slots=True)
class PreprocessTask:
    """Lazy task registration model for batch/dataset preprocessing."""
    name: str
    apply_fn: Callable[..., FlatBatchResult] | Callable[..., Any]
    scope: TaskScope = "batch"
    kwargs: dict[str, Any] = field(default_factory=dict)

def _normalize_output_path(input_meta: InputMetaData) -> str:
    """Helper function to ensure output paths are properly set based on input metadata."""
    input_meta.temp_dir = input_meta.temp_dir or tempfile.gettempdir()

    suffix = ".iib" if input_meta.result_format == "ion_image_binary" else ".h5"

    if input_meta.result_path:
        path = Path(input_meta.result_path)

        if path.suffix != suffix:
            path = path.with_suffix(suffix)
    else:
        path = Path(input_meta.temp_dir) / f"ion_image_output{suffix}"

    return str(path)

@dataclass(slots=True)
class OutputMetaData:
    """Structured metadata for output datasets."""
    output_format: str
    spectrum_type: SpectrumType
    file_path: Optional[str]

@dataclass(slots=True)
class Segment:
    """Represents a segment of tasks to be executed together, either at the batch level or dataset level."""

    scope: str  # "batch" or "dataset"
    tasks: list[PreprocessTask]
    input_meta: InputMetaData
    output_meta: OutputMetaData

def _resolve_input_meta(input_meta: InputMetaData) -> InputMetaData:
    """Helper function to resolve and validate input metadata."""
    input_meta.result_path = _normalize_output_path(input_meta)
    input_meta.numba_max_threads = input_meta.numba_max_threads or None
    input_meta.output_max_workers = max(1, input_meta.output_max_workers)
    return input_meta

def _arrange_tasks(tasks, spectrum_type) -> list[PreprocessTask]:

    def sort_key(task: PreprocessTask) -> int:
        return _OPERATION_ORDER.get(task.name, 10_000)

    ordered_tasks = sorted(tasks, key=sort_key)

    if spectrum_type == "centroid":
        dropped_task_names = [task.name for task in ordered_tasks if task.name in _PROFILE_ONLY_TASKS]
        if dropped_task_names:
            logger.warning(f"Detected centroid spectra; dropping unsupported tasks: {dropped_task_names}")
            ordered_tasks = [task for task in ordered_tasks if task.name not in _PROFILE_ONLY_TASKS]

    if ordered_tasks:
        task_names = [f"{task.name}:{task.scope}" for task in ordered_tasks]
        logger.info(f"async_pipeline_task_order: {task_names}")
    else:
        logger.warning("No supported preprocessing tasks remain after task arrangement.")

    return ordered_tasks

def make_task_segments(tasks, input_spectrum_type, input_meta: InputMetaData) -> tuple[list[Segment], InputMetaData]:
    """Execute registered tasks and return new data manager with processed data."""
    ordered_tasks = _arrange_tasks(tasks, input_spectrum_type)
    input_meta = _resolve_input_meta(input_meta)
    task_segments = []
    current_spectrum_type = input_spectrum_type
    task_cursor = 0

    while task_cursor < len(ordered_tasks):
        batch_tasks: list[PreprocessTask] = []

        while task_cursor < len(ordered_tasks) and ordered_tasks[task_cursor].scope == "batch":
            batch_tasks.append(ordered_tasks[task_cursor])
            current_spectrum_type = "centroid" if ordered_tasks[task_cursor].name == "peak_pick" else current_spectrum_type
            task_cursor += 1
        if batch_tasks:
            output_meta = OutputMetaData(output_format="imzml", spectrum_type=current_spectrum_type, file_path=None)
            task_segments.append(Segment(scope="batch", tasks=batch_tasks, input_meta=input_meta, output_meta=output_meta))

        if task_cursor >= len(ordered_tasks):
            break

        dataset_task = ordered_tasks[task_cursor]
        task_cursor += 1
        output_meta = OutputMetaData(output_format="imzml", spectrum_type=current_spectrum_type, file_path=None)
        task_segments.append(Segment(scope="dataset", tasks=[dataset_task], input_meta=input_meta, output_meta=output_meta))

    if task_segments:
        task_segments[-1].output_meta.output_format = input_meta.result_format
        task_segments[-1].output_meta.file_path = input_meta.result_path

    return task_segments, input_meta
