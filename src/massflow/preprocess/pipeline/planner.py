from dataclasses import replace
from pathlib import Path
import tempfile

from massflow.preprocess.pipeline.types import (
    ION_IMAGE_OUTPUT_FORMATS,
    OUTPUT_FORMAT_IMZML,
    OUTPUT_FORMAT_SUFFIXES,
    InputMetaData,
    OutputMetaData,
    PreprocessTask,
    Segment,
)
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

def _normalize_output_path(input_meta: InputMetaData, temp_dir: str) -> str | None:
    """Helper function to ensure output paths are properly set based on input metadata."""
    if input_meta.result_format not in ION_IMAGE_OUTPUT_FORMATS:
        return input_meta.result_path

    suffix = OUTPUT_FORMAT_SUFFIXES[input_meta.result_format]
    if input_meta.result_path:
        path = Path(input_meta.result_path)

        if path.suffix != suffix:
            path = path.with_suffix(suffix)
    else:
        path = Path(temp_dir) / f"ion_image_output{suffix}"

    return str(path)

def _resolve_input_meta(input_meta: InputMetaData) -> InputMetaData:
    """Helper function to resolve and validate input metadata."""
    temp_dir = input_meta.temp_dir or tempfile.gettempdir()
    return replace(
        input_meta,
        temp_dir=temp_dir,
        result_path=_normalize_output_path(input_meta, temp_dir),
        numba_max_threads=input_meta.numba_max_threads or None,
        output_max_workers=max(1, input_meta.output_max_workers),
    )

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
            output_meta = OutputMetaData(output_format=OUTPUT_FORMAT_IMZML, spectrum_type=current_spectrum_type, file_path=None)
            task_segments.append(Segment(scope="batch", tasks=batch_tasks, input_meta=input_meta, output_meta=output_meta))

        if task_cursor >= len(ordered_tasks):
            break

        dataset_task = ordered_tasks[task_cursor]
        task_cursor += 1
        if dataset_task.name == "peak_pick":
            current_spectrum_type = "centroid"
        output_meta = OutputMetaData(output_format=OUTPUT_FORMAT_IMZML, spectrum_type=current_spectrum_type, file_path=None)
        task_segments.append(Segment(scope="dataset", tasks=[dataset_task], input_meta=input_meta, output_meta=output_meta))

    if task_segments:
        last_segment = task_segments[-1]
        task_segments[-1] = replace(
            last_segment,
            output_meta=replace(
                last_segment.output_meta,
                output_format=input_meta.result_format,
                file_path=input_meta.result_path,
            ),
        )

    return task_segments, input_meta
