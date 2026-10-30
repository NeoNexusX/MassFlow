from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import tempfile
import uuid

from massflow.preprocess.pipeline.types import (
    ION_IMAGE_OUTPUT_FORMATS,
    OUTPUT_FORMAT_IMZML,
    OUTPUT_FORMAT_SPECTRA_ZARR,
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
    "filter_peaks_by_isotopes": 50,
    "reference_computer": 60,
    "peak_align": 70,
    "filter_peaks_by_adducts": 80,
    "filter_peaks_by_stats": 90,
    }
_PROFILE_ONLY_TASKS = frozenset({"baseline_correction", "noise_reduction", "peak_pick"})

def _normalize_output_path(input_meta: InputMetaData, temp_dir: str) -> str | None:
    """Helper function to ensure output paths are properly set based on input metadata."""
    suffix = OUTPUT_FORMAT_SUFFIXES[input_meta.result_format]
    if input_meta.result_path:
        path = Path(input_meta.result_path)

        if path.suffix != suffix:
            path = path.with_suffix(suffix)
    else:
        output_name = "ion_image_output" if input_meta.result_format in ION_IMAGE_OUTPUT_FORMATS else "preprocess_output"
        path = Path(temp_dir) / uuid.uuid4().hex / f"{output_name}{suffix}"

    return str(path)

def _resolve_input_meta(input_meta: InputMetaData) -> InputMetaData:
    """Helper function to resolve and validate input metadata."""
    temp_dir = input_meta.temp_dir or tempfile.gettempdir()
    return replace(
        input_meta,
        temp_dir=temp_dir,
        result_path=_normalize_output_path(input_meta, temp_dir),
    )

def _arrange_tasks(tasks, spectrum_type) -> list[PreprocessTask]:

    def sort_key(task: PreprocessTask) -> int:
        return _OPERATION_ORDER.get(task.name, 10_000)

    ordered_tasks = sorted(tasks, key=sort_key)

    if spectrum_type == "centroid":
        dropped_task_names = [task.name for task in ordered_tasks if task.name in _PROFILE_ONLY_TASKS]
        if dropped_task_names:
            logger.warning(f"Detected centroid spectra; dropping unsupported tasks: {dropped_task_names}") # pylint: disable=W1203
            ordered_tasks = [task for task in ordered_tasks if task.name not in _PROFILE_ONLY_TASKS]

    if ordered_tasks:
        task_names = [f"{task.name}:{task.scope}" for task in ordered_tasks]
        logger.info(f"async_pipeline_task_order: {task_names}") # pylint: disable=W1203
    else:
        logger.warning("No supported preprocessing tasks remain after task arrangement.")

    return ordered_tasks

def make_task_segments(tasks, input_spectrum_type, input_meta: InputMetaData) -> tuple[list[Segment], InputMetaData]:
    """Execute registered tasks and return new data manager with processed data."""
    # 规划工作在任务副本上：后续 planner 注入 output_path、runner 回写
    # reference/tolerance 都只落在副本，不污染用户注册的原始任务，
    # 保证对同一个 Preprocessor 重复调用 start() 时行为一致。
    ordered_tasks = _arrange_tasks([deepcopy(task) for task in tasks], input_spectrum_type)
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

    if input_meta.result_format == OUTPUT_FORMAT_SPECTRA_ZARR and (
        not task_segments or task_segments[-1].scope == "dataset"
    ):
        task_segments.append(
            Segment(
                scope="batch",
                tasks=[],
                input_meta=input_meta,
                output_meta=OutputMetaData(
                    output_format=OUTPUT_FORMAT_SPECTRA_ZARR,
                    spectrum_type=current_spectrum_type,
                    file_path=input_meta.result_path,
                ),
            )
        )
    elif task_segments:
        last_segment = task_segments[-1]
        task_segments[-1] = replace(
            last_segment,
            output_meta=replace(
                last_segment.output_meta,
                output_format=input_meta.result_format,
                file_path=input_meta.result_path,
            ),
        )
        # dataset-scope tasks build their own output data manager internally, so
        # the executor cannot route them through output_meta.file_path the way
        # batch tasks do. Inject result_path as the task's `output_path` kwarg
        # so the task writes directly to the user-requested location.
        if last_segment.scope == "dataset" and input_meta.result_format == OUTPUT_FORMAT_IMZML:
            task_segments[-1].tasks[0].kwargs["output_path"] = input_meta.result_path

    return task_segments, input_meta
