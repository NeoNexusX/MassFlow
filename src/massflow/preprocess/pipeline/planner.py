from __future__ import annotations

from massflow.tools.infer_spectrum_type import SpectrumType
from massflow.tools.logger import get_logger

from massflow.preprocess.pipeline.models import (
    ExecutionSegment,
    ION_IMAGE_FORMATS,
    PROFILE_ONLY_TASKS,
    OutputFormat,
    PipelinePlan,
    PreprocessTask,
)

logger = get_logger("massflow.preprocess.pipeline.planner")


class TaskPlanner:
    """Build executable preprocessing segments from registered tasks."""

    def __init__(
        self,
        tasks: list[PreprocessTask],
        *,
        output_format: OutputFormat = "imzml",
    ):
        self.tasks = list(tasks)
        self.output_format = output_format

    @staticmethod
    def resolve_output_spectrum_type(
        *,
        input_spectrum_type: SpectrumType,
        tasks: list[PreprocessTask],
    ) -> SpectrumType:
        return "centroid" if any(task.name == "peak_pick" for task in tasks) else input_spectrum_type

    def arrange_tasks(self, spectrum_type: SpectrumType) -> list[PreprocessTask]:
        tasks = list(self.tasks)

        if spectrum_type == "centroid":
            dropped_task_names = [task.name for task in tasks if task.name in PROFILE_ONLY_TASKS]
            if dropped_task_names:
                logger.warning(f"Detected centroid spectra; dropping unsupported tasks: {dropped_task_names}")
                tasks = [task for task in tasks if task.name not in PROFILE_ONLY_TASKS]

        if tasks:
            task_names = [
                f"{task.name}:{task.scope}{':prepare' if task.requires_dataset_prepare else ''}"
                for task in tasks
            ]
            logger.info(f"preprocess_task_order: {task_names}")
        else:
            logger.warning("No supported preprocessing tasks remain after task arrangement.")

        return tasks

    def plan(self, input_spectrum_type: SpectrumType) -> PipelinePlan:
        tasks = self.arrange_tasks(input_spectrum_type)
        segments: list[ExecutionSegment] = []
        current_spectrum_type = input_spectrum_type
        task_cursor = 0
        wants_ion_image = self.output_format in ION_IMAGE_FORMATS

        while task_cursor < len(tasks):
            batch_tasks: list[PreprocessTask] = []

            while (
                task_cursor < len(tasks)
                and tasks[task_cursor].scope == "batch"
                and not tasks[task_cursor].requires_dataset_prepare
            ):
                batch_tasks.append(tasks[task_cursor])
                task_cursor += 1

            if batch_tasks:
                output_spectrum_type = self.resolve_output_spectrum_type(
                    input_spectrum_type=current_spectrum_type,
                    tasks=batch_tasks,
                )
                final_output_format = self.output_format if wants_ion_image and task_cursor >= len(tasks) else None
                segments.append(
                    ExecutionSegment(
                        kind="batch",
                        tasks=batch_tasks,
                        output_spectrum_type=output_spectrum_type,
                        final_output_format=final_output_format,
                    )
                )
                current_spectrum_type = output_spectrum_type

            if task_cursor >= len(tasks):
                break

            prepared_batch_task = tasks[task_cursor]
            if prepared_batch_task.scope == "batch" and prepared_batch_task.requires_dataset_prepare:
                task_cursor += 1
                output_spectrum_type = self.resolve_output_spectrum_type(
                    input_spectrum_type=current_spectrum_type,
                    tasks=[prepared_batch_task],
                )
                final_output_format = self.output_format if wants_ion_image and task_cursor >= len(tasks) else None
                segments.append(
                    ExecutionSegment(
                        kind="batch",
                        tasks=[prepared_batch_task],
                        output_spectrum_type=output_spectrum_type,
                        final_output_format=final_output_format,
                    )
                )
                current_spectrum_type = output_spectrum_type
                continue

            dataset_task = tasks[task_cursor]
            task_cursor += 1
            output_spectrum_type = self.resolve_output_spectrum_type(
                input_spectrum_type=current_spectrum_type,
                tasks=[dataset_task],
            )
            final_output_format = self.output_format if wants_ion_image and task_cursor >= len(tasks) else None
            segments.append(
                ExecutionSegment(
                    kind="dataset",
                    tasks=[dataset_task],
                    output_spectrum_type=output_spectrum_type,
                    final_output_format=final_output_format,
                )
            )
            current_spectrum_type = output_spectrum_type

        return PipelinePlan(
            input_spectrum_type=input_spectrum_type,
            output_format=self.output_format,
            segments=segments,
        )
