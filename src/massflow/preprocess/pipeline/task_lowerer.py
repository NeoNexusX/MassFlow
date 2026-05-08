from __future__ import annotations

from massflow.data_manager import MSDataManager
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.preprocess.helper.peak_align_helper import reference_computer

from massflow.preprocess.pipeline.models import ExecutionSegment, PipelineConfig, PreprocessTask


class TaskLowerer:
    """Prepare dataset-aware tasks that can execute on the flat batch engine."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def lower_segment(self, *, data_manager: MSDataManager, segment: ExecutionSegment) -> ExecutionSegment:
        if segment.kind != "batch":
            return segment

        prepared_tasks: list[PreprocessTask] = []
        changed = False
        for task in segment.tasks:
            if task.requires_dataset_prepare:
                prepared_tasks.append(self.prepare_flat_peak_align_task(data_manager=data_manager, task=task))
                changed = True
            else:
                prepared_tasks.append(task)

        if not changed:
            return segment

        return ExecutionSegment(
            kind=segment.kind,
            tasks=prepared_tasks,
            output_spectrum_type=segment.output_spectrum_type,
            final_output_format=segment.final_output_format,
        )

    def prepare_flat_peak_align_task(self, *, data_manager: MSDataManager, task: PreprocessTask) -> PreprocessTask:
        if task.name != "peak_align" or task.apply_fn is not FlatPreprocess.peak_align_flat:
            raise NotImplementedError(f"Dataset preparation is not supported for task {task.name!r}.")

        kwargs = dict(task.kwargs)

        units = kwargs.get("units", "ppm")
        reference = kwargs.get("reference")
        tolerance = kwargs.get("tolerance")

        if reference is None or tolerance is None:
            reference, tolerance = reference_computer(
                data_manager=data_manager,
                reference=reference,
                tolerance=tolerance,
                units=units,
                binfun=kwargs.get("binfun", "median"),
                binratio=kwargs.get("binratio", 2.0),
                min_count=kwargs.get("min_count"),
                min_frequency=kwargs.get("min_frequency"),
                batch_size=self.config.batch_size,
            )
        else:
            tolerance = tolerance * 1e-6 if units == "ppm" else tolerance

        kwargs["reference"] = reference
        kwargs["tolerance"] = tolerance

        kwargs.pop("binfun", None)
        kwargs.pop("binratio", None)
        kwargs.pop("min_count", None)
        kwargs.pop("min_frequency", None)

        return PreprocessTask(
            name=task.name,
            scope="batch",
            apply_fn=task.apply_fn,
            requires_dataset_prepare=False,
            kwargs=kwargs,
        )
