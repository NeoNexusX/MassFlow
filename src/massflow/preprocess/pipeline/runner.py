from __future__ import annotations

from typing import Any, cast

from massflow.data_manager import MSDataManager, MSDataManagerImzML

from massflow.preprocess.pipeline.dataset_executor import DatasetTaskExecutor
from massflow.preprocess.pipeline.flat_executor import FlatPipelineExecutor
from massflow.preprocess.pipeline.models import PipelineConfig, PipelinePlan, PipelineResult
from massflow.preprocess.pipeline.sinks import OutputSinkFactory, close_if_intermediate
from massflow.preprocess.pipeline.task_lowerer import TaskLowerer


class PipelineRunner:
    """Execute a planned preprocessing pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.sink_factory = OutputSinkFactory(config)
        self.flat_executor = FlatPipelineExecutor(config)
        self.dataset_executor = DatasetTaskExecutor()
        self.task_lowerer = TaskLowerer(config)

    def run(self, *, data_manager: MSDataManager, plan: PipelinePlan) -> PipelineResult | MSDataManager:
        current_data_manager: MSDataManager = data_manager

        if not plan.segments:
            if plan.wants_ion_image:
                return self.sink_factory.export_ion_image(current_data_manager)
            return current_data_manager

        for planned_segment in plan.segments:
            segment = self.task_lowerer.lower_segment(
                data_manager=current_data_manager,
                segment=planned_segment,
            )

            if segment.kind == "batch":
                sink = self.sink_factory.make_sink(
                    data_manager=current_data_manager,
                    output_spectrum_type=segment.output_spectrum_type,
                    final_output_format=segment.final_output_format,
                )
                next_data_manager = self.flat_executor.run(
                    data_manager=current_data_manager,
                    tasks=segment.tasks,
                    sink=sink,
                )
                close_if_intermediate(current_data_manager, data_manager)
                if segment.final_output_format is not None:
                    return next_data_manager
                current_data_manager = cast(MSDataManager, next_data_manager)
                continue

            if segment.kind == "dataset":
                task = segment.tasks[0]
                next_data_manager = self.dataset_executor.run(
                    data_manager=current_data_manager,
                    task=task,
                    output_spectrum_type=segment.output_spectrum_type,
                )
                close_if_intermediate(current_data_manager, data_manager)
                if segment.final_output_format is not None and not isinstance(next_data_manager, MSDataManagerImzML):
                    return next_data_manager
                current_data_manager = cast(MSDataManager, next_data_manager)
                continue

            raise ValueError(f"Unsupported execution segment kind: {segment.kind!r}")

        if plan.wants_ion_image:
            output_dm = self.sink_factory.export_ion_image(current_data_manager)
            close_if_intermediate(current_data_manager, data_manager)
            return output_dm

        return cast(MSDataManagerImzML, current_data_manager)


def close_result_if_needed(result: Any, original_data_manager: MSDataManager) -> None:
    close_if_intermediate(result, original_data_manager)
