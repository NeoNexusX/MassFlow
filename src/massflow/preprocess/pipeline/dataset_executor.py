from __future__ import annotations

from typing import Callable, cast

from massflow.data_manager import MSDataManager, MSDataManagerImzML
from massflow.preprocess.flat_pre_fun import FlatPreprocess
from massflow.preprocess.helper.peak_align_helper import reference_computer
from massflow.tools.infer_spectrum_type import SpectrumType, set_spectrum_type_metadata

from massflow.preprocess.pipeline.flat_executor import FlatPipelineExecutor
from massflow.preprocess.pipeline.models import PipelineConfig, PipelineResult, PreprocessTask
from massflow.preprocess.pipeline.sinks import OutputSinkFactory


class DatasetTaskExecutor:
    """Run dataset-scoped preprocessing tasks."""

    def __init__(
        self,
        config: PipelineConfig,
        flat_executor: FlatPipelineExecutor,
        sink_factory: OutputSinkFactory,
    ):
        self.config = config
        self.flat_executor = flat_executor
        self.sink_factory = sink_factory

    def prepare_flat_peak_align_task(self, *, data_manager: MSDataManager, task: PreprocessTask) -> PreprocessTask:
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
            kwargs=kwargs,
        )

    def run(
        self,
        *,
        data_manager: MSDataManager,
        task: PreprocessTask,
        output_spectrum_type: SpectrumType,
        final_output_format: str | None = None,
    ) -> PipelineResult:
        if task.name == "peak_align" and task.apply_fn is FlatPreprocess.peak_align_flat:
            prepared_task = self.prepare_flat_peak_align_task(data_manager=data_manager, task=task)
            sink = self.sink_factory.make_sink(
                data_manager=data_manager,
                output_spectrum_type=output_spectrum_type,
                final_output_format=final_output_format,
            )
            return self.flat_executor.run(data_manager=data_manager, tasks=[prepared_task], sink=sink)

        if task.name == "peak_align" or task.name == "peak_pick":
            dataset_apply_fn = cast(Callable[..., MSDataManagerImzML], task.apply_fn)
            result = dataset_apply_fn(data_manager=data_manager, **task.kwargs)
            set_spectrum_type_metadata(result.ms.meta, output_spectrum_type)
            return result

        raise NotImplementedError(f"Dataset task {task.name} is not supported.")
