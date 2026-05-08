from __future__ import annotations

from typing import Callable, cast

from massflow.data_manager import MSDataManager, MSDataManagerImzML
from massflow.tools.infer_spectrum_type import SpectrumType, set_spectrum_type_metadata

from massflow.preprocess.pipeline.models import PipelineResult, PreprocessTask


class DatasetTaskExecutor:
    """Run dataset-scoped preprocessing tasks."""

    def run(
        self,
        *,
        data_manager: MSDataManager,
        task: PreprocessTask,
        output_spectrum_type: SpectrumType,
    ) -> PipelineResult:
        if task.name == "peak_align" or task.name == "peak_pick":
            dataset_apply_fn = cast(Callable[..., MSDataManagerImzML], task.apply_fn)
            result = dataset_apply_fn(data_manager=data_manager, **task.kwargs)
            set_spectrum_type_metadata(result.ms.meta, output_spectrum_type)
            return result

        raise NotImplementedError(f"Dataset task {task.name} is not supported.")
