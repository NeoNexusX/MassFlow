import numpy as np
from numpy.typing import NDArray

from massflow.data_manager import MSDataManager, MSDataManagerImzML, IonImageDataManagerBinary, IonImageDataManagerH5
from massflow.preprocess.pipeline.flat_executor import run_flat_task
from massflow.preprocess.pipeline.dataset_executor import run_dataset_task
from massflow.preprocess.pipeline.planner import Segment
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.pipeline.runner")


def _del_last_dm(
    original_dm: MSDataManager,
    del_dm: MSDataManager,
) -> None:
    """Helper function to delete intermediate data managers after a segment is done."""
    if del_dm is not original_dm:
        del_dm.close()
        del del_dm

def _flat_segment(
    data_manager: MSDataManager,
    segment: Segment,
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5:
    """Run a batch-level segment and return the resulting data manager."""
    return run_flat_task(
        data_manager=data_manager, # type: ignore
        tasks=segment.tasks,
        input_meta=segment.input_meta,
        output_meta=segment.output_meta,
    )

def _dataset_segment(
    data_manager: MSDataManager,
    segment: Segment,
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5 | tuple[NDArray[np.float64], float]:
    """Run a dataset-level segment and return the resulting data manager or reference data."""
    return run_dataset_task(
        data_manager=data_manager, # type: ignore
        task=segment.tasks[0],
        input_meta=segment.input_meta,
        output_meta=segment.output_meta,
    )

def run_segments(
    data_manager: MSDataManager,
    task_segments: list[Segment],
):
    """Run the planned task segments and return the final data manager."""
    current_dm = data_manager
    for i, segment in enumerate(task_segments):
        logger.info(f"Running segment with scope '{segment.scope}' and tasks {[task.name for task in segment.tasks]}")
        if segment.scope == "batch":
            next_dm = _flat_segment(
                data_manager=current_dm, # type: ignore
                segment=segment,
            )
            _del_last_dm(data_manager, current_dm) # type: ignore
            current_dm = next_dm

        elif segment.scope == "dataset":
            result = _dataset_segment(
                data_manager=current_dm, # type: ignore
                segment=segment,
            )
            if isinstance(result, tuple):
                reference, tolerance = result
                kwargs = dict(task_segments[i+1].tasks[0].kwargs)
                if kwargs.get("reference") is None:
                    kwargs["reference"] = reference
                if kwargs.get("tolerance") is None:
                    kwargs["tolerance"] = tolerance
                task_segments[i + 1].tasks[0].kwargs = kwargs
            else:
                next_dm = result
                _del_last_dm(data_manager, current_dm) # type: ignore
                current_dm = next_dm

    return current_dm
