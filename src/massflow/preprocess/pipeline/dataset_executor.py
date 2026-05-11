from __future__ import annotations
from typing import cast, Callable

import numpy as np
from numpy.typing import NDArray

from massflow.data_manager import MSDataManager, MSDataManagerImzML, IonImageDataManagerBinary, IonImageDataManagerH5
from massflow.preprocess.pipeline.types import (
    ION_IMAGE_OUTPUT_FORMATS,
    OUTPUT_FORMAT_ION_IMAGE_BINARY,
    OUTPUT_FORMAT_ION_IMAGE_H5,
    InputMetaData,
    OutputFormat,
    OutputMetaData,
    PreprocessTask,
)
from massflow.tools.infer_spectrum_type import set_spectrum_type_metadata
from massflow.tools.logger import get_logger
logger = get_logger("massflow.preprocess.pipeline.dataset_executor")

def _write_through_generator(
        data_manager: MSDataManager,
        output_dm: IonImageDataManagerBinary | IonImageDataManagerH5,
        input_meta: InputMetaData,
) -> IonImageDataManagerBinary | IonImageDataManagerH5:
    """Use a generator to write data from an imzML manager to an ion image manager."""
    for mz_axis, intensity_flat, lengths, coordinates in data_manager.flat_generator(
        batch_size=input_meta.batch_size,
        include_mz=True,
        max_threads=input_meta.output_max_workers,
    ):
        output_dm.write_flat_batch(
            intensity_flat=intensity_flat,
            lengths=lengths,
            coordinates=coordinates,
            mz_axis=mz_axis,
        )
    output_dm.finalize()
    return output_dm

def imzml_to_ion_image_data_manager(
    imzml_manager: MSDataManagerImzML,
    input_meta: InputMetaData,
    output_format: OutputFormat,
    file_path: str,
) -> IonImageDataManagerBinary | IonImageDataManagerH5:
    """Convert an imzML manager to an ion image data manager by writing through a generator."""
    if output_format == OUTPUT_FORMAT_ION_IMAGE_BINARY:
        ion_image_cls = IonImageDataManagerBinary
    elif output_format == OUTPUT_FORMAT_ION_IMAGE_H5:
        ion_image_cls = IonImageDataManagerH5
    else:
        raise ValueError(f"Unsupported output format: {output_format}")

    output_dm = ion_image_cls.from_data_manager(
        data_manager=imzml_manager,
        filepath=file_path,
        temp_dir=input_meta.temp_dir,
        temp_workers=input_meta.output_max_workers,
        merge_workers=input_meta.output_max_workers,
    )

    return _write_through_generator(imzml_manager, output_dm, input_meta)

def _task_to_data_manager(
    data_manager: MSDataManager,
    task: PreprocessTask,
    input_meta: InputMetaData,
    output_meta: OutputMetaData,
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5:
    """Apply a task to a data manager and return an imzML manager."""
    apply_fn = cast(Callable[..., MSDataManagerImzML], task.apply_fn)
    processed_data = apply_fn(data_manager=data_manager, **task.kwargs)
    set_spectrum_type_metadata(processed_data.ms.meta, output_meta.spectrum_type)

    if output_meta.output_format in ION_IMAGE_OUTPUT_FORMATS:
        ion_dm = imzml_to_ion_image_data_manager(
            imzml_manager=processed_data,
            input_meta=input_meta,
            output_format=output_meta.output_format,
            file_path=output_meta.file_path, # type: ignore
        )
        processed_data.close()
        return ion_dm

    return processed_data

def _task_to_array(
    data_manager: MSDataManager,
    task: PreprocessTask,
) -> tuple[NDArray[np.float64], float]:
    """Apply a task to a data manager and return an array."""
    apply_fn = cast(Callable[..., tuple[NDArray[np.float64], float]], task.apply_fn)
    return apply_fn(data_manager=data_manager, **task.kwargs)

def run_dataset_task(
    data_manager: MSDataManager,
    task: PreprocessTask,
    input_meta: InputMetaData,
    output_meta: OutputMetaData,
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5 | tuple[NDArray[np.float64], float]:
    """Execute whole-dataset task."""

    if task.name == "reference_computer":
        return _task_to_array(data_manager, task)

    return _task_to_data_manager(data_manager, task, input_meta, output_meta)
