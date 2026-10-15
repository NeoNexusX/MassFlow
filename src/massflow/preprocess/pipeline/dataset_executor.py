from __future__ import annotations
from typing import Any, cast, Callable

import numpy as np
from numpy.typing import NDArray

from massflow.data_manager import (
    MSDataManager,
    MSDataManagerImzML,
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
)
from massflow.preprocess.pipeline.types import (
    ION_IMAGE_OUTPUT_FORMATS,
    OUTPUT_FORMAT_ION_IMAGE_BINARY,
    OUTPUT_FORMAT_ION_IMAGE_H5,
    OUTPUT_FORMAT_ION_IMAGE_ZARR,
    InputMetaData,
    OutputFormat,
    OutputMetaData,
    PreprocessTask,
)
from massflow.tools.infer_spectrum_type import set_spectrum_type_metadata
from massflow.tools.logger import get_logger
logger = get_logger("massflow.preprocess.pipeline.dataset_executor")

IonImageDataManagerType = type[IonImageDataManagerBinary] | type[IonImageDataManagerH5] | type[IonImageDataManagerZarr]

_ION_IMAGE_CLASSES: dict[OutputFormat, IonImageDataManagerType] = {
    OUTPUT_FORMAT_ION_IMAGE_BINARY: IonImageDataManagerBinary,
    OUTPUT_FORMAT_ION_IMAGE_H5: IonImageDataManagerH5,
    OUTPUT_FORMAT_ION_IMAGE_ZARR: IonImageDataManagerZarr,
}


def ion_image_writer(
    output_format: OutputFormat,
    input_meta: InputMetaData,
) -> tuple[IonImageDataManagerType, dict[str, Any]]:
    """Return the ion-image manager class for ``output_format`` and its writer options."""
    ion_image_cls = _ION_IMAGE_CLASSES.get(output_format)
    if ion_image_cls is None:
        raise ValueError(f"Unsupported ion image output format: {output_format}")
    options: dict[str, Any] = {
        "temp_dir": input_meta.temp_dir,
        "temp_workers": input_meta.output_max_workers,
        "merge_workers": input_meta.output_max_workers,
        "ion_block": input_meta.ion_block,
        "compute_mean_spectrum": input_meta.compute_mean_spectrum,
    }
    if output_format == OUTPUT_FORMAT_ION_IMAGE_ZARR:
        options["chunk_size"] = input_meta.zarr_chunk_size
        options["compute_rms"] = input_meta.compute_rms
    return ion_image_cls, options

def _write_through_generator(
        data_manager: MSDataManager,
        output_dm: IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr,
        input_meta: InputMetaData,
) -> IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr:
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
) -> IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr:
    """Convert an imzML manager to an ion image data manager by writing through a generator."""
    ion_image_cls, options = ion_image_writer(output_format, input_meta)
    output_dm = ion_image_cls.from_data_manager(
        data_manager=imzml_manager,
        filepath=file_path,
        **options,
    )

    return _write_through_generator(imzml_manager, output_dm, input_meta)

def _task_to_data_manager(
    data_manager: MSDataManager,
    task: PreprocessTask,
    input_meta: InputMetaData,
    output_meta: OutputMetaData,
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr:
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
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr | tuple[NDArray[np.float64], float]:
    """Execute whole-dataset task."""

    if task.name == "reference_computer":
        return _task_to_array(data_manager, task)

    return _task_to_data_manager(data_manager, task, input_meta, output_meta)
