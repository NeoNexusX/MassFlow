from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

from massflow.tools.infer_spectrum_type import SpectrumType

TaskScope = Literal["batch", "dataset"]
SegmentScope = Literal["batch", "dataset"]
OutputFormat = Literal["imzml", "ion_image_binary", "ion_image_h5", "ion_image_zarr"]
IonImageOutputFormat = Literal["ion_image_binary", "ion_image_h5", "ion_image_zarr"]

OUTPUT_FORMAT_IMZML: OutputFormat = "imzml"
OUTPUT_FORMAT_ION_IMAGE_BINARY: IonImageOutputFormat = "ion_image_binary"
OUTPUT_FORMAT_ION_IMAGE_H5: IonImageOutputFormat = "ion_image_h5"
OUTPUT_FORMAT_ION_IMAGE_ZARR: IonImageOutputFormat = "ion_image_zarr"
ION_IMAGE_OUTPUT_FORMATS: frozenset[IonImageOutputFormat] = frozenset(
    {OUTPUT_FORMAT_ION_IMAGE_BINARY, OUTPUT_FORMAT_ION_IMAGE_H5, OUTPUT_FORMAT_ION_IMAGE_ZARR}
)
OUTPUT_FORMATS: frozenset[OutputFormat] = frozenset(
    {OUTPUT_FORMAT_IMZML, OUTPUT_FORMAT_ION_IMAGE_BINARY, OUTPUT_FORMAT_ION_IMAGE_H5, OUTPUT_FORMAT_ION_IMAGE_ZARR}
)
OUTPUT_FORMAT_SUFFIXES: dict[IonImageOutputFormat, str] = {
    OUTPUT_FORMAT_ION_IMAGE_BINARY: ".iib",
    OUTPUT_FORMAT_ION_IMAGE_H5: ".h5",
    OUTPUT_FORMAT_ION_IMAGE_ZARR: ".zarr",
}


@dataclass(frozen=True, slots=True)
class InputMetaData:
    """Structured runtime settings for pipeline execution."""

    batch_size: int
    temp_dir: str | None
    numba_max_threads: int | None
    output_max_workers: int
    queue_ab_size: int
    queue_bc_size: int
    result_path: str | None
    result_format: OutputFormat
    ion_block: int = 32


@dataclass(slots=True)
class PreprocessTask:
    """Lazy task registration model for batch/dataset preprocessing."""

    name: str
    apply_fn: Callable[..., Any]
    scope: TaskScope = "batch"
    kwargs: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class OutputMetaData:
    """Structured metadata for output datasets."""

    output_format: OutputFormat
    spectrum_type: SpectrumType
    file_path: str | None


@dataclass(frozen=True, slots=True)
class Segment:
    """Represents a planned pipeline segment."""

    scope: SegmentScope
    tasks: list[PreprocessTask]
    input_meta: InputMetaData
    output_meta: OutputMetaData
