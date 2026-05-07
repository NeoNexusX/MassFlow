from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal, Optional

import numpy as np

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    MSDataManagerImzML,
)
from massflow.preprocess.api import TaskScope
from massflow.preprocess.flat_pre_fun import FlatBatchResult
from massflow.tools.infer_spectrum_type import SpectrumType

OutputFormat = Literal["imzml", "ion_image_binary", "ion_image_h5"]
ExecutionKind = Literal["batch", "dataset"]
FlatTaskCallable = Callable[..., FlatBatchResult]
DatasetTaskCallable = Callable[..., MSDataManagerImzML]
TaskCallable = FlatTaskCallable | DatasetTaskCallable
PipelineResult = MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5

ION_IMAGE_FORMATS: frozenset[str] = frozenset({"ion_image_binary", "ion_image_h5"})
SUPPORTED_OUTPUT_FORMATS: frozenset[str] = frozenset({"imzml", *ION_IMAGE_FORMATS})
PROFILE_ONLY_TASKS: frozenset[str] = frozenset({"baseline_correction", "noise_reduction", "peak_pick"})


@dataclass(slots=True)
class FlatChunk:
    """Transport object between flat pipeline stages."""

    mz_data: Optional[np.ndarray]
    intensity: np.ndarray
    lengths: np.ndarray
    coordinates: np.ndarray


@dataclass(slots=True)
class PreprocessTask:
    """Lazy task registration model for batch and dataset preprocessing."""

    name: str
    apply_fn: TaskCallable
    scope: TaskScope = "batch"
    kwargs: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class PipelineConfig:
    """Runtime configuration shared by pipeline components."""

    batch_size: int = 256
    queue_ab_size: int = 1
    queue_bc_size: int = 1
    temp_dir: str | None = "./temp"
    numba_max_threads: int | None = None
    output_format: OutputFormat = "imzml"
    output_path: str | None = None
    output_max_workers: int = 2

    def validate(self) -> None:
        if self.batch_size <= 0 or self.batch_size > 9056:
            raise ValueError("batch_size must be in the range [1, 9056].")
        if self.queue_ab_size <= 0 or self.queue_bc_size <= 0:
            raise ValueError("queue sizes must be positive.")
        if self.output_format not in SUPPORTED_OUTPUT_FORMATS:
            raise ValueError("output_format must be 'imzml', 'ion_image_binary', or 'ion_image_h5'.")
        if self.output_max_workers <= 0:
            raise ValueError("output_max_workers must be positive.")


@dataclass(slots=True)
class ExecutionSegment:
    """A contiguous executable section of the preprocessing plan."""

    kind: ExecutionKind
    tasks: list[PreprocessTask]
    output_spectrum_type: SpectrumType
    final_output_format: OutputFormat | None = None


@dataclass(slots=True)
class PipelinePlan:
    """Fully planned preprocessing execution."""

    input_spectrum_type: SpectrumType
    output_format: OutputFormat
    segments: list[ExecutionSegment]

    @property
    def wants_ion_image(self) -> bool:
        return self.output_format in ION_IMAGE_FORMATS
