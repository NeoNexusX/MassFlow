from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    MSDataManager,
    MSDataManagerImzML,
)
from massflow.module import MassSpectrumSet
from massflow.tools.infer_spectrum_type import SpectrumType, set_spectrum_type_metadata

from massflow.preprocess.pipeline.models import FlatChunk, OutputFormat, PipelineConfig, PipelineResult


class OutputSink(Protocol):
    """Write target used by the flat pipeline writer stage."""

    def write_chunk(self, chunk: FlatChunk) -> None:
        ...

    def finalize(self) -> PipelineResult:
        ...

    def close(self) -> None:
        ...


@dataclass(slots=True)
class ImzMLSink:
    data_manager: MSDataManagerImzML

    def write_chunk(self, chunk: FlatChunk) -> None:
        self.data_manager.write_flat_batch(
            intensity_flat=chunk.intensity,
            lengths=chunk.lengths,
            coordinates=chunk.coordinates,
            mz_axis=chunk.mz_data,
        )

    def finalize(self) -> MSDataManagerImzML:
        self.data_manager.close_writer()
        self.data_manager.load_head_data()
        return self.data_manager

    def close(self) -> None:
        self.data_manager.close_writer()
        self.data_manager.close()


@dataclass(slots=True)
class IonImageSink:
    data_manager: IonImageDataManagerBinary | IonImageDataManagerH5

    def write_chunk(self, chunk: FlatChunk) -> None:
        if chunk.mz_data is None:
            raise ValueError("ion image output requires mz_data from flat_generator.")
        self.data_manager.write_flat_batch(
            intensity_flat=chunk.intensity,
            lengths=chunk.lengths,
            coordinates=chunk.coordinates,
            mz_axis=chunk.mz_data,
        )

    def finalize(self) -> IonImageDataManagerBinary | IonImageDataManagerH5:
        self.data_manager.finalize()
        return self.data_manager

    def close(self) -> None:
        self.data_manager.close()


class OutputSinkFactory:
    """Create output sinks and export existing managers to ion image formats."""

    _ION_IMAGE: dict[str, tuple[type[IonImageDataManagerBinary] | type[IonImageDataManagerH5], str]] = {
        "ion_image_binary": (IonImageDataManagerBinary, ".iib"),
        "ion_image_h5": (IonImageDataManagerH5, ".h5"),
    }

    def __init__(self, config: PipelineConfig):
        self.config = config

    def make_sink(
        self,
        data_manager: MSDataManager,
        output_spectrum_type: SpectrumType,
        final_output_format: OutputFormat | None = None,
    ) -> OutputSink:
        if final_output_format in self._ION_IMAGE:
            return IonImageSink(self.make_ion_image_manager(data_manager, final_output_format))
        return self.make_imzml_sink(data_manager, output_spectrum_type)

    def make_imzml_sink(self, data_manager: MSDataManager, output_spectrum_type: SpectrumType) -> ImzMLSink:
        processed_ms = MassSpectrumSet()
        processed_data_manager = MSDataManagerImzML(processed_ms, temp_dir=self.config.temp_dir)
        processed_data_manager.copy_meta(data_manager)
        set_spectrum_type_metadata(processed_data_manager.ms.meta, output_spectrum_type)
        return ImzMLSink(processed_data_manager)

    def make_ion_image_manager(
        self,
        data_manager: MSDataManager,
        output_format: OutputFormat,
    ) -> IonImageDataManagerBinary | IonImageDataManagerH5:
        cls, suffix = self._ION_IMAGE[str(output_format)]
        return cls.from_data_manager(
            data_manager,
            self.default_ion_image_output_path(suffix),
            temp_dir=self.config.temp_dir,
            temp_workers=max(1, self.config.output_max_workers),
            merge_workers=max(1, self.config.output_max_workers),
        )

    def default_ion_image_output_path(self, suffix: str) -> str:
        if self.config.output_path is not None:
            return self.config.output_path
        base_dir = self.config.temp_dir if self.config.temp_dir is not None else tempfile.gettempdir()
        return str(Path(base_dir) / f"ion_image_output{suffix}")

    def export_ion_image(self, data_manager: MSDataManager) -> IonImageDataManagerBinary | IonImageDataManagerH5:
        output_dm = self.make_ion_image_manager(data_manager, self.config.output_format)
        try:
            for mz_axis, intensity_flat, lengths, coordinates in data_manager.flat_generator(
                batch_size=self.config.batch_size,
                include_mz=True,
                max_threads=max(2, self.config.output_max_workers),
            ):
                output_dm.write_flat_batch(
                    intensity_flat=intensity_flat,
                    lengths=lengths,
                    coordinates=coordinates,
                    mz_axis=mz_axis,
                )
            output_dm.finalize()
            return output_dm
        except BaseException:
            output_dm.close()
            raise


def close_if_intermediate(data_manager: Any, original_data_manager: MSDataManager) -> None:
    if data_manager is not original_data_manager:
        data_manager.close()
