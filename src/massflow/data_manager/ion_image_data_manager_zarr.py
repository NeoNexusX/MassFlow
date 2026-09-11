"""Ion-major MassFlow MSI Zarr data manager."""

from __future__ import annotations

import os
import uuid
from copy import deepcopy
from pathlib import Path

import numpy as np

from massflow.data_manager.ion_image_data_manager import IonImageDataManager
from massflow.msi_zarr import (
    MSIZarrReader,
    MSIZarrWriter,
    ZARR_MODE_CONTINUOUS,
    update_metadata_attrs,
)
from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.tools.infer_spectrum_type import SpectrumType, set_spectrum_type_metadata


class IonImageDataManagerZarr(IonImageDataManager):
    """Read and write continuous ion-major MSI Zarr stores."""

    def __init__(
        self,
        filepath: str | os.PathLike[str] | None = None,
        *,
        mode: str = "r",
        width: int | None = None,
        height: int | None = None,
        temp_dir: str | os.PathLike[str] | None = None,
        mz_dtype: np.dtype | type = np.float64,
        intensity_dtype: np.dtype | type = np.float32,
        ion_block: int = 256,
        temp_workers: int = 4,
        merge_workers: int = 4,
        chunk_size: int = 1,
        metadata=None,
        shared_mz_axis: np.ndarray | None = None,
        compute_mean_spectrum: bool = False,
        compute_rms: bool = False,
        include_spectra: bool = True,
    ) -> None:
        del merge_workers
        if mode not in {"r", "w"}:
            raise ValueError("mode must be 'r' or 'w'.")
        if filepath is None:
            if mode == "r":
                raise ValueError("filepath is required in read mode.")
            filepath = Path("output") / f"{uuid.uuid4().hex}.zarr"
        path = Path(filepath)
        self.filepath = path if path.suffix == ".zarr" else Path(f"{path}.zarr")
        self.mode = mode
        self._width = width
        self._height = height
        self.mz_dtype = np.dtype(mz_dtype).newbyteorder("<")
        if self.mz_dtype.kind != "f" or self.mz_dtype.itemsize not in {4, 8}:
            raise ValueError("mz_dtype must be float32 or float64.")
        self.intensity_dtype = np.dtype(intensity_dtype)
        self.ion_block = max(1, int(ion_block))
        self.meta = MetaDataZarr()
        if metadata is not None:
            self.meta.copy(metadata, exclude=False)
        self._reader: MSIZarrReader | None = None
        self._writer: MSIZarrWriter | None = None
        self._finalized = False
        self._shared_mz_axis = (
            None
            if shared_mz_axis is None
            else np.asarray(shared_mz_axis, dtype=self.mz_dtype)
        )
        if mode == "r":
            self._open_for_read()
        else:
            self._validate_metadata(metadata, self._shared_mz_axis)
            self._writer = MSIZarrWriter(
                self.filepath,
                row_axis="ion",
                encoding=ZARR_MODE_CONTINUOUS,
                metadata=self.meta,
                width=width,
                height=height,
                mz_dtype=self.mz_dtype,
                intensity_dtype=self.intensity_dtype,
                target_chunk_bytes=max(1, int(chunk_size)) * 1024 * 1024,
                temp_dir=temp_dir,
                temp_workers=temp_workers,
                ion_block=self.ion_block,
                compute_mean_spectrum=compute_mean_spectrum,
                compute_rms=compute_rms,
                include_spectra=include_spectra,
            )

    def _require_reader(self) -> MSIZarrReader:
        if self._reader is None:
            self._open_for_read()
        assert self._reader is not None
        return self._reader

    def _require_writer(self) -> MSIZarrWriter:
        if self._writer is None:
            raise ValueError("Writer is only available in write mode.")
        return self._writer

    @property
    def width(self) -> int:
        if self._width is None:
            self._open_for_read()
        assert self._width is not None
        return int(self._width)

    @property
    def height(self) -> int:
        if self._height is None:
            self._open_for_read()
        assert self._height is not None
        return int(self._height)

    @property
    def n_ions(self) -> int:
        return self._require_reader().n_ions

    @property
    def n_pixels(self) -> int:
        return self._require_reader().n_pixels

    @property
    def mz_axis(self) -> np.ndarray:
        return self._require_reader().read_shared_mz_axis()

    @property
    def coordinates(self) -> np.ndarray:
        return self._require_reader().read_coordinates()

    @property
    def mean_spectrum(self) -> np.ndarray | None:
        return self._require_reader().read_mean_spectrum()

    @property
    def rms(self) -> np.ndarray | None:
        return self._require_reader().read_rms()

    def write_flat_batch(
        self,
        intensity_flat: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
        mz_axis: np.ndarray,
    ) -> None:
        self._require_writer().write_flat_batch(
            intensity_flat=intensity_flat,
            lengths=lengths,
            coordinates=coordinates,
            mz_axis=mz_axis,
        )

    def add_metadata(self, metadata: dict) -> None:
        """Upsert key/value pairs into the on-disk metadata attrs."""
        update_metadata_attrs(self.filepath, metadata)

    def finalize(self) -> None:
        if self.mode != "w" or self._finalized:
            return
        self._require_writer().finalize()
        self._finalized = True
        self._open_for_read()

    def get_ion_image(self, ion_idx: int) -> np.ndarray:
        return self._require_reader().read_ion_image(ion_idx)

    def get_ion_vector(self, ion_idx: int) -> np.ndarray:
        return self._require_reader().read_ion_vector(ion_idx)

    def get_ion_block(self, start: int, end: int) -> np.ndarray:
        return self._require_reader().read_rows(start, end)

    def get_pixel_spectrum(self, pixel_idx: int) -> np.ndarray:
        raise NotImplementedError(
            "Pixel spectrum access is not efficient for ion-major Zarr stores."
        )

    def close(self) -> None:
        if self._writer is not None:
            self._writer.close()
        if self._reader is not None:
            self._reader.close()
        self._writer = None
        self._reader = None

    def _open_for_read(self) -> None:
        if self._reader is not None:
            return
        self._reader = MSIZarrReader(self.filepath)
        if self._reader.row_axis != "ion":
            raise ValueError("IonImageDataManagerZarr requires an ion-major store.")
        self._height, self._width = self._reader.spatial_shape
        self.intensity_dtype = self._reader.intensity_dtype
        self.mz_dtype = self._reader.mz_dtype
        self.meta = self._reader.metadata

    @staticmethod
    def _validate_metadata(
        metadata,
        shared_mz_axis: np.ndarray | None,
    ) -> None:
        if metadata is None:
            return
        if getattr(metadata, "continuous", None) is not True:
            raise ValueError("Ion-major Zarr output requires continuous=True.")
        if getattr(metadata, "centroid_spectrum", None) is False:
            raise ValueError("Ion-major Zarr output requires centroid data.")
        if shared_mz_axis is None or shared_mz_axis.size == 0:
            raise ValueError("Ion-major Zarr output requires a shared m/z axis.")

    @classmethod
    def from_data_manager(
        cls,
        data_manager,
        filepath,
        *,
        temp_dir=None,
        mz_dtype=np.float64,
        intensity_dtype=None,
        ion_block=256,
        temp_workers=2,
        merge_workers=2,
        chunk_size=1,
        compute_mean_spectrum=False,
        compute_rms=False,
        include_spectra=True,
    ) -> "IonImageDataManagerZarr":
        meta = data_manager.ms.meta
        if (
            getattr(meta, "continuous", None) is not True
            or data_manager.ms.shared_mz_list is None
            or len(data_manager.ms.shared_mz_list) == 0
        ):
            raise ValueError(
                "from_data_manager requires continuous=True and a shared m/z axis."
            )
        return cls(
            filepath=filepath,
            mode="w",
            width=int(meta.max_count_of_pixels_x),
            height=int(meta.max_count_of_pixels_y),
            temp_dir=temp_dir,
            mz_dtype=mz_dtype,
            intensity_dtype=(
                data_manager.intensity_dtype
                if intensity_dtype is None
                else intensity_dtype
            ),
            ion_block=ion_block,
            temp_workers=temp_workers,
            merge_workers=merge_workers,
            chunk_size=chunk_size,
            metadata=meta,
            shared_mz_axis=data_manager.ms.shared_mz_list,
            compute_mean_spectrum=compute_mean_spectrum,
            compute_rms=compute_rms,
            include_spectra=include_spectra,
        )

    @classmethod
    def from_flat_output_contract(
        cls,
        data_manager,
        filepath,
        *,
        shared_mz_axis: np.ndarray,
        output_spectrum_type: SpectrumType = "centroid",
        temp_dir=None,
        mz_dtype=np.float64,
        intensity_dtype=None,
        ion_block=256,
        temp_workers=2,
        merge_workers=2,
        chunk_size=1,
        compute_mean_spectrum=False,
        compute_rms=False,
        include_spectra=True,
    ) -> "IonImageDataManagerZarr":
        output_meta = deepcopy(data_manager.ms.meta)
        output_meta.continuous = True
        output_meta.processed = None
        set_spectrum_type_metadata(output_meta, output_spectrum_type)
        return cls(
            filepath=filepath,
            mode="w",
            width=int(output_meta.max_count_of_pixels_x),
            height=int(output_meta.max_count_of_pixels_y),
            temp_dir=temp_dir,
            mz_dtype=mz_dtype,
            intensity_dtype=(
                data_manager.intensity_dtype
                if intensity_dtype is None
                else intensity_dtype
            ),
            ion_block=ion_block,
            temp_workers=temp_workers,
            merge_workers=merge_workers,
            chunk_size=chunk_size,
            metadata=output_meta,
            shared_mz_axis=np.asarray(shared_mz_axis, dtype="<f8"),
            compute_mean_spectrum=compute_mean_spectrum,
            compute_rms=compute_rms,
            include_spectra=include_spectra,
        )
