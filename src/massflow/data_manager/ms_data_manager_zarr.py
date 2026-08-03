"""MS data manager for MassFlow spectra Zarr stores."""

from __future__ import annotations

from typing import Any, Generator

import numpy as np

from massflow.data_manager.ms_data_manager import MSDataManager
from massflow.msi_zarr import (
    MSIZarrReader,
    MSIZarrWriter,
    ZARR_MODE_AUTO,
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
    ZarrStorageMode,
    update_metadata_attrs,
)
from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.module.spectrum_zarr import SpectrumZarr
from massflow.tools.logger import get_logger

logger = get_logger("massflow.data_manager.ms_data_manager_zarr")


class MSDataManagerZarr(MSDataManager):
    """MassFlow spectra Zarr data manager with lazy spectrum loading."""

    def __init__(
        self,
        ms=None,
        target_locs=None,
        filepath=None,
        max_threads: int = 8,
        temp_dir=None,
        mz_dtype: Any = np.float64,
        intensity_dtype: Any = np.float32,
        mode: str = "r",
        storage_mode: ZarrStorageMode | None = None,
        chunk_size: int = 1,
        compute_tic: bool = False,
    ):
        if mode not in {"r", "w"}:
            raise ValueError("mode must be 'r' or 'w'")
        if storage_mode is not None and storage_mode not in {
            ZARR_MODE_CONTINUOUS,
            ZARR_MODE_PROCESSED,
        }:
            raise ValueError("storage_mode must be 'continuous' or 'processed'.")
        super().__init__(
            ms=ms,
            target_locs=target_locs,
            filepath=filepath,
            max_threads=max_threads,
            mz_dtype=mz_dtype,
            intensity_dtype=intensity_dtype,
            temp_dir=temp_dir,
        )
        self.mode = mode
        self._temp_dir = temp_dir
        self._chunk_size = max(1, int(chunk_size))
        self._compute_tic = bool(compute_tic)
        self.ms.meta = MetaDataZarr()
        self._storage_mode: ZarrStorageMode | None = storage_mode
        self._finalized = False
        self._reader: MSIZarrReader | None = None
        self._writer: MSIZarrWriter | None = None

    @property
    def zarr_filepath(self) -> str:
        """Return full .zarr path derived from the base path."""
        return f"{self.file_base_path}.zarr"

    @property
    def storage_mode(self) -> ZarrStorageMode:
        if self.mode == "r":
            self._open_arrays_for_read()
            assert self._storage_mode is not None
            return self._storage_mode
        encoding = self._require_writer().encoding
        if encoding == ZARR_MODE_AUTO:
            raise ValueError("Zarr storage mode has not been determined.")
        return encoding

    def _require_writer(self) -> MSIZarrWriter:
        """Return the internal streaming Zarr writer."""
        if self.ms.meta is None:
            raise ValueError(
                "MS meta data is None. Please load or copy meta data first."
            )
        if self._writer is None:
            width = getattr(self.ms.meta, "max_count_of_pixels_x", None)
            height = getattr(self.ms.meta, "max_count_of_pixels_y", None)
            if width is not None and int(width) <= 0:
                width = None
            if height is not None and int(height) <= 0:
                height = None
            self._writer = MSIZarrWriter(
                self.zarr_filepath,
                row_axis="pixel",
                encoding=(
                    ZARR_MODE_AUTO
                    if self._storage_mode is None
                    else self._storage_mode
                ),
                metadata=self.ms.meta, # type: ignore
                width=width,
                height=height,
                mz_dtype=self.mz_dtype,
                intensity_dtype=self.intensity_dtype,
                target_chunk_bytes=self._chunk_size * 1024 * 1024,
                temp_dir=self._temp_dir,
                compute_tic=self._compute_tic,
            )
        return self._writer

    def _require_reader(self) -> MSIZarrReader:
        if self._reader is None:
            self._reader = MSIZarrReader(self.zarr_filepath)
        return self._reader

    def copy_meta(self, source_dm: MSDataManager):
        """Copy metadata from another data manager, dropping fields that may change here."""
        self.ms.meta.copy(source_dm.ms.meta)

    def extract_metadata(self):
        """Load metadata from the Zarr metadata group."""
        self.ms.meta = self._require_reader().metadata

    @property
    def tic(self) -> np.ndarray | None:
        """Return the per-pixel TIC array from stats/tic, if present."""
        return self._require_reader().read_tic()

    def preload_hook(self, *args, **kwargs):
        """Load metadata and shared m/z axis before building spectrum placeholders."""
        self.extract_metadata()
        self._open_arrays_for_read()
        if self.storage_mode == ZARR_MODE_CONTINUOUS and self.ms.shared_mz_list is None:
            self.ms.shared_mz_list = self._require_reader().read_shared_mz_axis()

    def loading_hook(self, *args, **kwargs):
        """Hook invoked while spectra are being loaded."""

    def loaded_hook(self):
        """Finalize metadata after loading spectra."""
        if self.ms.meta.mask is None:
            self.create_ms_meta_mask()

    def load_head_data(self):
        """Load Zarr metadata and add lazy SpectrumZarr placeholders into `ms`."""
        logger.info(f"Loading data from Zarr file: {self.zarr_filepath}")
        self.preload_hook()
        reader = self._require_reader()
        coordinates = reader.read_coordinates()
        for index, coord in enumerate(coordinates):
            x, y = int(coord[0]), int(coord[1])
            if self.target_locs[0][0] <= x <= self.target_locs[1][0] and self.target_locs[0][1] <= y <= self.target_locs[1][1]:
                spectrum = SpectrumZarr(
                    reader=reader,
                    index=index,
                    coordinates=coord.tolist(),
                    shared_mz_list=self.ms.shared_mz_list,
                )
                self.ms.add_spectrum(spectrum)
                self.current_spectrum_num += 1
        self.loaded_hook()

    def write_flat_batch(
        self,
        intensity_flat: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
        mz_axis: np.ndarray | None,
    ) -> None:
        """Write a flat spectra batch through the Zarr writer."""
        if self.mode != "w":
            raise ValueError("write_flat_batch is only available in write mode.")
        if self._finalized:
            raise ValueError("Cannot write after finalize().")
        self._require_writer().write_flat_batch(
            intensity_flat=intensity_flat,
            lengths=lengths,
            coordinates=coordinates,
            mz_axis=mz_axis,
        )

    def finalize(self) -> None:
        """Finalize the Zarr writer."""
        if self.mode != "w" or self._finalized:
            return
        self._require_writer().finalize()
        self._finalized = True

    def add_metadata(self, metadata: dict) -> None:
        """Upsert key/value pairs into the on-disk metadata attrs."""
        update_metadata_attrs(self.zarr_filepath, metadata)

    def matrix_generator(
        self,
        batch_size: int = 256,
        include_mz: bool = True,
        max_threads: int = 0,
    ) -> Generator[tuple, Any, None]:
        """Yield matrix spectra batches from the Zarr arrays."""
        for mz_data, intensity_flat, lengths, coordinates in self.flat_generator(
            batch_size=batch_size,
            include_mz=include_mz,
            max_threads=max_threads,
        ):
            max_len = int(np.max(lengths))
            if self.storage_mode == ZARR_MODE_CONTINUOUS:
                intensity = np.asarray(intensity_flat).reshape((-1, max_len))
                yield mz_data, intensity, lengths, coordinates
                continue

            intensity = np.zeros((len(lengths), max_len), dtype=self.intensity_dtype)
            mz_matrix = (
                np.zeros((len(lengths), max_len), dtype=self.mz_dtype)
                if include_mz
                else None
            )
            starts = np.empty(len(lengths), dtype=np.int64)
            starts[0] = 0
            if len(lengths) > 1:
                np.cumsum(lengths[:-1], dtype=np.int64, out=starts[1:])
            for row, start in enumerate(starts):
                length = int(lengths[row])
                end = start + length
                intensity[row, :length] = intensity_flat[start:end]
                if mz_matrix is not None:
                    mz_matrix[row, :length] = mz_data[start:end]
            yield mz_matrix, intensity, lengths, coordinates

    def flat_generator(
        self,
        batch_size: int = 256,
        include_mz: bool = True,
        max_threads: int = 0,  # max_threads 用不到
    ) -> Generator[tuple, Any, None]:
        """Yield flat spectra batches from the Zarr arrays."""
        self._open_arrays_for_read()
        indices = self._target_indices()
        total = int(indices.size)
        if total == 0:
            return

        for batch_start in range(0, total, batch_size):
            batch_indices = indices[batch_start : batch_start + batch_size]
            mz_data, intensity, lengths, coordinates = (
                self._require_reader().read_flat_batch(
                    batch_indices,
                    include_mz=include_mz,
                )
            )
            if self.storage_mode == ZARR_MODE_CONTINUOUS:
                self.ms.shared_mz_list = mz_data
            yield mz_data, intensity, lengths, coordinates

    def close(self):
        """Close the Zarr reader and writer."""
        if self._reader is not None:
            self._reader.close()
            self._reader = None
        if self._writer is not None:
            self._writer.close()
            self._writer = None
        super().close()

    def _open_arrays_for_read(self) -> None:
        reader = self._require_reader()
        if reader.row_axis != "pixel":
            raise ValueError("MSDataManagerZarr requires a pixel-major Zarr store.")
        if reader.encoding not in {ZARR_MODE_CONTINUOUS, ZARR_MODE_PROCESSED}:
            raise ValueError(f"Unsupported Zarr storage mode: {reader.encoding}.")
        self._storage_mode = reader.encoding
        self.intensity_dtype = reader.intensity_dtype
        self.mz_dtype = reader.mz_dtype

    def _target_indices(self) -> np.ndarray:
        coords = self._require_reader().read_coordinates()
        xs = coords[:, 0].astype(np.int64, copy=False)
        ys = coords[:, 1].astype(np.int64, copy=False)
        mask = (
            (xs >= self.target_locs[0][0])
            & (xs <= self.target_locs[1][0])
            & (ys >= self.target_locs[0][1])
            & (ys <= self.target_locs[1][1])
        )
        return np.nonzero(mask)[0].astype(np.int64, copy=False)
