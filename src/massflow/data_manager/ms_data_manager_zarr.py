"""MS data manager for MassFlow spectra Zarr stores."""

from __future__ import annotations

from os import PathLike
from typing import Any, Generator

import numpy as np

from massflow.data_manager.ms_data_manager import MSDataManager
from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.module.spectrum_zarr import SpectrumZarr
from massflow.tools.zarr_io import (
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
    ZarrStorageMode,
    MSIZarrReader,
    MSIZarrWriter,
    storage_mode_from_metadata,
)
from massflow.tools.infer_spectrum_type import (
    SpectrumType,
    set_spectrum_type_metadata,
)
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
        chunk_size: int = 256,
        coordinate_base: int | None = None,
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
        self._input_coordinate_base = coordinate_base
        self._coordinate_base = "zero" if mode == "r" else None
        self.ms.meta = MetaDataZarr()
        # storage_mode lives on meta only; encode the constructor request into meta now
        # so write paths can recover it via storage_mode_from_metadata, and copy_meta
        # can preserve it across a meta replacement.
        if storage_mode is not None:
            _set_storage_mode_metadata(self.ms.meta, storage_mode)
        self._storage_mode: str | None = None
        self._finalized = False
        self._reader: MSIZarrReader | None = None

    @property
    def zarr_filepath(self) -> str:
        """Return full .zarr path derived from the base path."""
        return f"{self.file_base_path}.zarr"

    @property
    def storage_mode(self) -> str:
        if self._storage_mode is None:
            self._open_arrays_for_read()
        if self._storage_mode is None:
            raise ValueError("Zarr storage mode has not been initialized.")
        return self._storage_mode

    @property
    def writer(self) -> MSIZarrWriter:
        """Return the streaming Zarr writer."""
        if self.ms.meta is None:
            raise ValueError(
                "MS meta data is None. Please load or copy meta data first."
            )
        if self._writer is None:
            self._writer = MSIZarrWriter(
                self.zarr_filepath,
                row_axis="pixel",
                encoding=storage_mode_from_metadata(self.ms.meta),
                metadata=self.ms.meta,
                width=getattr(self.ms.meta, "max_count_of_pixels_x", None),
                height=getattr(self.ms.meta, "max_count_of_pixels_y", None),
                intensity_dtype=self.intensity_dtype,
                target_chunk_bytes=4 * 1024 * 1024,
                add_buffer_size=self._chunk_size,
                temp_dir=self._temp_dir,
                coordinate_base=self._input_coordinate_base,
            )
        return self._writer

    @property
    def reader(self) -> MSIZarrReader:
        if self._reader is None:
            self._reader = MSIZarrReader(self.zarr_filepath)
        return self._reader

    def copy_meta(self, source_dm: MSDataManager):
        """Copy metadata from another data manager, dropping fields that may change here."""
        # Preserve a storage_mode that was already declared on this manager (typically
        # set in __init__ via the storage_mode argument); the source's value is dropped
        # below because this manager owns its own write/finalize.
        declared_mode = _meta_storage_mode(self.ms.meta) if self.ms.meta else None
        self.ms.meta = MetaDataZarr.from_metadata(source_dm.ms.meta)
        # These fields are determined by this manager's own write/finalize, not the source.
        for key in ("continuous", "processed", "spectrum_count_num"):
            self.ms.meta._meta.pop(key, None)
        if declared_mode is not None:
            _set_storage_mode_metadata(self.ms.meta, declared_mode)

    def extract_metadata(self):
        """Load metadata from the Zarr metadata group."""
        self.ms.meta = self.reader.metadata

    def preload_hook(self, *args, **kwargs):
        """Load metadata and shared m/z axis before building spectrum placeholders."""
        self.extract_metadata()
        self._open_arrays_for_read()
        if self.storage_mode == ZARR_MODE_CONTINUOUS and self.ms.shared_mz_list is None:
            assert self.reader.shared_mz_axis is not None
            self.ms.shared_mz_list = np.asarray(
                self.reader.shared_mz_axis[:], dtype=self.mz_dtype
            )

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
        coordinates = np.asarray(self.reader.coordinates[:])
        for index, coord in enumerate(coordinates):
            x, y = int(coord[0]), int(coord[1])
            if (
                self.target_locs[0][0] <= x <= self.target_locs[1][0]
                and self.target_locs[0][1] <= y <= self.target_locs[1][1]
            ):
                spectrum = SpectrumZarr(
                    reader=self.reader,
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
        self.writer.write_flat_batch(
            intensity_flat=intensity_flat,
            lengths=lengths,
            coordinates=coordinates,
            mz_axis=mz_axis,
        )

    def finalize(self) -> None:
        """Finalize the Zarr writer."""
        if self.mode != "w" or self._finalized:
            return
        self.writer.finalize()
        self._finalized = True

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
            lengths = self._batch_lengths(batch_indices)
            coordinates = np.asarray(
                self.reader.coordinates.get_orthogonal_selection(
                    (batch_indices, slice(None))
                )
            )
            intensity = self._read_flat_values(self.reader.intensity, batch_indices)

            if self.storage_mode == ZARR_MODE_CONTINUOUS:
                mz_data = self._shared_mz_axis()
            elif include_mz:
                assert self.reader.point_mz is not None
                mz_data = self._read_flat_values(self.reader.point_mz, batch_indices)
            else:
                mz_data = None

            yield mz_data, intensity, lengths, coordinates

    def close(self):
        """Close the Zarr reader and writer."""
        if self._reader is not None:
            self._reader.close()
        super().close()

    @classmethod
    def from_flat_output_contract(
        cls,
        data_manager,
        filepath: str | PathLike[str],
        *,
        output_spectrum_type: SpectrumType = "centroid",
        storage_mode: ZarrStorageMode,
        temp_dir: str | PathLike[str] | None = None,
        mz_dtype: Any | None = None,
        intensity_dtype: Any | None = None,
        chunk_size: int = 256,
    ) -> "MSDataManagerZarr":
        output_dm = cls(
            filepath=filepath,
            mode="w",
            storage_mode=storage_mode,
            temp_dir=temp_dir,
            mz_dtype=data_manager.mz_dtype if mz_dtype is None else mz_dtype,
            intensity_dtype=(
                data_manager.intensity_dtype
                if intensity_dtype is None
                else intensity_dtype
            ),
            chunk_size=chunk_size,
            coordinate_base=_manager_coordinate_base(data_manager),
        )
        output_dm.copy_meta(data_manager)
        set_spectrum_type_metadata(output_dm.ms.meta, output_spectrum_type)
        return output_dm

    @classmethod
    def from_data_manager(
        cls,
        data_manager,
        filepath: str | PathLike[str],
        *,
        temp_dir: str | PathLike[str] | None = None,
        mz_dtype: Any | None = None,
        intensity_dtype: Any | None = None,
        chunk_size: int = 256,
        batch_size: int = 256,
        storage_mode: ZarrStorageMode | None = None,
    ) -> "MSDataManagerZarr":
        if storage_mode is None:
            storage_mode = storage_mode_from_metadata(data_manager.ms.meta)
        output_dm = cls(
            filepath=filepath,
            mode="w",
            storage_mode=storage_mode,
            temp_dir=temp_dir,
            mz_dtype=data_manager.mz_dtype if mz_dtype is None else mz_dtype,
            intensity_dtype=(
                data_manager.intensity_dtype
                if intensity_dtype is None
                else intensity_dtype
            ),
            chunk_size=chunk_size,
            coordinate_base=_manager_coordinate_base(data_manager),
        )
        output_dm.copy_meta(data_manager)
        for (
            mz_axis,
            intensity_flat,
            lengths,
            coordinates,
        ) in data_manager.flat_generator(
            batch_size=batch_size,
            include_mz=True,
        ):
            output_dm.write_flat_batch(
                intensity_flat=intensity_flat,
                lengths=lengths,
                coordinates=coordinates,
                mz_axis=mz_axis,
            )
        output_dm.finalize()
        return output_dm

    def _shared_mz_axis(self) -> np.ndarray:
        if self.ms.shared_mz_list is None:
            assert self.reader.shared_mz_axis is not None
            self.ms.shared_mz_list = np.asarray(
                self.reader.shared_mz_axis[:], dtype=self.mz_dtype
            )
        return np.asarray(self.ms.shared_mz_list, dtype=self.mz_dtype)

    def _open_arrays_for_read(self) -> None:
        if self.reader.row_axis != "pixel":
            raise ValueError("MSDataManagerZarr requires a pixel-major Zarr store.")
        self._storage_mode = self.reader.encoding
        self.intensity_dtype = np.dtype(self.reader.intensity.dtype)
        if self.reader.shared_mz_axis is not None:
            self.mz_dtype = np.dtype(self.reader.shared_mz_axis.dtype)
        elif self.reader.point_mz is not None:
            self.mz_dtype = np.dtype(self.reader.point_mz.dtype)

    def _target_indices(self) -> np.ndarray:
        coords = np.asarray(self.reader.coordinates[:])
        xs = coords[:, 0].astype(np.int64, copy=False)
        ys = coords[:, 1].astype(np.int64, copy=False)
        mask = (
            (xs >= self.target_locs[0][0])
            & (xs <= self.target_locs[1][0])
            & (ys >= self.target_locs[0][1])
            & (ys <= self.target_locs[1][1])
        )
        return np.nonzero(mask)[0].astype(np.int64, copy=False)

    def _batch_lengths(self, indices: np.ndarray) -> np.ndarray:
        if _is_contiguous_indices(indices):
            start = int(indices[0])
            end = int(indices[-1]) + 2
            return np.diff(np.asarray(self.reader.offsets[start:end])).astype("<u4")
        starts = np.asarray(self.reader.offsets.get_orthogonal_selection((indices,)))
        ends = np.asarray(self.reader.offsets.get_orthogonal_selection((indices + 1,)))
        return (ends - starts).astype("<u4")

    def _read_flat_values(self, array: Any, indices: np.ndarray) -> np.ndarray:
        if _is_contiguous_indices(indices):
            start = int(self.reader.offsets[int(indices[0])])
            end = int(self.reader.offsets[int(indices[-1]) + 1])
            return np.asarray(array[start:end])
        parts = []
        for index in indices:
            start = int(self.reader.offsets[int(index)])
            end = int(self.reader.offsets[int(index) + 1])
            parts.append(np.asarray(array[start:end]))
        if not parts:
            return np.array([], dtype=array.dtype)
        return np.concatenate(parts)


def _is_contiguous_indices(indices: np.ndarray) -> bool:
    return bool(indices.size > 0 and np.all(np.diff(indices) == 1))


def _set_storage_mode_metadata(meta: MetaDataZarr, storage_mode: str) -> None:
    if storage_mode == ZARR_MODE_CONTINUOUS:
        meta.continuous = True
        meta.processed = None
        return
    if storage_mode == ZARR_MODE_PROCESSED:
        meta.continuous = None
        meta.processed = True
        return
    raise ValueError("storage_mode must be 'continuous' or 'processed'.")


def _meta_storage_mode(meta: Any) -> ZarrStorageMode | None:
    """Read storage_mode declared on a meta object, or None if not set."""
    if meta is None:
        return None
    if getattr(meta, "continuous", None) is True:
        return ZARR_MODE_CONTINUOUS
    if getattr(meta, "processed", None) is True:
        return ZARR_MODE_PROCESSED
    return None


def _manager_coordinate_base(data_manager: Any) -> int | None:
    base = getattr(data_manager, "_coordinate_base", None)
    if base in {0, "zero"}:
        return 0
    if base in {1, "one"}:
        return 1
    return None
