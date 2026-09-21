"""Writer for the MassFlow MSI Zarr schema."""

# pylint: disable=too-many-lines

from __future__ import annotations

import os
from concurrent.futures import Future
from contextlib import nullcontext
from pathlib import Path
from typing import cast

import numpy as np
import zarr

from massflow.concurrency import get_effective_budget
from massflow.concurrency.runtime import io_executor

from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.msi_zarr.schema import (
    ZARR_FORMAT_NAME,
    ZARR_FORMAT_VERSION,
    ZARR_GROUP_ION_IMAGE,
    ZARR_GROUP_SPECTRA,
    ZARR_MODE_AUTO,
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
    ZarrRowAxis,
    ZarrStorageMode,
    ZarrWriteMode,
)
from massflow.tools.ion_major_temp_store import BatchInfo, IonMajorTempStore
from massflow.tools.logger import get_logger

logger = get_logger("massflow.msi_zarr.zarr_writer")


class MSIZarrWriter:
    """Write pixel-major or ion-major MassFlow MSI Zarr stores."""

    def __init__(
        self,
        path: str | os.PathLike[str],
        *,
        row_axis: ZarrRowAxis,
        encoding: ZarrWriteMode,
        metadata: MetaDataZarr | None,
        width: int | None = None,
        height: int | None = None,
        mz_dtype: np.dtype | type = np.float64,
        intensity_dtype: np.dtype | type = np.float32,
        target_chunk_bytes: int = 4 * 1024 * 1024,
        add_buffer_size: int = 256,
        temp_dir: str | os.PathLike[str] | None = None,
        temp_workers: int = 4,
        ion_block: int = 256,
        compute_mean_spectrum: bool = False,
        compute_tic: bool = False,
        compute_rms: bool = False,
        include_spectra: bool = False,
    ) -> None:
        if row_axis not in {"pixel", "ion"}:
            raise ValueError("row_axis must be 'pixel' or 'ion'.")
        if encoding not in {
            ZARR_MODE_AUTO,
            ZARR_MODE_CONTINUOUS,
            ZARR_MODE_PROCESSED,
        }:
            raise ValueError(
                "encoding must be 'auto', 'continuous', or 'processed'."
            )
        if include_spectra and (
            row_axis != "ion" or encoding != ZARR_MODE_CONTINUOUS
        ):
            raise ValueError(
                "include_spectra=True requires row_axis='ion' and "
                "continuous encoding."
            )
        if row_axis == "ion" and encoding != ZARR_MODE_CONTINUOUS:
            raise ValueError("ion-major output requires continuous encoding.")
        if include_spectra:
            compute_tic = True
            compute_mean_spectrum = True
        elif row_axis == "ion" and compute_tic:
            raise ValueError("compute_tic is only supported for pixel-major output.")
        elif row_axis == "ion" and compute_rms:
            raise ValueError("compute_rms is only supported for pixel-major output.")
        self.path = _normalize_zarr_path(path)
        self.row_axis: ZarrRowAxis = row_axis
        self.encoding: ZarrWriteMode = encoding
        self.include_spectra = bool(include_spectra)
        self.meta = MetaDataZarr() if metadata is None else metadata
        self.width = None if width is None else int(width)
        self.height = None if height is None else int(height)
        self.mz_dtype = _normalize_mz_dtype(mz_dtype)
        self.intensity_dtype = np.dtype(intensity_dtype)
        self.target_chunk_bytes = max(1, int(target_chunk_bytes))
        self.add_buffer_size = max(1, int(add_buffer_size))
        self.ion_block = max(1, int(ion_block))
        self.compute_mean_spectrum = bool(compute_mean_spectrum)
        self.compute_tic = bool(compute_tic)
        self.compute_rms = bool(compute_rms)
        self._temp_dir = temp_dir
        self._temp_workers = max(1, int(temp_workers))
        self._root: zarr.Group | None = None
        self._data: zarr.Group | None = None
        self._intensity: zarr.Array | None = None
        self._offsets: zarr.Array | None = None
        self._point_mz: zarr.Array | None = None
        self._spectra_intensity: zarr.Array | None = None
        self._spectra_offsets: zarr.Array | None = None
        self._shared_mz: np.ndarray | None = None
        self._coordinates: list[np.ndarray] = []
        self._temp_store: IonMajorTempStore | None = None
        self._total_points = 0
        self._spectra_total_points = 0
        self._n_input_pixels = 0
        self._finalized = False
        self._add_intensity: list[np.ndarray] = []
        self._add_mz: list[np.ndarray] = []
        self._add_coordinates: list[np.ndarray] = []
        self._pending_batches: list[
            tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        ] = []
        self._tic_batches: list[np.ndarray] = []
        self._rms_batches: list[np.ndarray] = []

    # --- Public write lifecycle ----------------------------------------------

    def add_spectrum(
        self,
        intensity: np.ndarray,
        coordinates: np.ndarray | list[int] | tuple[int, ...],
        mz: np.ndarray,
    ) -> None:
        intensity = np.asarray(intensity, dtype=self.intensity_dtype)
        mz = np.asarray(mz, dtype=self.mz_dtype)
        if intensity.ndim != 1 or mz.ndim != 1 or intensity.size != mz.size:
            raise ValueError("intensity and mz must be equal-length 1-D arrays.")
        self._add_intensity.append(intensity.copy())
        self._add_mz.append(mz.copy())
        self._add_coordinates.append(np.asarray(coordinates).copy())
        if len(self._add_intensity) >= self.add_buffer_size:
            self._flush_add_buffer()

    def write_flat_batch(
        self,
        intensity_flat: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
        mz_axis: np.ndarray | None,
    ) -> None:
        if self._finalized:
            raise ValueError("Cannot write after finalize().")
        lengths = np.asarray(lengths, dtype=np.int64)
        coordinates = _ensure_xyz_coordinates(coordinates)
        intensity = np.asarray(intensity_flat, dtype=self.intensity_dtype)
        mz = np.asarray(mz_axis, dtype=self.mz_dtype)
        if lengths.ndim != 1 or intensity.ndim != 1 or mz.ndim != 1:
            raise ValueError("intensity_flat, lengths, and mz_axis must be 1-D.")
        if lengths.size == 0:
            return
        if np.any(lengths <= 0):
            raise ValueError("Writer does not accept empty spectra.")
        if coordinates.shape[0] != lengths.size:
            raise ValueError("coordinates row count must equal lengths size.")
        if intensity.size != int(np.sum(lengths, dtype=np.int64)):
            raise ValueError("intensity_flat size must equal sum(lengths).")

        if self.encoding == ZARR_MODE_AUTO:
            detected = _infer_encoding(mz, lengths)
            if detected is None:
                if self._pending_batches:
                    first_mz = self._pending_batches[0][3]
                    self.encoding = (
                        ZARR_MODE_CONTINUOUS
                        if np.array_equal(mz, first_mz)
                        else ZARR_MODE_PROCESSED
                    )
                    self._flush_pending_batches()
                else:
                    self._pending_batches.append(
                        (
                            intensity.copy(),
                            lengths.copy(),
                            coordinates.copy(),
                            mz.copy(),
                        )
                    )
                    return
            else:
                self.encoding = detected
                self._flush_pending_batches()

        self._set_or_validate_shared_mz(mz, lengths)
        self._coordinates.append(coordinates.copy())
        self._n_input_pixels += int(lengths.size)
        if self.row_axis == "pixel":
            self._write_pixel_batch(intensity, lengths, mz)
        else:
            self._write_ion_batch(intensity, lengths)

    def finalize(self) -> None:
        if self._finalized:
            return
        self._flush_add_buffer()
        if self.encoding == ZARR_MODE_AUTO:
            if not self._pending_batches:
                raise ValueError("No spectra were written.")
            self.encoding = ZARR_MODE_PROCESSED
            self._flush_pending_batches()
        if self._n_input_pixels == 0:
            raise ValueError("No spectra were written.")
        coordinates = np.vstack(self._coordinates)
        if self.row_axis == "pixel":
            coordinate_index = prepare_pixel_coordinate_index(
                coordinates,
                width=self.width,
                height=self.height,
            )
        else:
            coordinate_index = prepare_coordinate_index(
                coordinates,
                width=self.width,
                height=self.height,
            )
        if coordinate_index.source_idx.size != coordinates.shape[0]:
            logger.warning(
                "Removed duplicate coordinates while keeping the last pixel: "
                f"input={coordinates.shape[0]}, output={coordinate_index.source_idx.size}, "
                f"removed={coordinates.shape[0] - coordinate_index.source_idx.size}."
            )
        self.width = coordinate_index.width
        self.height = coordinate_index.height

        if self.row_axis == "pixel":
            self._finalize_pixel(coordinate_index)
        else:
            self._finalize_ion(coordinate_index)
        assert self._root is not None
        self.meta.max_count_of_pixels_x = self.width
        self.meta.max_count_of_pixels_y = self.height
        self.meta.spectrum_count_num = int(coordinate_index.coordinates.shape[0])
        if self.encoding == ZARR_MODE_CONTINUOUS:
            self.meta.continuous = True
            self.meta.processed = None
        else:
            self.meta.continuous = None
            self.meta.processed = True
        self.meta.to_zarr(self._root)
        self._write_root_attrs(state="writing")
        self._validate_written_structure()
        self._root.attrs["write_state"] = "complete"
        self._finalized = True
        if self._temp_store is not None:
            self._temp_store.close()
            self._temp_store = None

    def close(self) -> None:
        if self._temp_store is not None:
            self._temp_store.close()
            self._temp_store = None
        self._root = None
        self._data = None
        self._intensity = None
        self._offsets = None
        self._point_mz = None
        self._spectra_intensity = None
        self._spectra_offsets = None
        self._pending_batches.clear()
        self._tic_batches.clear()
        self._rms_batches.clear()

    # --- Batch ingestion ------------------------------------------------------

    def _write_pixel_batch(
        self,
        intensity: np.ndarray,
        lengths: np.ndarray,
        mz: np.ndarray,
    ) -> None:
        self._ensure_pixel_arrays(lengths)
        assert self._intensity is not None and self._offsets is not None
        self._intensity.append(intensity, axis=0)
        if self.encoding == ZARR_MODE_PROCESSED:
            assert self._point_mz is not None
            self._point_mz.append(mz, axis=0)
        self._total_points += int(intensity.size)
        batch_offsets = self._total_points - intensity.size + np.cumsum(
            lengths, dtype=np.int64
        )
        self._offsets.append(batch_offsets.astype("<i8", copy=False), axis=0)
        self._accumulate_spectrum_stats(intensity, lengths)

    def _write_ion_batch(
        self,
        intensity: np.ndarray,
        lengths: np.ndarray,
    ) -> None:
        assert self._shared_mz is not None
        if self.include_spectra:
            self._append_spectra_batch(intensity, lengths)
        if self._temp_store is None:
            self._temp_store = IonMajorTempStore(
                temp_dir=self._temp_dir,
                dtype=self.intensity_dtype,
                workers=self._temp_workers,
            )
        temp_store = self._temp_store
        temp_store.append_batch(
            intensity,
            n_pixels=int(lengths.size),
            n_ions=int(self._shared_mz.size),
        )

    # --- Combined /spectra streaming -----------------------------------------

    def _ensure_spectra_arrays(self, lengths: np.ndarray) -> None:
        """Create the v1.1 root and the streaming /spectra group on first write."""
        if self._spectra_intensity is not None:
            return
        if self._root is None:
            Path(self.path).parent.mkdir(parents=True, exist_ok=True)
            self._root = zarr.open_group(self.path, mode="w")
            self._write_root_attrs(state="writing")
        root = self._root
        root.require_group("axes")
        spectra = root.create_group(ZARR_GROUP_SPECTRA)
        spectra.attrs.update(
            {"row_axis": "pixel", "encoding": ZARR_MODE_CONTINUOUS}
        )
        assert self._shared_mz is not None
        row_length = int(self._shared_mz.size)
        chunk_points = _chunk_points(
            row_length,
            self.intensity_dtype.itemsize,
            self.target_chunk_bytes,
        )
        self._spectra_intensity = spectra.create_array(
            "intensity",
            shape=(0,),
            chunks=(chunk_points,),
            dtype=self.intensity_dtype,
        )
        self._spectra_offsets = spectra.create_array(
            "offsets",
            data=np.array([0], dtype="<i8"),
            chunks=(max(2, min(4096, int(lengths.size) + 1)),),
        )

    def _append_spectra_batch(
        self,
        intensity: np.ndarray,
        lengths: np.ndarray,
    ) -> None:
        """Append one pixel-major batch and accumulate enabled spectrum stats."""
        self._ensure_spectra_arrays(lengths)
        assert self._spectra_intensity is not None
        assert self._spectra_offsets is not None
        self._spectra_intensity.append(intensity, axis=0)
        self._spectra_total_points += int(intensity.size)
        batch_offsets = (
            self._spectra_total_points - intensity.size
            + np.cumsum(lengths, dtype=np.int64)
        )
        self._spectra_offsets.append(
            batch_offsets.astype("<i8", copy=False), axis=0
        )
        self._accumulate_spectrum_stats(intensity, lengths)

    def _finalize_spectra(self, source_idx: np.ndarray) -> None:
        """Compact streamed /spectra and write enabled spectrum stats."""
        if source_idx.size != self._n_input_pixels:
            self._compact_spectra_arrays(source_idx)
            tic = (
                np.concatenate(self._tic_batches)[source_idx]
                if self._tic_batches
                else np.empty(0, dtype=np.float64)
            )
            rms = (
                np.concatenate(self._rms_batches)[source_idx]
                if self._rms_batches
                else np.empty(0, dtype=np.float64)
            )
        else:
            tic = None
            rms = None
        self._write_tic(tic)
        self._write_rms(rms)

    def _compact_spectra_arrays(self, source_idx: np.ndarray) -> None:
        assert self._spectra_intensity is not None
        assert self._spectra_offsets is not None
        # Single bulk read from zarr, gather in numpy, single write back.
        # Avoids per-row zarr sync overhead (and potential deadlocks with
        # zarr v3's asyncio-based sync in constrained runtimes).
        old_offsets = np.asarray(self._spectra_offsets[:], dtype=np.int64)
        starts = old_offsets[source_idx]
        lengths = old_offsets[source_idx + 1] - starts
        total = int(lengths.sum())
        gather_idx = np.repeat(starts, lengths) + np.concatenate(
            [np.arange(n) for n in lengths]
        )
        all_intensity = np.asarray(self._spectra_intensity[:])
        compacted = all_intensity[gather_idx]
        new_offsets = np.empty(source_idx.size + 1, dtype="<i8")
        np.cumsum(lengths, out=new_offsets[1:])
        new_offsets[0] = 0
        self._spectra_intensity.resize((total,))
        self._spectra_intensity[:] = compacted
        self._spectra_offsets.resize(new_offsets.shape)
        self._spectra_offsets[:] = new_offsets
        self._spectra_total_points = total

    # --- Pixel-major finalization --------------------------------------------

    def _ensure_pixel_arrays(self, lengths: np.ndarray) -> None:
        if self._root is not None:
            return
        Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        root = zarr.open_group(self.path, mode="w")
        self._root = root
        self._write_root_attrs(state="writing")
        axes = root.create_group("axes")
        data = root.create_group(ZARR_GROUP_SPECTRA)
        self._data = data
        data.attrs.update(
            {"row_axis": self.row_axis, "encoding": self.encoding}
        )
        row_length = (
            int(self._shared_mz.size)
            if self.encoding == ZARR_MODE_CONTINUOUS
            and self._shared_mz is not None
            else None
        )
        chunk_points = _chunk_points(
            row_length,
            self.intensity_dtype.itemsize,
            self.target_chunk_bytes,
        )
        self._intensity = data.create_array(
            "intensity",
            shape=(0,),
            chunks=(chunk_points,),
            dtype=self.intensity_dtype,
        )
        self._offsets = data.create_array(
            "offsets",
            data=np.array([0], dtype="<i8"),
            chunks=(max(2, min(4096, int(lengths.size) + 1)),),
        )
        if self.encoding == ZARR_MODE_CONTINUOUS:
            assert self._shared_mz is not None
            mz_array = axes.create_array("mz", data=self._shared_mz)
            mz_array.attrs["unit"] = "m/z"
        else:
            point_mz = data.create_array(
                "mz",
                shape=(0,),
                chunks=(chunk_points,),
                dtype=self.mz_dtype,
            )
            point_mz.attrs["unit"] = "m/z"
            self._point_mz = point_mz

    def _finalize_pixel(self, coordinate_index: "CoordinateIndex") -> None:
        assert self._root is not None
        if coordinate_index.source_idx.size != self._n_input_pixels:
            self._compact_pixel_arrays(coordinate_index.source_idx)
        axes = cast(zarr.Group, self._root["axes"])
        axes.create_array(
            "coordinates",
            data=coordinate_index.coordinates.astype("<u4", copy=False),
            overwrite=True,
        )
        self._write_tic()
        self._write_rms()

    def _accumulate_spectrum_stats(
        self,
        intensity: np.ndarray,
        lengths: np.ndarray,
    ) -> None:
        if not self.compute_tic and not self.compute_rms:
            return
        starts = np.empty(lengths.size, dtype=np.int64)
        starts[0] = 0
        if lengths.size > 1:
            np.cumsum(lengths[:-1], dtype=np.int64, out=starts[1:])
        work = intensity.astype(np.float64)
        if self.compute_tic:
            self._tic_batches.append(np.add.reduceat(work, starts))
        if self.compute_rms:
            np.square(work, out=work)
            sum_squares = np.add.reduceat(work, starts)
            self._rms_batches.append(np.sqrt(sum_squares / lengths))

    def _write_tic(self, tic: np.ndarray | None = None) -> None:
        if not self.compute_tic:
            return
        assert self._root is not None
        if tic is None:
            tic = (
                np.concatenate(self._tic_batches)
                if self._tic_batches
                else np.empty(0, dtype=np.float64)
            )
        stats = self._root.require_group("stats")
        array = stats.create_array("tic", data=tic, overwrite=True)
        array.attrs.update(
            {
                "source": f"/{ZARR_GROUP_SPECTRA}",
                "reduced_axis": "mz",
            }
        )

    def _write_rms(self, rms: np.ndarray | None = None) -> None:
        if not self.compute_rms:
            return
        assert self._root is not None
        if rms is None:
            rms = (
                np.concatenate(self._rms_batches)
                if self._rms_batches
                else np.empty(0, dtype=np.float64)
            )
        stats = self._root.require_group("stats")
        array = stats.create_array("rms", data=rms, overwrite=True)
        array.attrs.update(
            {
                "source": f"/{ZARR_GROUP_SPECTRA}",
                "reduced_axis": "mz",
                "denominator": "spectrum_points",
            }
        )

    def _compact_pixel_arrays(self, source_idx: np.ndarray) -> None:
        assert self._intensity is not None and self._offsets is not None
        old_offsets = np.asarray(self._offsets[:], dtype=np.int64)
        starts = old_offsets[source_idx]
        lengths = old_offsets[source_idx + 1] - starts
        total = int(lengths.sum())
        gather_idx = np.repeat(starts, lengths) + np.concatenate(
            [np.arange(n) for n in lengths]
        )
        new_offsets = np.empty(source_idx.size + 1, dtype="<i8")
        np.cumsum(lengths, out=new_offsets[1:])
        new_offsets[0] = 0
        # Bulk read → gather → write for intensity
        all_intensity = np.asarray(self._intensity[:])
        self._intensity.resize((total,))
        self._intensity[:] = all_intensity[gather_idx]
        # Same for mz if present
        if self._point_mz is not None:
            all_mz = np.asarray(self._point_mz[:])
            self._point_mz.resize((total,))
            self._point_mz[:] = all_mz[gather_idx]
        self._offsets.resize(new_offsets.shape)
        self._offsets[:] = new_offsets
        self._total_points = total

    # --- Ion-major writing ----------------------------------------------------

    def _finalize_ion(self, coordinate_index: "CoordinateIndex") -> None:
        if self._temp_store is None or self._shared_mz is None:
            raise ValueError("Ion-major writer has no temporary data.")
        temp_store = self._temp_store
        shared_mz = self._shared_mz
        batches = temp_store.finish()
        if self._root is None:
            Path(self.path).parent.mkdir(parents=True, exist_ok=True)
            self._root = zarr.open_group(self.path, mode="w")
            self._write_root_attrs(state="writing")
        root = self._root
        axes = root.require_group("axes")
        axes.create_array(
            "coordinates",
            data=coordinate_index.coordinates.astype("<u4", copy=False),
            overwrite=True,
        )
        mz_array = axes.create_array("mz", data=shared_mz, overwrite=True)
        mz_array.attrs["unit"] = "m/z"
        data = root.require_group(ZARR_GROUP_ION_IMAGE)
        self._data = data
        data.attrs.update(
            {"row_axis": self.row_axis, "encoding": self.encoding}
        )
        n_pixels = int(coordinate_index.coordinates.shape[0])
        n_ions = int(shared_mz.size)
        rows_per_chunk = _rows_per_chunk(
            n_pixels,
            self.intensity_dtype.itemsize,
            self.target_chunk_bytes,
        )
        chunk_points = max(1, rows_per_chunk * n_pixels)
        intensity_array = data.create_array(
            "intensity",
            shape=(n_ions * n_pixels,),
            chunks=(chunk_points,),
            dtype=self.intensity_dtype,
            overwrite=True,
        )
        self._intensity = intensity_array
        offsets = np.arange(n_ions + 1, dtype=np.int64) * n_pixels
        self._offsets = data.create_array(
            "offsets", data=offsets.astype("<i8"), overwrite=True
        )
        selections = _split_source_indices(coordinate_index.source_idx, batches)
        merge_rows = max(rows_per_chunk, self.ion_block)
        merge_rows = ((merge_rows + rows_per_chunk - 1) // rows_per_chunk) * rows_per_chunk
        mean_spectrum = (
            np.empty(n_ions, dtype=np.float64)
            if self.compute_mean_spectrum
            else None
        )

        def assemble_band(ion_start: int, ion_end: int) -> np.ndarray:
            band = np.empty((ion_end - ion_start, n_pixels), dtype=self.intensity_dtype)
            output_start = 0
            for batch, local_indices in zip(batches, selections):
                source = temp_store.read_ion_block(batch, ion_start, ion_end)
                output_end = output_start + local_indices.size
                band[:, output_start:output_end] = source[:, local_indices]
                output_start = output_end
            return band

        starts = list(range(0, n_ions, merge_rows))
        prefetch = get_effective_budget() != 1
        with io_executor(max_workers=1) if prefetch else nullcontext() as executor:
            future: Future[np.ndarray] | None = None
            for index, ion_start in enumerate(starts):
                ion_end = min(n_ions, ion_start + merge_rows)
                band = (
                    assemble_band(ion_start, ion_end)
                    if future is None
                    else future.result()
                )
                future = None
                if executor is not None and index + 1 < len(starts):
                    next_start = starts[index + 1]
                    next_end = min(n_ions, next_start + merge_rows)
                    future = executor.submit(assemble_band, next_start, next_end)
                point_start = ion_start * n_pixels
                point_end = ion_end * n_pixels
                intensity_array[point_start:point_end] = band.reshape(-1)
                if mean_spectrum is not None:
                    mean_spectrum[ion_start:ion_end] = np.mean(
                        band, axis=1, dtype=np.float64
                    )
        self._write_mean_spectrum(mean_spectrum)
        if self.include_spectra:
            self._finalize_spectra(coordinate_index.source_idx)

    def _write_mean_spectrum(
        self,
        mean_spectrum: np.ndarray | None,
    ) -> None:
        if mean_spectrum is None:
            return
        assert self._root is not None
        stats = self._root.require_group("stats")
        array = stats.create_array(
            "mean_spectrum", data=mean_spectrum, overwrite=True
        )
        array.attrs.update(
            {
                "reduced_axis": "pixel",
                "denominator": "sampled_pixels",
                "source": f"/{ZARR_GROUP_ION_IMAGE}",
            }
        )

    # --- Buffered input and shared-axis handling -----------------------------

    def _set_or_validate_shared_mz(
        self,
        mz: np.ndarray,
        lengths: np.ndarray,
    ) -> None:
        if self.encoding == ZARR_MODE_PROCESSED:
            if self.row_axis != "pixel":
                raise ValueError("Processed encoding is only valid for pixel-major.")
            if mz.size != int(np.sum(lengths, dtype=np.int64)):
                raise ValueError("Processed mz size must equal sum(lengths).")
            return
        if self._shared_mz is None:
            self._shared_mz = mz.copy()
        elif not np.array_equal(mz, self._shared_mz):
            raise ValueError("All continuous batches must share the same m/z axis.")
        assert self._shared_mz is not None
        if np.any(lengths != self._shared_mz.size):
            raise ValueError("Every continuous spectrum must match axes/mz length.")

    def _flush_add_buffer(self) -> None:
        if not self._add_intensity:
            return
        lengths = np.asarray([row.size for row in self._add_intensity], dtype="<u4")
        intensity = np.concatenate(self._add_intensity)
        coordinates = _ensure_xyz_coordinates(np.vstack(self._add_coordinates))
        shared_mz = all(
            np.array_equal(row, self._add_mz[0]) for row in self._add_mz[1:]
        )
        if self.encoding == ZARR_MODE_CONTINUOUS:
            if not shared_mz:
                raise ValueError("Continuous add_spectrum calls must share m/z.")
            mz = self._add_mz[0]
        elif self.encoding == ZARR_MODE_AUTO and shared_mz:
            mz = self._add_mz[0]
        else:
            mz = np.concatenate(self._add_mz)
        self._add_intensity.clear()
        self._add_mz.clear()
        self._add_coordinates.clear()
        self.write_flat_batch(intensity, lengths, coordinates, mz)

    def _flush_pending_batches(self) -> None:
        pending = self._pending_batches
        self._pending_batches = []
        for intensity, lengths, coordinates, mz in pending:
            self.write_flat_batch(intensity, lengths, coordinates, mz)

    # --- Metadata and structural validation ----------------------------------

    def _write_root_attrs(self, *, state: str) -> None:
        assert self._root is not None
        if self.width is None or self.height is None:
            spatial_shape = [0, 0]
        else:
            spatial_shape = [int(self.height), int(self.width)]
        self._root.attrs.update(
            {
                "format": ZARR_FORMAT_NAME,
                "format_version": ZARR_FORMAT_VERSION,
                "write_state": state,
                "coordinate_order": ["x", "y", "z"],
                "coordinate_base": 1,
                "spatial_shape": spatial_shape,
            }
        )

    def _validate_written_structure(self) -> None:
        assert self._root is not None
        axes = cast(zarr.Group, self._root["axes"])
        coordinates = cast(zarr.Array, axes["coordinates"])
        if self.row_axis == "pixel":
            data = cast(zarr.Group, self._root[ZARR_GROUP_SPECTRA])
            offsets = cast(zarr.Array, data["offsets"])
            intensity = cast(zarr.Array, data["intensity"])
            if (
                _read_array_int(offsets, 0) != 0
                or _read_array_int(offsets, -1) != int(intensity.shape[0])
            ):
                raise ValueError("Written offsets do not match intensity.")
            n_rows = int(offsets.shape[0]) - 1
            if n_rows != int(coordinates.shape[0]):
                raise ValueError(
                    "Written pixel row count does not match coordinates."
                )
            if self.encoding == ZARR_MODE_PROCESSED:
                point_mz = cast(zarr.Array, data["mz"])
                if int(point_mz.shape[0]) != int(intensity.shape[0]):
                    raise ValueError(
                        "Written spectra/mz does not match intensity."
                    )
            return
        data = cast(zarr.Group, self._root[ZARR_GROUP_ION_IMAGE])
        offsets = cast(zarr.Array, data["offsets"])
        intensity = cast(zarr.Array, data["intensity"])
        shared_mz = cast(zarr.Array, axes["mz"])
        if (
            _read_array_int(offsets, 0) != 0
            or _read_array_int(offsets, -1) != int(intensity.shape[0])
        ):
            raise ValueError("Written offsets do not match intensity.")
        n_rows = int(offsets.shape[0]) - 1
        if n_rows != int(shared_mz.shape[0]):
            raise ValueError("Written ion row count does not match axes/mz.")
        if self.include_spectra:
            spectra = cast(zarr.Group, self._root[ZARR_GROUP_SPECTRA])
            if (
                spectra.attrs.get("row_axis") != "pixel"
                or spectra.attrs.get("encoding") != ZARR_MODE_CONTINUOUS
            ):
                raise ValueError(
                    "Written spectra/ attrs must be pixel-major continuous."
                )
            s_offsets = cast(zarr.Array, spectra["offsets"])
            s_intensity = cast(zarr.Array, spectra["intensity"])
            if (
                _read_array_int(s_offsets, 0) != 0
                or _read_array_int(s_offsets, -1)
                != int(s_intensity.shape[0])
            ):
                raise ValueError("Written spectra offsets do not match intensity.")
            if int(s_offsets.shape[0]) - 1 != int(coordinates.shape[0]):
                raise ValueError(
                    "Written spectra row count does not match coordinates."
                )
            if int(s_intensity.shape[0]) != int(coordinates.shape[0]) * int(
                shared_mz.shape[0]
            ):
                raise ValueError(
                    "Written spectra intensity does not match pixel*ion size."
                )


# ---------------------------------------------------------------------------
# Coordinate preparation
# ---------------------------------------------------------------------------


class CoordinateIndex:
    """Normalized coordinates and retained input rows."""

    def __init__(
        self,
        *,
        width: int,
        height: int,
        coordinates: np.ndarray,
        source_idx: np.ndarray,
    ) -> None:
        self.width = int(width)
        self.height = int(height)
        self.coordinates = coordinates
        self.source_idx = source_idx


def prepare_coordinate_index(
    coordinates: np.ndarray,
    *,
    width: int | None,
    height: int | None,
) -> CoordinateIndex:
    """Validate one-based coordinates and keep the last duplicate occurrence."""
    coords = _ensure_xyz_coordinates(coordinates).astype(np.int64, copy=False)
    xy = coords[:, :2]
    if width is None:
        width = int(np.max(xy[:, 0]))
    if height is None:
        height = int(np.max(xy[:, 1]))
    width = int(width)
    height = int(height)
    valid = (
        (xy[:, 0] >= 1)
        & (xy[:, 0] <= width)
        & (xy[:, 1] >= 1)
        & (xy[:, 1] <= height)
    )
    if not bool(np.all(valid)):
        raise ValueError("coordinates contains out-of-bounds values.")
    normalized = coords.copy()
    normalized[:, 2] = 1
    last: dict[tuple[int, int, int], int] = {}
    for index, row in enumerate(normalized):
        last[(int(row[0]), int(row[1]), int(row[2]))] = index
    source_idx = np.asarray(sorted(last.values()), dtype=np.intp)
    retained = normalized[source_idx].astype("<u4", copy=False)
    return CoordinateIndex(
        width=width,
        height=height,
        coordinates=retained,
        source_idx=source_idx,
    )


def prepare_pixel_coordinate_index(
    coordinates: np.ndarray,
    *,
    width: int | None,
    height: int | None,
) -> CoordinateIndex:
    """Validate pixel-major coordinates without dropping rows or rewriting z."""
    coords = _ensure_xyz_coordinates(coordinates).astype(np.int64, copy=False)
    xy = coords[:, :2]
    if width is None:
        width = int(np.max(xy[:, 0]))
    if height is None:
        height = int(np.max(xy[:, 1]))
    width = int(width)
    height = int(height)
    valid = (
        (xy[:, 0] >= 1)
        & (xy[:, 0] <= width)
        & (xy[:, 1] >= 1)
        & (xy[:, 1] <= height)
    )
    if not bool(np.all(valid)):
        raise ValueError("coordinates contains out-of-bounds values.")
    return CoordinateIndex(
        width=width,
        height=height,
        coordinates=coords.astype("<u4", copy=False),
        source_idx=np.arange(coords.shape[0], dtype=np.intp),
    )


# ---------------------------------------------------------------------------
# Encoding and dtype helpers
# ---------------------------------------------------------------------------


def _infer_encoding(
    mz_axis: np.ndarray,
    lengths: np.ndarray,
) -> ZarrStorageMode | None:
    total_points = int(np.sum(lengths, dtype=np.int64))
    shared_axis = bool(
        mz_axis.size == int(lengths[0])
        and np.all(lengths == lengths[0])
    )
    point_axis = mz_axis.size == total_points
    if shared_axis and point_axis:
        return None
    if shared_axis:
        return ZARR_MODE_CONTINUOUS
    if point_axis:
        return ZARR_MODE_PROCESSED
    raise ValueError(
        "Cannot infer Zarr encoding from flat batch: "
        f"mz_size={mz_axis.size}, total_points={total_points}, "
        f"lengths={lengths.tolist()}."
    )


def update_metadata_attrs(
    path: str | os.PathLike[str],
    items: dict,
) -> None:
    """Upsert key/value pairs into a finalized MSI Zarr store's metadata group.

    Existing keys are overwritten; new keys are inserted. Values must be
    JSON-serializable scalars accepted by Zarr attrs.
    """
    if not items:
        return
    root = zarr.open_group(_normalize_zarr_path(path), mode="r+")
    if "metadata" not in root:
        raise ValueError("Zarr store has no metadata group.")
    metadata = cast(zarr.Group, root["metadata"])
    metadata.attrs.update(items)


def _normalize_zarr_path(path: str | os.PathLike[str]) -> str:
    value = str(path)
    return value if value.endswith(".zarr") else f"{value}.zarr"


def _read_array_int(array: zarr.Array, index: int) -> int:
    """Read one integer scalar without exposing Zarr's broad selection type."""
    return int(np.asarray(array[index]).item())


def _normalize_mz_dtype(dtype: np.dtype | type) -> np.dtype:
    normalized = np.dtype(dtype)
    if normalized.kind != "f" or normalized.itemsize not in {4, 8}:
        raise ValueError("mz_dtype must be float32 or float64.")
    return normalized.newbyteorder("<")


def _ensure_xyz_coordinates(coordinates: np.ndarray) -> np.ndarray:
    coords = np.asarray(coordinates)
    if coords.ndim == 1:
        coords = coords.reshape(1, -1)
    if coords.ndim != 2 or coords.shape[1] not in {2, 3}:
        raise ValueError("coordinates must have shape (n, 2) or (n, 3).")
    if coords.shape[1] == 2:
        coords = np.column_stack(
            (coords, np.ones(coords.shape[0], dtype=coords.dtype))
        )
    return coords


# ---------------------------------------------------------------------------
# Chunking and batch helpers
# ---------------------------------------------------------------------------


def _rows_per_chunk(
    row_points: int,
    itemsize: int,
    target_chunk_bytes: int,
) -> int:
    row_bytes = max(1, int(row_points) * int(itemsize))
    return max(1, int(target_chunk_bytes) // row_bytes)


def _chunk_points(
    row_points: int | None,
    itemsize: int,
    target_chunk_bytes: int,
) -> int:
    if row_points is None:
        return max(1, int(target_chunk_bytes) // int(itemsize))
    return max(1, _rows_per_chunk(row_points, itemsize, target_chunk_bytes) * row_points)


def _split_source_indices(
    source_idx: np.ndarray,
    batches: list[BatchInfo],
) -> list[np.ndarray]:
    selections: list[np.ndarray] = []
    for batch in batches:
        start = batch.global_pixel_start
        end = start + batch.n_pixels
        selected = source_idx[(source_idx >= start) & (source_idx < end)] - start
        selections.append(selected.astype(np.intp, copy=False))
    if sum(selection.size for selection in selections) != source_idx.size:
        raise ValueError("Coordinate selection does not match temporary batches.")
    return selections
