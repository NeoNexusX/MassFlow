"""Unified reader and writer for the MassFlow MSI Zarr schema."""

from __future__ import annotations

import os
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from typing import Any, Literal

import numpy as np
import zarr

from massflow.module.ms_meta_data import MetaDataBase
from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.tools.ion_major_temp_store import BatchInfo, IonMajorTempStore
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tools.zarr_io")

ZARR_FORMAT_NAME = "massflow.msi_zarr"
ZARR_FORMAT_VERSION = "1.0"
ZARR_MODE_CONTINUOUS = "continuous"
ZARR_MODE_PROCESSED = "processed"
ZarrStorageMode = Literal["continuous", "processed"]
ZarrRowAxis = Literal["pixel", "ion"]
CoordinateBase = Literal[0, 1]


class MSIZarrReader:
    """Read one finalized MassFlow MSI Zarr store."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = _normalize_zarr_path(path)
        self.root = zarr.open_group(self.path, mode="r")
        self._validate_structure()
        self.data = self.root["data"]
        self.intensity = self.data["intensity"]
        self.offsets = self.data["offsets"]
        self.coordinates = self.root["axes/coordinates"]
        self.shared_mz_axis = (
            self.root["axes/mz"] if "mz" in self.root["axes"] else None
        )
        self.point_mz = self.data["mz"] if "mz" in self.data else None
        self.metadata = MetaDataZarr.from_zarr(self.root)

    @property
    def row_axis(self) -> str:
        return str(self.data.attrs["row_axis"])

    @property
    def encoding(self) -> str:
        return str(self.data.attrs["encoding"])

    @property
    def n_rows(self) -> int:
        return int(self.offsets.shape[0]) - 1

    @property
    def n_pixels(self) -> int:
        return int(self.coordinates.shape[0])

    @property
    def n_ions(self) -> int:
        if self.shared_mz_axis is None:
            raise ValueError("This store does not have a shared m/z axis.")
        return int(self.shared_mz_axis.shape[0])

    @property
    def spatial_shape(self) -> tuple[int, int]:
        height, width = self.root.attrs["spatial_shape"]
        return int(height), int(width)

    def read_row(self, index: int) -> np.ndarray:
        index = self._normalize_row_index(index)
        start = int(self.offsets[index])
        end = int(self.offsets[index + 1])
        return np.asarray(self.intensity[start:end])

    def read_rows(self, start: int, end: int) -> np.ndarray:
        start = max(0, int(start))
        end = min(self.n_rows, int(end))
        if end < start:
            raise ValueError("end must be greater than or equal to start.")
        point_start = int(self.offsets[start])
        point_end = int(self.offsets[end])
        values = np.asarray(self.intensity[point_start:point_end])
        if start == end:
            return values.reshape((0, 0))
        lengths = np.diff(np.asarray(self.offsets[start : end + 1], dtype=np.int64))
        if not np.all(lengths == lengths[0]):
            raise ValueError("read_rows requires uniform row lengths.")
        return values.reshape(end - start, int(lengths[0]))

    def read_spectrum(self, pixel_index: int) -> tuple[np.ndarray, np.ndarray]:
        if self.row_axis != "pixel":
            raise NotImplementedError(
                "Pixel spectrum access is not efficient for ion-major stores."
            )
        intensity = self.read_row(pixel_index)
        if self.encoding == ZARR_MODE_CONTINUOUS:
            assert self.shared_mz_axis is not None
            mz = np.asarray(self.shared_mz_axis[:], dtype=np.float64)
        else:
            assert self.point_mz is not None
            pixel_index = self._normalize_row_index(pixel_index)
            start = int(self.offsets[pixel_index])
            end = int(self.offsets[pixel_index + 1])
            mz = np.asarray(self.point_mz[start:end], dtype=np.float64)
        return mz, intensity

    def read_ion_vector(self, ion_index: int) -> np.ndarray:
        if self.row_axis != "ion":
            raise NotImplementedError(
                "Ion-vector access is not efficient for pixel-major stores."
            )
        return self.read_row(ion_index)

    def read_ion_image(self, ion_index: int) -> np.ndarray:
        values = self.read_ion_vector(ion_index)
        height, width = self.spatial_shape
        image = np.zeros((height, width), dtype=self.intensity.dtype)
        coordinates = np.asarray(self.coordinates[:], dtype=np.int64)
        image[coordinates[:, 1], coordinates[:, 0]] = values
        return image

    def validate_full(self) -> None:
        offsets = np.asarray(self.offsets[:], dtype=np.int64)
        if np.any(np.diff(offsets) < 0):
            raise ValueError("data/offsets must be monotonic.")
        coordinates = np.asarray(self.coordinates[:], dtype=np.int64)
        height, width = self.spatial_shape
        valid = (
            (coordinates[:, 0] >= 0)
            & (coordinates[:, 0] < width)
            & (coordinates[:, 1] >= 0)
            & (coordinates[:, 1] < height)
        )
        if not bool(np.all(valid)):
            raise ValueError("axes/coordinates contains out-of-bounds values.")
        if np.unique(coordinates, axis=0).shape[0] != coordinates.shape[0]:
            raise ValueError("axes/coordinates contains duplicates.")
        if self.row_axis == "ion":
            lengths = np.diff(offsets)
            if np.any(lengths != self.n_pixels):
                raise ValueError("Every ion row must contain n_pixels values.")
        elif self.encoding == ZARR_MODE_CONTINUOUS:
            assert self.shared_mz_axis is not None
            lengths = np.diff(offsets)
            if np.any(lengths != self.shared_mz_axis.shape[0]):
                raise ValueError("Every continuous spectrum must match axes/mz.")
        if "stats" in self.root:
            stats = self.root["stats"]
            if "mean_spectrum" in stats:
                if self.shared_mz_axis is None:
                    raise ValueError(
                        "stats/mean_spectrum requires a shared m/z axis."
                    )
                if stats["mean_spectrum"].ndim != 1 or int(
                    stats["mean_spectrum"].shape[0]
                ) != int(self.shared_mz_axis.shape[0]):
                    raise ValueError(
                        "stats/mean_spectrum must match axes/mz length."
                    )
            if "tic" in stats and (
                stats["tic"].ndim != 1
                or int(stats["tic"].shape[0]) != self.n_pixels
            ):
                raise ValueError("stats/tic must match the coordinate count.")

    def close(self) -> None:
        self.root = None

    def _normalize_row_index(self, index: int) -> int:
        index = int(index)
        if index < 0:
            index += self.n_rows
        if index < 0 or index >= self.n_rows:
            raise IndexError("row index out of range")
        return index

    def _validate_structure(self) -> None:
        attrs = self.root.attrs
        required_attrs = {
            "format",
            "format_version",
            "write_state",
            "coordinate_order",
            "coordinate_base",
            "spatial_shape",
        }
        missing_attrs = required_attrs - set(attrs)
        if missing_attrs:
            raise ValueError(
                f"Missing required root attrs: {sorted(missing_attrs)}"
            )
        if attrs.get("format") != ZARR_FORMAT_NAME:
            raise ValueError("Not a MassFlow MSI Zarr store.")
        if attrs.get("format_version") != ZARR_FORMAT_VERSION:
            raise ValueError(
                f"Unsupported MassFlow MSI Zarr version: {attrs.get('format_version')!r}"
            )
        if attrs.get("write_state") != "complete":
            raise ValueError("Zarr store is not finalized.")
        if list(attrs["coordinate_order"]) != ["x", "y", "z"]:
            raise ValueError("coordinate_order must be ['x', 'y', 'z'].")
        if attrs["coordinate_base"] != 0:
            raise ValueError("coordinate_base must be 0.")
        spatial_shape = attrs["spatial_shape"]
        if (
            not isinstance(spatial_shape, (list, tuple))
            or len(spatial_shape) != 2
            or any(
                isinstance(value, bool)
                or not isinstance(value, (int, np.integer))
                or int(value) <= 0
                for value in spatial_shape
            )
        ):
            raise ValueError("spatial_shape must be [height, width] with positive integers.")
        if "data" not in self.root or "axes" not in self.root:
            raise ValueError("Zarr store requires data/ and axes/.")
        if "metadata" not in self.root:
            raise ValueError("Zarr store requires metadata/.")
        data = self.root["data"]
        axes = self.root["axes"]
        if "intensity" not in data or "offsets" not in data:
            raise ValueError("data/intensity and data/offsets are required.")
        if "coordinates" not in axes:
            raise ValueError("axes/coordinates is required.")
        row_axis = data.attrs.get("row_axis")
        encoding = data.attrs.get("encoding")
        if row_axis not in {"pixel", "ion"}:
            raise ValueError(f"Unsupported row_axis: {row_axis!r}")
        if encoding not in {ZARR_MODE_CONTINUOUS, ZARR_MODE_PROCESSED}:
            raise ValueError(f"Unsupported encoding: {encoding!r}")
        if row_axis == "ion" and encoding != ZARR_MODE_CONTINUOUS:
            raise ValueError("ion-major stores must use continuous encoding.")
        if encoding == ZARR_MODE_CONTINUOUS:
            if "mz" not in axes or "mz" in data:
                raise ValueError(
                    "continuous stores require axes/mz and forbid data/mz."
                )
        elif "mz" not in data or "mz" in axes:
            raise ValueError("processed stores require data/mz and forbid axes/mz.")
        intensity = data["intensity"]
        offsets = data["offsets"]
        coordinates = axes["coordinates"]
        if offsets.ndim != 1 or intensity.ndim != 1:
            raise ValueError("data/intensity and data/offsets must be 1-D.")
        if not np.can_cast(offsets.dtype, np.dtype(np.int64), casting="safe"):
            raise ValueError("data/offsets must use an integer dtype safe for int64.")
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("axes/coordinates must have shape (n_pixels, 3).")
        if np.dtype(coordinates.dtype) != np.dtype("<u4"):
            raise ValueError("axes/coordinates must use <u4.")
        if encoding == ZARR_MODE_CONTINUOUS:
            shared_mz = axes["mz"]
            if shared_mz.ndim != 1 or np.dtype(shared_mz.dtype) != np.dtype("<f8"):
                raise ValueError("axes/mz must be a 1-D <f8 array.")
            if shared_mz.attrs.get("unit") != "m/z":
                raise ValueError("axes/mz attrs.unit must be 'm/z'.")
        else:
            point_mz = data["mz"]
            if point_mz.ndim != 1 or np.dtype(point_mz.dtype) != np.dtype("<f8"):
                raise ValueError("data/mz must be a 1-D <f8 array.")
            if point_mz.attrs.get("unit") != "m/z":
                raise ValueError("data/mz attrs.unit must be 'm/z'.")
        n_rows = int(offsets.shape[0]) - 1
        if n_rows < 0:
            raise ValueError("data/offsets must contain at least one value.")
        if int(offsets[0]) != 0 or int(offsets[-1]) != int(intensity.shape[0]):
            raise ValueError("data/offsets boundaries do not match intensity.")
        if row_axis == "pixel" and n_rows != int(coordinates.shape[0]):
            raise ValueError("pixel-major row count must equal coordinate count.")
        if row_axis == "ion" and n_rows != int(axes["mz"].shape[0]):
            raise ValueError("ion-major row count must equal axes/mz length.")
        if encoding == ZARR_MODE_PROCESSED and (
            int(data["mz"].shape[0]) != int(intensity.shape[0])
        ):
            raise ValueError("data/mz length must equal data/intensity length.")
        metadata_attrs = self.root["metadata"].attrs
        continuous = metadata_attrs.get("continuous")
        processed = metadata_attrs.get("processed")
        if encoding == ZARR_MODE_CONTINUOUS and (
            continuous is not True or processed is True
        ):
            raise ValueError(
                "metadata continuous/processed conflicts with continuous encoding."
            )
        if encoding == ZARR_MODE_PROCESSED and (
            processed is not True or continuous is True
        ):
            raise ValueError(
                "metadata continuous/processed conflicts with processed encoding."
            )


class MSIZarrWriter:
    """Write pixel-major or ion-major MassFlow MSI Zarr stores."""

    def __init__(
        self,
        path: str | os.PathLike[str],
        *,
        row_axis: ZarrRowAxis,
        encoding: ZarrStorageMode,
        metadata: MetaDataBase | None,
        width: int | None = None,
        height: int | None = None,
        intensity_dtype: np.dtype | type = np.float32,
        target_chunk_bytes: int = 4 * 1024 * 1024,
        add_buffer_size: int = 256,
        temp_dir: str | os.PathLike[str] | None = None,
        temp_workers: int = 4,
        ion_block: int = 256,
        compute_mean_spectrum: bool = False,
        compute_tic: bool = False,
        coordinate_base: CoordinateBase | None = None,
    ) -> None:
        if row_axis not in {"pixel", "ion"}:
            raise ValueError("row_axis must be 'pixel' or 'ion'.")
        if encoding not in {ZARR_MODE_CONTINUOUS, ZARR_MODE_PROCESSED}:
            raise ValueError("encoding must be 'continuous' or 'processed'.")
        if row_axis == "ion" and encoding != ZARR_MODE_CONTINUOUS:
            raise ValueError("ion-major output requires continuous encoding.")
        self.path = _normalize_zarr_path(path)
        self.row_axis = row_axis
        self.encoding = encoding
        self.meta = MetaDataZarr.from_metadata(metadata)
        self.width = None if width is None else int(width)
        self.height = None if height is None else int(height)
        self.intensity_dtype = np.dtype(intensity_dtype)
        self.target_chunk_bytes = max(1, int(target_chunk_bytes))
        self.add_buffer_size = max(1, int(add_buffer_size))
        self.ion_block = max(1, int(ion_block))
        self.compute_mean_spectrum = bool(compute_mean_spectrum)
        self.compute_tic = bool(compute_tic)
        if coordinate_base not in {None, 0, 1}:
            raise ValueError("coordinate_base must be 0, 1, or None.")
        self.coordinate_base = coordinate_base
        self._temp_dir = temp_dir
        self._temp_workers = max(1, int(temp_workers))
        self._root: Any | None = None
        self._data: Any | None = None
        self._intensity: Any | None = None
        self._offsets: Any | None = None
        self._point_mz: Any | None = None
        self._shared_mz: np.ndarray | None = None
        self._coordinates: list[np.ndarray] = []
        self._temp_store: IonMajorTempStore | None = None
        self._total_points = 0
        self._n_input_pixels = 0
        self._finalized = False
        self._add_intensity: list[np.ndarray] = []
        self._add_mz: list[np.ndarray] = []
        self._add_coordinates: list[np.ndarray] = []

    def add_spectrum(
        self,
        intensity: np.ndarray,
        coordinates: np.ndarray | list[int] | tuple[int, ...],
        mz: np.ndarray,
    ) -> None:
        intensity = np.asarray(intensity, dtype=self.intensity_dtype)
        mz = np.asarray(mz, dtype=np.float64)
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
        mz_axis: np.ndarray,
    ) -> None:
        if self._finalized:
            raise ValueError("Cannot write after finalize().")
        lengths = np.asarray(lengths, dtype=np.int64)
        coordinates = _ensure_xyz_coordinates(coordinates)
        intensity = np.asarray(intensity_flat, dtype=self.intensity_dtype)
        mz = np.asarray(mz_axis, dtype=np.float64)
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
        if self._n_input_pixels == 0:
            raise ValueError("No spectra were written.")
        coordinates = np.vstack(self._coordinates)
        coordinate_index = prepare_coordinate_index(
            coordinates,
            width=self.width,
            height=self.height,
            coordinate_base=self.coordinate_base,
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

    def _write_ion_batch(
        self,
        intensity: np.ndarray,
        lengths: np.ndarray,
    ) -> None:
        assert self._shared_mz is not None
        if self._temp_store is None:
            self._temp_store = IonMajorTempStore(
                temp_dir=self._temp_dir,
                dtype=self.intensity_dtype,
                workers=self._temp_workers,
            )
        self._temp_store.append_batch(
            intensity,
            n_pixels=int(lengths.size),
            n_ions=int(self._shared_mz.size),
        )

    def _ensure_pixel_arrays(self, lengths: np.ndarray) -> None:
        if self._root is not None:
            return
        Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self._root = zarr.open_group(self.path, mode="w")
        self._write_root_attrs(state="writing")
        axes = self._root.create_group("axes")
        self._data = self._root.create_group("data")
        self._data.attrs.update(
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
        self._intensity = self._data.create_array(
            "intensity",
            shape=(0,),
            chunks=(chunk_points,),
            dtype=self.intensity_dtype,
        )
        self._offsets = self._data.create_array(
            "offsets",
            data=np.array([0], dtype="<i8"),
            chunks=(max(2, min(4096, int(lengths.size) + 1)),),
        )
        if self.encoding == ZARR_MODE_CONTINUOUS:
            assert self._shared_mz is not None
            mz_array = axes.create_array("mz", data=self._shared_mz)
            mz_array.attrs["unit"] = "m/z"
        else:
            self._point_mz = self._data.create_array(
                "mz",
                shape=(0,),
                chunks=(chunk_points,),
                dtype="<f8",
            )
            self._point_mz.attrs["unit"] = "m/z"

    def _finalize_pixel(self, coordinate_index: "CoordinateIndex") -> None:
        assert self._root is not None
        if coordinate_index.source_idx.size != self._n_input_pixels:
            self._compact_pixel_arrays(coordinate_index.source_idx)
        axes = self._root["axes"]
        axes.create_array(
            "coordinates",
            data=coordinate_index.coordinates.astype("<u4", copy=False),
            overwrite=True,
        )

    def _compact_pixel_arrays(self, source_idx: np.ndarray) -> None:
        assert self._intensity is not None and self._offsets is not None
        old_offsets = np.asarray(self._offsets[:], dtype=np.int64)
        destination = 0
        new_offsets = np.empty(source_idx.size + 1, dtype="<i8")
        new_offsets[0] = 0
        for output_row, source_row in enumerate(source_idx):
            start = int(old_offsets[source_row])
            end = int(old_offsets[source_row + 1])
            intensity = np.asarray(self._intensity[start:end])
            next_destination = destination + intensity.size
            self._intensity[destination:next_destination] = intensity
            if self._point_mz is not None:
                mz = np.asarray(self._point_mz[start:end])
                self._point_mz[destination:next_destination] = mz
            destination = next_destination
            new_offsets[output_row + 1] = destination
        self._intensity.resize((destination,))
        if self._point_mz is not None:
            self._point_mz.resize((destination,))
        self._offsets.resize(new_offsets.shape)
        self._offsets[:] = new_offsets
        self._total_points = destination

    def _finalize_ion(self, coordinate_index: "CoordinateIndex") -> None:
        if self._temp_store is None or self._shared_mz is None:
            raise ValueError("Ion-major writer has no temporary data.")
        batches = self._temp_store.finish()
        Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self._root = zarr.open_group(self.path, mode="w")
        self._write_root_attrs(state="writing")
        axes = self._root.create_group("axes")
        axes.create_array(
            "coordinates",
            data=coordinate_index.coordinates.astype("<u4", copy=False),
        )
        mz_array = axes.create_array("mz", data=self._shared_mz)
        mz_array.attrs["unit"] = "m/z"
        self._data = self._root.create_group("data")
        self._data.attrs.update(
            {"row_axis": self.row_axis, "encoding": self.encoding}
        )
        n_pixels = int(coordinate_index.coordinates.shape[0])
        n_ions = int(self._shared_mz.size)
        rows_per_chunk = _rows_per_chunk(
            n_pixels,
            self.intensity_dtype.itemsize,
            self.target_chunk_bytes,
        )
        chunk_points = max(1, rows_per_chunk * n_pixels)
        self._intensity = self._data.create_array(
            "intensity",
            shape=(n_ions * n_pixels,),
            chunks=(chunk_points,),
            dtype=self.intensity_dtype,
        )
        offsets = np.arange(n_ions + 1, dtype=np.int64) * n_pixels
        self._offsets = self._data.create_array("offsets", data=offsets.astype("<i8"))
        selections = _split_source_indices(coordinate_index.source_idx, batches)
        merge_rows = max(rows_per_chunk, self.ion_block)
        merge_rows = ((merge_rows + rows_per_chunk - 1) // rows_per_chunk) * rows_per_chunk
        mean_spectrum = (
            np.empty(n_ions, dtype=np.float64)
            if self.compute_mean_spectrum
            else None
        )
        tic = np.zeros(n_pixels, dtype=np.float64) if self.compute_tic else None

        def assemble_band(ion_start: int, ion_end: int) -> np.ndarray:
            band = np.empty((ion_end - ion_start, n_pixels), dtype=self.intensity_dtype)
            output_start = 0
            for batch, local_indices in zip(batches, selections):
                source = self._temp_store.read_ion_block(batch, ion_start, ion_end)
                output_end = output_start + local_indices.size
                band[:, output_start:output_end] = source[:, local_indices]
                output_start = output_end
            return band

        starts = list(range(0, n_ions, merge_rows))
        with ThreadPoolExecutor(max_workers=1) as executor:
            future: Future[np.ndarray] | None = None
            for index, ion_start in enumerate(starts):
                ion_end = min(n_ions, ion_start + merge_rows)
                band = (
                    assemble_band(ion_start, ion_end)
                    if future is None
                    else future.result()
                )
                future = None
                if index + 1 < len(starts):
                    next_start = starts[index + 1]
                    next_end = min(n_ions, next_start + merge_rows)
                    future = executor.submit(assemble_band, next_start, next_end)
                point_start = ion_start * n_pixels
                point_end = ion_end * n_pixels
                self._intensity[point_start:point_end] = band.reshape(-1)
                if mean_spectrum is not None:
                    mean_spectrum[ion_start:ion_end] = np.mean(
                        band, axis=1, dtype=np.float64
                    )
                if tic is not None:
                    tic += np.sum(band, axis=0, dtype=np.float64)
        self._write_stats(mean_spectrum=mean_spectrum, tic=tic)

    def _write_stats(
        self,
        *,
        mean_spectrum: np.ndarray | None,
        tic: np.ndarray | None,
    ) -> None:
        if mean_spectrum is None and tic is None:
            return
        assert self._root is not None
        stats = self._root.create_group("stats")
        if mean_spectrum is not None:
            array = stats.create_array("mean_spectrum", data=mean_spectrum)
            array.attrs.update(
                {
                    "reduced_axis": "pixel",
                    "denominator": "sampled_pixels",
                    "source": "/data",
                }
            )
        if tic is not None:
            array = stats.create_array("tic", data=tic)
            array.attrs.update({"reduced_axis": "mz", "source": "/data"})

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
            self._shared_mz = mz.astype("<f8", copy=True)
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
        if self.encoding == ZARR_MODE_CONTINUOUS:
            mz = self._add_mz[0]
            if any(not np.array_equal(row, mz) for row in self._add_mz[1:]):
                raise ValueError("Continuous add_spectrum calls must share m/z.")
        else:
            mz = np.concatenate(self._add_mz)
        self._add_intensity.clear()
        self._add_mz.clear()
        self._add_coordinates.clear()
        self.write_flat_batch(intensity, lengths, coordinates, mz)

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
                "coordinate_base": 0,
                "spatial_shape": spatial_shape,
            }
        )

    def _validate_written_structure(self) -> None:
        assert self._root is not None
        data = self._root["data"]
        axes = self._root["axes"]
        offsets = data["offsets"]
        intensity = data["intensity"]
        coordinates = axes["coordinates"]
        if int(offsets[0]) != 0 or int(offsets[-1]) != int(intensity.shape[0]):
            raise ValueError("Written offsets do not match intensity.")
        n_rows = int(offsets.shape[0]) - 1
        if self.row_axis == "pixel" and n_rows != int(coordinates.shape[0]):
            raise ValueError("Written pixel row count does not match coordinates.")
        if self.row_axis == "ion" and n_rows != int(axes["mz"].shape[0]):
            raise ValueError("Written ion row count does not match axes/mz.")
        if self.encoding == ZARR_MODE_PROCESSED and (
            int(data["mz"].shape[0]) != int(intensity.shape[0])
        ):
            raise ValueError("Written data/mz does not match intensity.")


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
    coordinate_base: CoordinateBase | None = None,
) -> CoordinateIndex:
    """Normalize coordinates and keep the last duplicate occurrence."""
    coords = _ensure_xyz_coordinates(coordinates).astype(np.int64, copy=False)
    xy = coords[:, :2]
    if coordinate_base is None:
        coordinate_base = 0 if bool(np.any(xy == 0)) else 1
    xy0 = xy.copy() if coordinate_base == 0 else xy - 1
    z = np.zeros(coords.shape[0], dtype=np.int64)
    if width is None:
        width = int(np.max(xy0[:, 0])) + 1
    if height is None:
        height = int(np.max(xy0[:, 1])) + 1
    width = int(width)
    height = int(height)
    valid = (
        (xy0[:, 0] >= 0)
        & (xy0[:, 0] < width)
        & (xy0[:, 1] >= 0)
        & (xy0[:, 1] < height)
    )
    if not bool(np.all(valid)):
        raise ValueError("coordinates contains out-of-bounds values.")
    normalized = np.column_stack((xy0, z))
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


def storage_mode_from_metadata(metadata: MetaDataBase | None) -> ZarrStorageMode:
    if metadata is not None and metadata.continuous is True:
        return ZARR_MODE_CONTINUOUS
    if metadata is not None and metadata.processed is True:
        return ZARR_MODE_PROCESSED
    raise ValueError("Metadata must declare continuous=True or processed=True.")


def _normalize_zarr_path(path: str | os.PathLike[str]) -> str:
    value = str(path)
    return value if value.endswith(".zarr") else f"{value}.zarr"


def _ensure_xyz_coordinates(coordinates: np.ndarray) -> np.ndarray:
    coords = np.asarray(coordinates)
    if coords.ndim == 1:
        coords = coords.reshape(1, -1)
    if coords.ndim != 2 or coords.shape[1] not in {2, 3}:
        raise ValueError("coordinates must have shape (n, 2) or (n, 3).")
    if coords.shape[1] == 2:
        coords = np.column_stack(
            (coords, np.zeros(coords.shape[0], dtype=coords.dtype))
        )
    return coords


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
