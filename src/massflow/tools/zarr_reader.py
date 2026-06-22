"""Reader for the MassFlow MSI Zarr schema."""

from __future__ import annotations

import os
from typing import Any, cast

import numpy as np
import zarr

from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.tools.zarr_writer import (
    ZARR_FORMAT_NAME,
    ZARR_FORMAT_VERSION,
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
    ZarrRowAxis,
    ZarrStorageMode,
)


class MSIZarrReader:
    """Read one finalized MassFlow MSI Zarr store."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = _normalize_zarr_path(path)
        self.root: zarr.Group = zarr.open_group(self.path, mode="r")
        self._validate_structure()
        self.data = cast(zarr.Group, self.root["data"])
        axes = cast(zarr.Group, self.root["axes"])
        self.intensity = cast(zarr.Array, self.data["intensity"])
        self.offsets = cast(zarr.Array, self.data["offsets"])
        self.coordinates = cast(zarr.Array, axes["coordinates"])
        self.shared_mz_axis = (
            cast(zarr.Array, axes["mz"]) if "mz" in axes else None
        )
        self.point_mz = (
            cast(zarr.Array, self.data["mz"]) if "mz" in self.data else None
        )
        self.metadata = MetaDataZarr.from_zarr(self.root)

    # --- Store metadata and dimensions ---------------------------------------

    @property
    def row_axis(self) -> ZarrRowAxis:
        return cast(ZarrRowAxis, str(self.data.attrs["row_axis"]))

    @property
    def encoding(self) -> ZarrStorageMode:
        return cast(ZarrStorageMode, str(self.data.attrs["encoding"]))

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
        spatial_shape = self.root.attrs["spatial_shape"]
        if not isinstance(spatial_shape, (list, tuple)) or len(spatial_shape) != 2:
            raise ValueError("spatial_shape must contain height and width.")
        height, width = cast(tuple[int, int], tuple(spatial_shape))
        return int(height), int(width)

    @property
    def intensity_dtype(self) -> np.dtype:
        return np.dtype(self.intensity.dtype)

    @property
    def mz_dtype(self) -> np.dtype:
        array = self.shared_mz_axis if self.shared_mz_axis is not None else self.point_mz
        if array is None:
            raise ValueError("This store does not contain m/z values.")
        return np.dtype(array.dtype)

    def read_coordinates(self, indices: np.ndarray | None = None) -> np.ndarray:
        if indices is None:
            return np.asarray(self.coordinates[:])
        indices = np.asarray(indices, dtype=np.int64)
        return np.asarray(
            self.coordinates.get_orthogonal_selection((indices, slice(None)))
        )

    def read_shared_mz_axis(self) -> np.ndarray:
        if self.shared_mz_axis is None:
            raise ValueError("This store does not have a shared m/z axis.")
        return np.asarray(self.shared_mz_axis[:])

    def read_mean_spectrum(self) -> np.ndarray | None:
        if "stats" not in self.root:
            return None
        stats = cast(zarr.Group, self.root["stats"])
        if "mean_spectrum" not in stats:
            return None
        array = cast(zarr.Array, stats["mean_spectrum"])
        return np.asarray(array[:])

    # --- Pixel-major batch access --------------------------------------------

    def read_flat_batch(
        self,
        indices: np.ndarray,
        *,
        include_mz: bool = True,
    ) -> tuple[np.ndarray | None, np.ndarray, np.ndarray, np.ndarray]:
        if self.row_axis != "pixel":
            raise ValueError("read_flat_batch requires a pixel-major store.")
        indices = np.asarray(indices, dtype=np.int64)
        if indices.ndim != 1:
            raise ValueError("indices must be one-dimensional.")
        coordinates = self.read_coordinates(indices)
        if indices.size == 0:
            empty = np.array([], dtype=self.intensity_dtype)
            return None, empty, np.array([], dtype="<u4"), coordinates

        starts = self._read_offsets(indices)
        ends = self._read_offsets(indices + 1)
        lengths = (ends - starts).astype("<u4")
        intensity = self._read_indexed_values(self.intensity, indices, starts, ends)
        if self.encoding == ZARR_MODE_CONTINUOUS:
            mz = self.read_shared_mz_axis()
        elif include_mz:
            assert self.point_mz is not None
            mz = self._read_indexed_values(self.point_mz, indices, starts, ends)
        else:
            mz = None
        return mz, intensity, lengths, coordinates

    # --- Row and spectrum access ---------------------------------------------

    def read_row(self, index: int) -> np.ndarray:
        index = self._normalize_row_index(index)
        start = _read_array_int(self.offsets, index)
        end = _read_array_int(self.offsets, index + 1)
        return np.asarray(self.intensity[start:end])

    def read_rows(self, start: int, end: int) -> np.ndarray:
        start = max(0, int(start))
        end = min(self.n_rows, int(end))
        if end < start:
            raise ValueError("end must be greater than or equal to start.")
        point_start = _read_array_int(self.offsets, start)
        point_end = _read_array_int(self.offsets, end)
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
            mz = np.asarray(self.shared_mz_axis[:])
        else:
            assert self.point_mz is not None
            pixel_index = self._normalize_row_index(pixel_index)
            start = _read_array_int(self.offsets, pixel_index)
            end = _read_array_int(self.offsets, pixel_index + 1)
            mz = np.asarray(self.point_mz[start:end])
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
        image[coordinates[:, 1] - 1, coordinates[:, 0] - 1] = values
        return image

    # --- Validation and lifecycle --------------------------------------------

    def validate_full(self) -> None:
        offsets = np.asarray(self.offsets[:], dtype=np.int64)
        if np.any(np.diff(offsets) < 0):
            raise ValueError("data/offsets must be monotonic.")
        coordinates = np.asarray(self.coordinates[:], dtype=np.int64)
        height, width = self.spatial_shape
        valid = (
            (coordinates[:, 0] >= 1)
            & (coordinates[:, 0] <= width)
            & (coordinates[:, 1] >= 1)
            & (coordinates[:, 1] <= height)
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
            stats = cast(zarr.Group, self.root["stats"])
            if "mean_spectrum" in stats:
                mean_spectrum = cast(zarr.Array, stats["mean_spectrum"])
                if self.shared_mz_axis is None:
                    raise ValueError(
                        "stats/mean_spectrum requires a shared m/z axis."
                    )
                if mean_spectrum.ndim != 1 or int(
                    mean_spectrum.shape[0]
                ) != int(self.shared_mz_axis.shape[0]):
                    raise ValueError(
                        "stats/mean_spectrum must match axes/mz length."
                    )

    def close(self) -> None:
        """Release reader resources.

        Zarr groups do not hold a closeable file handle, so no explicit action
        is required.
        """

    # --- Internal selection helpers ------------------------------------------

    def _read_offsets(self, indices: np.ndarray) -> np.ndarray:
        if indices.size > 1 and bool(np.all(np.diff(indices) == 1)):
            start = int(indices[0])
            end = int(indices[-1]) + 1
            return np.asarray(self.offsets[start:end], dtype=np.int64)
        return np.asarray(
            self.offsets.get_orthogonal_selection((indices,)),
            dtype=np.int64,
        )

    def _read_indexed_values(
        self,
        array: Any,
        indices: np.ndarray,
        starts: np.ndarray,
        ends: np.ndarray,
    ) -> np.ndarray:
        if indices.size > 1 and bool(np.all(np.diff(indices) == 1)):
            return np.asarray(array[int(starts[0]) : int(ends[-1])])
        return np.concatenate(
            [
                np.asarray(array[int(start) : int(end)])
                for start, end in zip(starts, ends)
            ]
        )

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
        coordinate_order = attrs["coordinate_order"]
        if (
            not isinstance(coordinate_order, (list, tuple))
            or list(coordinate_order) != ["x", "y", "z"]
        ):
            raise ValueError("coordinate_order must be ['x', 'y', 'z'].")
        if attrs["coordinate_base"] != 1:
            raise ValueError("coordinate_base must be 1.")
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
        data = cast(zarr.Group, self.root["data"])
        axes = cast(zarr.Group, self.root["axes"])
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
        intensity = cast(zarr.Array, data["intensity"])
        offsets = cast(zarr.Array, data["offsets"])
        coordinates = cast(zarr.Array, axes["coordinates"])
        if offsets.ndim != 1 or intensity.ndim != 1:
            raise ValueError("data/intensity and data/offsets must be 1-D.")
        if not np.can_cast(offsets.dtype, np.dtype(np.int64), casting="safe"):
            raise ValueError("data/offsets must use an integer dtype safe for int64.")
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("axes/coordinates must have shape (n_pixels, 3).")
        if np.dtype(coordinates.dtype) != np.dtype("<u4"):
            raise ValueError("axes/coordinates must use <u4.")
        if encoding == ZARR_MODE_CONTINUOUS:
            shared_mz = cast(zarr.Array, axes["mz"])
            if shared_mz.ndim != 1 or not _is_supported_mz_dtype(shared_mz.dtype):
                raise ValueError("axes/mz must be a 1-D <f4 or <f8 array.")
            if shared_mz.attrs.get("unit") != "m/z":
                raise ValueError("axes/mz attrs.unit must be 'm/z'.")
        else:
            point_mz = cast(zarr.Array, data["mz"])
            if point_mz.ndim != 1 or not _is_supported_mz_dtype(point_mz.dtype):
                raise ValueError("data/mz must be a 1-D <f4 or <f8 array.")
            if point_mz.attrs.get("unit") != "m/z":
                raise ValueError("data/mz attrs.unit must be 'm/z'.")
        n_rows = int(offsets.shape[0]) - 1
        if n_rows < 0:
            raise ValueError("data/offsets must contain at least one value.")
        if (
            _read_array_int(offsets, 0) != 0
            or _read_array_int(offsets, -1) != int(intensity.shape[0])
        ):
            raise ValueError("data/offsets boundaries do not match intensity.")
        if row_axis == "pixel" and n_rows != int(coordinates.shape[0]):
            raise ValueError("pixel-major row count must equal coordinate count.")
        if row_axis == "ion":
            shared_mz = cast(zarr.Array, axes["mz"])
            if n_rows != int(shared_mz.shape[0]):
                raise ValueError("ion-major row count must equal axes/mz length.")
        if encoding == ZARR_MODE_PROCESSED:
            point_mz = cast(zarr.Array, data["mz"])
            if int(point_mz.shape[0]) != int(intensity.shape[0]):
                raise ValueError("data/mz length must equal data/intensity length.")
        metadata = cast(zarr.Group, self.root["metadata"])
        metadata_attrs = metadata.attrs
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


# ---------------------------------------------------------------------------
# Module helpers
# ---------------------------------------------------------------------------


def _normalize_zarr_path(path: str | os.PathLike[str]) -> str:
    value = str(path)
    return value if value.endswith(".zarr") else f"{value}.zarr"


def _read_array_int(array: zarr.Array, index: int) -> int:
    """Read one integer scalar without exposing Zarr's broad selection type."""
    return int(np.asarray(array[index]).item())


def _is_supported_mz_dtype(dtype: np.dtype | type) -> bool:
    return np.dtype(dtype) in {np.dtype("<f4"), np.dtype("<f8")}
