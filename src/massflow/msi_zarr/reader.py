"""Reader for the MassFlow MSI Zarr schema."""

from __future__ import annotations

import os
from typing import Any, cast

import numpy as np
import zarr

from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.msi_zarr.schema import (
    ZARR_FORMAT_NAME,
    ZARR_FORMAT_VERSION,
    ZARR_GROUP_ION_IMAGE,
    ZARR_GROUP_SPECTRA,
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
    ZarrRowAxis,
    ZarrStorageMode,
)
from massflow.tools.logger import get_logger

logger = get_logger("massflow.msi_zarr.zarr_reader")


class MSIZarrReader:
    """Read one finalized MassFlow MSI Zarr store."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = _normalize_zarr_path(path)
        self.root: zarr.Group = zarr.open_group(self.path, mode="r")
        self._combined = False
        self._spectra_group: zarr.Group | None = None
        self._validate_structure()

    # --- Store metadata and dimensions ---------------------------------------

    @property
    def row_axis(self) -> ZarrRowAxis:
        return cast(ZarrRowAxis, str(self.data.attrs["row_axis"]))

    @property
    def encoding(self) -> ZarrStorageMode:
        return cast(ZarrStorageMode, str(self.data.attrs["encoding"]))

    @property
    def has_pixel_spectra(self) -> bool:
        """Return whether pixel-major spectra are available for batch access."""
        return self.row_axis == "pixel" or self._spectra_group is not None

    @property
    def spectra_encoding(self) -> ZarrStorageMode:
        """Return the encoding used by the pixel-major spectra view."""
        if self.row_axis == "pixel":
            return self.encoding
        if self._spectra_group is None:
            raise ValueError("This store does not contain pixel-major spectra.")
        return cast(
            ZarrStorageMode,
            str(self._spectra_group.attrs["encoding"]),
        )

    @property
    def _main_group_name(self) -> str:
        return (
            ZARR_GROUP_ION_IMAGE
            if self.row_axis == "ion"
            else ZARR_GROUP_SPECTRA
        )

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

    def read_tic(self) -> np.ndarray | None:
        if "stats" not in self.root:
            return None
        stats = cast(zarr.Group, self.root["stats"])
        if "tic" not in stats:
            return None
        array = cast(zarr.Array, stats["tic"])
        return np.asarray(array[:])

    def read_rms(self) -> np.ndarray | None:
        if "stats" not in self.root:
            return None
        stats = cast(zarr.Group, self.root["stats"])
        if "rms" not in stats:
            return None
        array = cast(zarr.Array, stats["rms"])
        return np.asarray(array[:])

    # --- Pixel-major batch access --------------------------------------------

    def read_flat_batch(
        self,
        indices: np.ndarray,
        *,
        include_mz: bool = True,
    ) -> tuple[np.ndarray | None, np.ndarray, np.ndarray, np.ndarray]:
        intensity_array, offsets_array, encoding, point_mz = (
            self._pixel_spectra_arrays()
        )
        indices = np.asarray(indices, dtype=np.int64)
        if indices.ndim != 1:
            raise ValueError("indices must be one-dimensional.")
        coordinates = self.read_coordinates(indices)
        if indices.size == 0:
            empty = np.array([], dtype=intensity_array.dtype)
            return None, empty, np.array([], dtype="<u4"), coordinates

        starts = self._read_offsets(offsets_array, indices)
        ends = self._read_offsets(offsets_array, indices + 1)
        lengths = (ends - starts).astype("<u4")
        intensity = self._read_indexed_values(
            intensity_array,
            indices,
            starts,
            ends,
        )
        if encoding == ZARR_MODE_CONTINUOUS:
            mz = self.read_shared_mz_axis()
        elif include_mz:
            assert point_mz is not None
            mz = self._read_indexed_values(point_mz, indices, starts, ends)
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
        intensity_array, offsets_array, encoding, point_mz = (
            self._pixel_spectra_arrays()
        )
        pixel_index = int(pixel_index)
        if pixel_index < 0:
            pixel_index += self.n_pixels
        if pixel_index < 0 or pixel_index >= self.n_pixels:
            raise IndexError("pixel index out of range")
        start = _read_array_int(offsets_array, pixel_index)
        end = _read_array_int(offsets_array, pixel_index + 1)
        intensity = np.asarray(intensity_array[start:end])
        if encoding == ZARR_MODE_CONTINUOUS:
            assert self.shared_mz_axis is not None
            mz = np.asarray(self.shared_mz_axis[:])
        else:
            assert point_mz is not None
            mz = np.asarray(point_mz[start:end])
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
            raise ValueError(f"{self._main_group_name}/offsets must be monotonic.")
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
        if (
            self.row_axis == "ion"
            and np.unique(coordinates, axis=0).shape[0] != coordinates.shape[0]
        ):
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
        if self._combined:
            assert self._spectra_group is not None
            spectra_offsets = cast(zarr.Array, self._spectra_group["offsets"])
            s_offsets = np.asarray(spectra_offsets[:], dtype=np.int64)
            if np.any(np.diff(s_offsets) < 0):
                raise ValueError("spectra/offsets must be monotonic.")
            if int(spectra_offsets.shape[0]) - 1 != int(
                self.coordinates.shape[0]
            ):
                raise ValueError(
                    "spectra row count must match axes/coordinates count."
                )
            if np.any(np.diff(s_offsets) != self.n_ions):
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
            if "tic" in stats:
                tic = cast(zarr.Array, stats["tic"])
                if tic.ndim != 1 or int(
                    tic.shape[0]
                ) != int(self.coordinates.shape[0]):
                    raise ValueError(
                        "stats/tic must match axes/coordinates length."
                    )
            if "rms" in stats:
                rms = cast(zarr.Array, stats["rms"])
                if rms.ndim != 1 or int(
                    rms.shape[0]
                ) != int(self.coordinates.shape[0]):
                    raise ValueError(
                        "stats/rms must match axes/coordinates length."
                    )

    def close(self) -> None:
        """Release reader resources.

        Zarr groups do not hold a closeable file handle, so no explicit action
        is required.
        """

    # --- Internal selection helpers ------------------------------------------

    def _pixel_spectra_arrays(
        self,
    ) -> tuple[zarr.Array, zarr.Array, ZarrStorageMode, zarr.Array | None]:
        if self.row_axis == "pixel":
            return self.intensity, self.offsets, self.encoding, self.point_mz
        if self._spectra_group is None:
            raise ValueError("This store does not contain pixel-major spectra.")
        return (
            cast(zarr.Array, self._spectra_group["intensity"]),
            cast(zarr.Array, self._spectra_group["offsets"]),
            self.spectra_encoding,
            None,
        )

    def _read_offsets(
        self,
        offsets: zarr.Array,
        indices: np.ndarray,
    ) -> np.ndarray:
        if indices.size > 1 and bool(np.all(np.diff(indices) == 1)):
            start = int(indices[0])
            end = int(indices[-1]) + 1
            return np.asarray(offsets[start:end], dtype=np.int64)
        return np.asarray(
            offsets.get_orthogonal_selection((indices,)),
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

        lengths = ends - starts
        output_offsets = np.concatenate(([0], np.cumsum(lengths)))
        result = np.empty(int(output_offsets[-1]), dtype=array.dtype)
        order = np.argsort(indices, kind="stable")
        order = order[lengths[order] > 0]
        if order.size == 0:
            return result
        sorted_starts, sorted_ends = starts[order], ends[order]
        chunk_size = array.chunks[0]
        chunk_start = int(sorted_starts[0]) // chunk_size * chunk_size
        first = 0
        # Read each touched chunk once, then scatter into the requested row order.
        while first < order.size:
            chunk_end = min(chunk_start + chunk_size, array.shape[0])
            last = int(np.searchsorted(sorted_starts, chunk_end, side="left"))
            values = np.asarray(array[chunk_start:chunk_end])
            for position in range(first, last):
                row = order[position]
                start = max(int(starts[row]), chunk_start)
                end = min(int(ends[row]), chunk_end)
                target = int(output_offsets[row]) + start - int(starts[row])
                result[target : target + end - start] = values[
                    start - chunk_start : end - chunk_start
                ]
            first = int(np.searchsorted(sorted_ends, chunk_end, side="right"))
            if first < order.size:
                chunk_start = max(
                    chunk_end, int(sorted_starts[first]) // chunk_size * chunk_size
                )
        return result

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
            raise ValueError(
                "spatial_shape must be [height, width] with positive integers."
            )

        # --- v1.1 layout detection -------------------------------------------
        has_ion_image = ZARR_GROUP_ION_IMAGE in self.root
        has_spectra = ZARR_GROUP_SPECTRA in self.root
        if "data" in self.root:
            raise ValueError("v1.1 stores must not contain a data/ group.")
        if not has_ion_image and not has_spectra:
            raise ValueError("Zarr store requires ion_image/ or spectra/.")
        if "metadata" not in self.root:
            raise ValueError("Zarr store requires metadata/.")
        if "axes" not in self.root:
            raise ValueError("Zarr store requires axes/.")
        axes = cast(zarr.Group, self.root["axes"])
        if "coordinates" not in axes:
            raise ValueError("axes/coordinates is required.")
        coordinates = cast(zarr.Array, axes["coordinates"])
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("axes/coordinates must have shape (n_pixels, 3).")
        if np.dtype(coordinates.dtype) != np.dtype("<u4"):
            raise ValueError("axes/coordinates must use <u4.")

        self._combined = has_ion_image and has_spectra
        if has_ion_image:
            self.data = cast(zarr.Group, self.root[ZARR_GROUP_ION_IMAGE])
            self._validate_ion_image_group(self.data, axes)
        else:
            self.data = cast(zarr.Group, self.root[ZARR_GROUP_SPECTRA])
            self._validate_spectra_group(
                self.data,
                axes,
                coordinates,
                allow_processed=True,
            )
        self.intensity = cast(zarr.Array, self.data["intensity"])
        self.offsets = cast(zarr.Array, self.data["offsets"])
        self.coordinates = coordinates
        self.shared_mz_axis = (
            cast(zarr.Array, axes["mz"]) if "mz" in axes else None
        )
        self.point_mz = (
            cast(zarr.Array, self.data["mz"]) if "mz" in self.data else None
        )
        self.metadata = MetaDataZarr.from_zarr(self.root)

        if self._combined:
            self._spectra_group = cast(
                zarr.Group, self.root[ZARR_GROUP_SPECTRA]
            )
            self._validate_spectra_group(
                self._spectra_group,
                axes,
                coordinates,
                allow_processed=False,
            )

        self._validate_stats_structure(require_both=self._combined)

        encoding = str(self.data.attrs["encoding"])
        metadata_attrs = cast(zarr.Group, self.root["metadata"]).attrs
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

    def _validate_ion_image_group(
        self,
        group: zarr.Group,
        axes: zarr.Group,
    ) -> None:
        """Validate one ion-major continuous layout group."""
        if group.attrs.get("row_axis") != "ion":
            raise ValueError("ion_image attrs.row_axis must be 'ion'.")
        if group.attrs.get("encoding") != ZARR_MODE_CONTINUOUS:
            raise ValueError("ion_image attrs.encoding must be 'continuous'.")
        if "intensity" not in group or "offsets" not in group:
            raise ValueError(
                "ion_image/intensity and ion_image/offsets are required."
            )
        if "mz" in group:
            raise ValueError("ion_image must not contain ion_image/mz.")
        if "mz" not in axes:
            raise ValueError("ion_image requires axes/mz.")
        shared_mz = cast(zarr.Array, axes["mz"])
        if shared_mz.ndim != 1 or not _is_supported_mz_dtype(shared_mz.dtype):
            raise ValueError("axes/mz must be a 1-D <f4 or <f8 array.")
        if shared_mz.attrs.get("unit") != "m/z":
            raise ValueError("axes/mz attrs.unit must be 'm/z'.")
        intensity = cast(zarr.Array, group["intensity"])
        offsets = cast(zarr.Array, group["offsets"])
        if offsets.ndim != 1 or intensity.ndim != 1:
            raise ValueError("ion_image/intensity and ion_image/offsets must be 1-D.")
        if not np.can_cast(offsets.dtype, np.dtype(np.int64), casting="safe"):
            raise ValueError(
                "ion_image/offsets must use an integer dtype safe for int64."
            )
        n_rows = int(offsets.shape[0]) - 1
        if n_rows < 0:
            raise ValueError("ion_image/offsets must contain at least one value.")
        if (
            _read_array_int(offsets, 0) != 0
            or _read_array_int(offsets, -1) != int(intensity.shape[0])
        ):
            raise ValueError(
                "ion_image/offsets boundaries do not match intensity."
            )
        if n_rows != int(shared_mz.shape[0]):
            raise ValueError("ion-major row count must equal axes/mz length.")

    def _validate_spectra_group(
        self,
        group: zarr.Group,
        axes: zarr.Group,
        coordinates: zarr.Array,
        *,
        allow_processed: bool,
    ) -> None:
        """Validate one pixel-major spectra layout group."""
        if group.attrs.get("row_axis") != "pixel":
            raise ValueError("spectra attrs.row_axis must be 'pixel'.")
        encoding = group.attrs.get("encoding")
        if encoding not in {ZARR_MODE_CONTINUOUS, ZARR_MODE_PROCESSED}:
            raise ValueError(
                "spectra attrs.encoding must be 'continuous' or 'processed'."
            )
        if not allow_processed and encoding != ZARR_MODE_CONTINUOUS:
            raise ValueError("combined stores only support continuous encoding.")
        if "intensity" not in group or "offsets" not in group:
            raise ValueError("spectra/intensity and spectra/offsets are required.")
        if encoding == ZARR_MODE_CONTINUOUS:
            if "mz" not in axes or "mz" in group:
                raise ValueError(
                    "continuous spectra require axes/mz and forbid spectra/mz."
                )
        elif "mz" not in group or "mz" in axes:
            raise ValueError(
                "processed spectra require spectra/mz and forbid axes/mz."
            )
        intensity = cast(zarr.Array, group["intensity"])
        offsets = cast(zarr.Array, group["offsets"])
        if offsets.ndim != 1 or intensity.ndim != 1:
            raise ValueError("spectra/intensity and spectra/offsets must be 1-D.")
        if not np.can_cast(offsets.dtype, np.dtype(np.int64), casting="safe"):
            raise ValueError(
                "spectra/offsets must use an integer dtype safe for int64."
            )
        n_rows = int(offsets.shape[0]) - 1
        if n_rows < 0:
            raise ValueError("spectra/offsets must contain at least one value.")
        if (
            _read_array_int(offsets, 0) != 0
            or _read_array_int(offsets, -1) != int(intensity.shape[0])
        ):
            raise ValueError(
                "spectra/offsets boundaries do not match intensity."
            )
        if n_rows != int(coordinates.shape[0]):
            raise ValueError(
                "spectra row count must equal axes/coordinates count."
            )
        if encoding == ZARR_MODE_CONTINUOUS:
            shared_mz = cast(zarr.Array, axes["mz"])
            if shared_mz.ndim != 1 or not _is_supported_mz_dtype(
                shared_mz.dtype
            ):
                raise ValueError("axes/mz must be a 1-D <f4 or <f8 array.")
            if shared_mz.attrs.get("unit") != "m/z":
                raise ValueError("axes/mz attrs.unit must be 'm/z'.")
        else:
            point_mz = cast(zarr.Array, group["mz"])
            if point_mz.ndim != 1 or not _is_supported_mz_dtype(point_mz.dtype):
                raise ValueError("spectra/mz must be a 1-D <f4 or <f8 array.")
            if point_mz.attrs.get("unit") != "m/z":
                raise ValueError("spectra/mz attrs.unit must be 'm/z'.")
            if int(point_mz.shape[0]) != int(intensity.shape[0]):
                raise ValueError(
                    "spectra/mz length must equal spectra/intensity length."
                )

    def _validate_stats_structure(self, *, require_both: bool) -> None:
        """Validate present stats arrays; combined stores require both."""
        stats = (
            cast(zarr.Group, self.root["stats"]) if "stats" in self.root else None
        )
        mean_spectrum = (
            cast(zarr.Array, stats["mean_spectrum"])
            if stats is not None and "mean_spectrum" in stats
            else None
        )
        tic = (
            cast(zarr.Array, stats["tic"])
            if stats is not None and "tic" in stats
            else None
        )
        rms = (
            cast(zarr.Array, stats["rms"])
            if stats is not None and "rms" in stats
            else None
        )
        if mean_spectrum is not None and ZARR_GROUP_ION_IMAGE not in self.root:
            raise ValueError("stats/mean_spectrum requires an ion_image/ group.")
        if tic is not None and ZARR_GROUP_SPECTRA not in self.root:
            raise ValueError("stats/tic requires a spectra/ group.")
        if rms is not None and ZARR_GROUP_SPECTRA not in self.root:
            raise ValueError("stats/rms requires a spectra/ group.")
        if require_both and (mean_spectrum is None or tic is None):
            raise ValueError(
                "combined stores require stats/mean_spectrum and stats/tic."
            )
        if mean_spectrum is not None:
            if self.shared_mz_axis is None:
                raise ValueError("stats/mean_spectrum requires a shared m/z axis.")
            if mean_spectrum.ndim != 1 or int(mean_spectrum.shape[0]) != int(
                self.shared_mz_axis.shape[0]
            ):
                raise ValueError("stats/mean_spectrum must match axes/mz length.")
            if np.dtype(mean_spectrum.dtype) != np.dtype(np.float64):
                raise ValueError("stats/mean_spectrum must be float64.")
            if mean_spectrum.attrs.get("source") != f"/{ZARR_GROUP_ION_IMAGE}":
                raise ValueError("stats/mean_spectrum attrs.source must be '/ion_image'.")
            if mean_spectrum.attrs.get("reduced_axis") != "pixel":
                raise ValueError(
                    "stats/mean_spectrum attrs.reduced_axis must be 'pixel'."
                )
            if mean_spectrum.attrs.get("denominator") != "sampled_pixels":
                raise ValueError(
                    "stats/mean_spectrum attrs.denominator must be 'sampled_pixels'."
                )
        if tic is not None:
            if tic.ndim != 1 or int(tic.shape[0]) != int(
                self.coordinates.shape[0]
            ):
                raise ValueError("stats/tic must match axes/coordinates length.")
            if np.dtype(tic.dtype) != np.dtype(np.float64):
                raise ValueError("stats/tic must be float64.")
            if tic.attrs.get("source") != f"/{ZARR_GROUP_SPECTRA}":
                raise ValueError("stats/tic attrs.source must be '/spectra'.")
            if tic.attrs.get("reduced_axis") != "mz":
                raise ValueError("stats/tic attrs.reduced_axis must be 'mz'.")
        if rms is not None:
            if rms.ndim != 1 or int(rms.shape[0]) != int(
                self.coordinates.shape[0]
            ):
                raise ValueError("stats/rms must match axes/coordinates length.")
            if np.dtype(rms.dtype) != np.dtype(np.float64):
                raise ValueError("stats/rms must be float64.")
            if rms.attrs.get("source") != f"/{ZARR_GROUP_SPECTRA}":
                raise ValueError("stats/rms attrs.source must be '/spectra'.")
            if rms.attrs.get("reduced_axis") != "mz":
                raise ValueError("stats/rms attrs.reduced_axis must be 'mz'.")
            if rms.attrs.get("denominator") != "spectrum_points":
                raise ValueError(
                    "stats/rms attrs.denominator must be 'spectrum_points'."
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
