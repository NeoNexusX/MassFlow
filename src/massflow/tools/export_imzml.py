"""Export full or mask-selected MassFlow datasets with the batch imzML writer."""

# pylint: disable=too-many-arguments,too-many-locals

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any, Iterator, Mapping

import numpy as np

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
    MSDataManagerZarr,
)
from massflow.msi_zarr import ZARR_MODE_CONTINUOUS
from massflow.tools.funs import strip_known_suffix
from massflow.tools.stream_imzml_writer import ImzMLWriter

IonImageSource = (
    IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr
)
ExportSource = IonImageSource | MSDataManagerZarr | MSDataManagerImzML

# ---------------------------------------------------------------------------
# Positioned I/O and tile layout
# ---------------------------------------------------------------------------

_OS_PWRITE = getattr(os, "pwrite", None)  # pylint: disable=invalid-name


def _write_at(fd: int, data: np.ndarray, offset: int) -> None:
    view = memoryview(data).cast("B")
    if _OS_PWRITE is not None:
        while view:
            written = _OS_PWRITE(fd, view, offset)
            if written == 0:
                raise OSError("Positioned write returned zero bytes.")
            view = view[written:]
            offset += written
        return

    os.lseek(fd, offset, os.SEEK_SET)
    while view:
        written = os.write(fd, view)
        if written == 0:
            raise OSError("Write returned zero bytes.")
        view = view[written:]


def _pixel_batches(
    n_pixels: int,
    batch_size: int,
) -> Iterator[tuple[int, int]]:
    for start in range(0, n_pixels, batch_size):
        yield start, min(n_pixels, start + batch_size)


def _tile_byte_offset(
    *,
    pixel_start: int,
    batch_pixels: int,
    n_ions: int,
    ion_start: int,
    itemsize: int,
) -> int:
    return (pixel_start * n_ions + ion_start * batch_pixels) * itemsize


# ---------------------------------------------------------------------------
# Public export API
# ---------------------------------------------------------------------------


def export_imzml(
    source: ExportSource,
    output_path: str | os.PathLike[str],
    *,
    temp_dir: str | os.PathLike[str] | None = None,
    batch_size: int = 256,
    mask: np.ndarray | None = None,
) -> Path:
    """Export pixels to ``.imzML``/``.ibd`` in batches and return the imzML path.

    Ion-image backends use a tiled transpose file, reading source ion blocks
    once. Spectra Zarr and uncompressed imzML stream through ``flat_generator``.
    An imzML source requires loaded mode metadata and, in continuous mode, a
    shared m/z axis (call ``preload_hook()`` or ``load_head_data()`` first).
    ``mask`` must match the original spatial canvas; coordinates keep their
    values and order. Source spectrum type, polarity and supported scan settings
    are passed to the existing writer; other global metadata is not exported.

    The project writer uses zero spectrum statistics and a placeholder SHA-1,
    writes uncompressed blocks, and retains O(output spectra) metadata for XML
    rendering. The source is never finalized or closed by this function.
    """
    batch_size = _positive_int(batch_size, "batch_size")
    _validate_source(source)
    imzml_path, ibd_path = _output_paths(output_path)
    if isinstance(source, MSDataManagerImzML):
        source_paths = {
            Path(source.imzml_filepath).resolve(),
            Path(source.ibd_filepath).resolve(),
        }
        if any(
            path.resolve() in source_paths
            or (
                path.exists()
                and any(path.samefile(original) for original in source_paths)
            )
            for path in (imzml_path, ibd_path)
        ):
            raise ValueError("imzML output must not overwrite the source file pair.")

    if mask is not None:
        width, height = _spatial_dimensions(source)
        mask = _normalize_mask(mask, width=width, height=height)

    if isinstance(source, (MSDataManagerZarr, MSDataManagerImzML)):
        writer_mode, coordinate_base = _spectra_export_context(source, mask)
        mz_dtype = np.dtype(source.mz_dtype)
        intensity_dtype = np.dtype(source.intensity_dtype)
    else:
        coordinates = _prepare_imzml_coordinates(source.coordinates, source.n_pixels)
        mz_axis = np.asarray(source.mz_axis)
        if mz_axis.ndim != 1 or mz_axis.size != source.n_ions:
            raise ValueError(
                "mz_axis must be one-dimensional and match n_ions: "
                f"shape={mz_axis.shape}, n_ions={source.n_ions}."
            )
        mz_dtype = mz_axis.dtype
        intensity_dtype = np.dtype(source.intensity_dtype)
        selected_indices = None
        if mask is not None:
            coordinate_base = (
                "one" if isinstance(source, IonImageDataManagerZarr) else None
            )
            selected_indices = _selected_indices(mask, coordinates, coordinate_base)
            coordinates = coordinates[selected_indices]
        if coordinates.shape[0] == 0 or source.n_ions == 0:
            raise ValueError("source selection contains no spectra.")

    # Validate writer encodings before opening or cleaning up output files.
    ImzMLWriter._np_type_to_name(mz_dtype)  # pylint: disable=protected-access
    ImzMLWriter._np_type_to_name(intensity_dtype)  # pylint: disable=protected-access
    writer_metadata = _writer_metadata_from_source(source)
    imzml_path.parent.mkdir(parents=True, exist_ok=True)
    if temp_dir is not None:
        Path(temp_dir).mkdir(parents=True, exist_ok=True)

    try:
        if isinstance(source, (MSDataManagerZarr, MSDataManagerImzML)):
            _export_spectra(
                source=source,
                imzml_path=imzml_path,
                batch_size=batch_size,
                writer_metadata=writer_metadata,
                writer_mode=writer_mode,  # type: ignore[possibly-undefined]
                mask=mask,
                coordinate_base=coordinate_base,  # type: ignore[possibly-undefined]
            )
        else:
            _export_ion_image(
                source=source,
                imzml_path=imzml_path,
                temp_dir=temp_dir,
                batch_size=batch_size,
                writer_metadata=writer_metadata,
                coordinates=coordinates,  # type: ignore[possibly-undefined]
                mz_axis=mz_axis,  # type: ignore[possibly-undefined]
                selected_indices=selected_indices,  # type: ignore[possibly-undefined]
            )
    except BaseException:
        _remove_incomplete_output(imzml_path, ibd_path)
        raise

    return imzml_path


def _validate_source(source: ExportSource) -> None:
    if not isinstance(
        source,
        (
            IonImageDataManagerBinary,
            IonImageDataManagerH5,
            IonImageDataManagerZarr,
            MSDataManagerZarr,
            MSDataManagerImzML,
        ),
    ):
        raise TypeError("source must be a supported MassFlow data manager.")

    if isinstance(source, MSDataManagerImzML):
        meta = source.ms.meta
        if bool(meta.continuous) == bool(meta.processed):
            raise ValueError(
                "imzML source requires exactly one loaded continuous/processed mode; "
                "call preload_hook() or load_head_data() first."
            )
        if meta.continuous and source.ms.shared_mz_list is None:
            raise ValueError(
                "continuous imzML source requires an initialized shared m/z axis; "
                "call preload_hook() or load_head_data() first."
            )
        return
    if source.mode == "w" and not source._finalized:  # pylint: disable=protected-access
        raise ValueError("source writer must be finalized before export.")


def _remove_incomplete_output(
    imzml_path: Path,
    ibd_path: Path,
) -> None:
    imzml_path.unlink(missing_ok=True)
    ibd_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Batch imzML writer output
# ---------------------------------------------------------------------------


def _export_ion_image(
    *,
    source: IonImageSource,
    imzml_path: Path,
    temp_dir: str | os.PathLike[str] | None,
    batch_size: int,
    writer_metadata: Mapping[str, Any],
    coordinates: np.ndarray,
    mz_axis: np.ndarray,
    selected_indices: np.ndarray | None,
) -> None:
    intensity_dtype = np.dtype(source.intensity_dtype)
    temp_path = _create_tiled_temp(
        source,
        temp_dir=temp_dir,
        batch_size=batch_size,
        selected_indices=selected_indices,
    )
    try:
        _write_ion_image_imzml(
            imzml_path=imzml_path,
            temp_path=temp_path,
            mz_axis=mz_axis,
            intensity_dtype=intensity_dtype,
            coordinates=coordinates,
            n_pixels=coordinates.shape[0],
            n_ions=source.n_ions,
            batch_size=batch_size,
            writer_metadata=writer_metadata,
        )
    finally:
        temp_path.unlink(missing_ok=True)


def _write_ion_image_imzml(
    *,
    imzml_path: Path,
    temp_path: Path,
    mz_axis: np.ndarray,
    intensity_dtype: np.dtype,
    coordinates: np.ndarray,
    n_pixels: int,
    n_ions: int,
    batch_size: int,
    writer_metadata: Mapping[str, Any],
) -> None:
    with ImzMLWriter(
        str(imzml_path),
        mz_dtype=mz_axis.dtype,
        intensity_dtype=intensity_dtype,
        mode="continuous",
        compute_stats=False,
        **writer_metadata,
    ) as writer:
        for pixel_start, pixel_end in _pixel_batches(
            n_pixels,
            batch_size,
        ):
            batch_pixels = pixel_end - pixel_start
            batch_offset = pixel_start * n_ions * intensity_dtype.itemsize
            ion_major_batch = np.memmap(
                temp_path,
                mode="r",
                dtype=intensity_dtype,
                offset=batch_offset,
                shape=(n_ions, batch_pixels),
            )
            spectra = np.ascontiguousarray(ion_major_batch.T)
            try:
                writer.swap_flat_out2disk(
                    mz_flat=mz_axis,
                    intensity_flat=spectra.reshape(-1),
                    lengths=np.full(batch_pixels, n_ions, dtype=np.int64),
                    coords=coordinates[pixel_start:pixel_end],
                )
            finally:
                del spectra
                del ion_major_batch


# ---------------------------------------------------------------------------
# Direct spectra stream export
# ---------------------------------------------------------------------------


def _spectra_export_context(
    source: MSDataManagerZarr | MSDataManagerImzML,
    mask: np.ndarray | None,
) -> tuple[str, str | None]:
    if isinstance(source, MSDataManagerZarr):
        if source._target_indices().size == 0:  # pylint: disable=protected-access
            raise ValueError("source selection contains no spectra.")
        mode = (
            "continuous" if source.storage_mode == ZARR_MODE_CONTINUOUS else "processed"
        )
        return mode, "one"
    mode = "continuous" if source.ms.meta.continuous else "processed"
    coordinate_base = None
    if mask is not None:
        xy = np.asarray(source.parser.coordinates)[:, :2]
        coordinate_base = "zero" if bool(np.any(xy == 0)) else "one"
    return mode, coordinate_base


def _export_spectra(
    *,
    source: MSDataManagerZarr | MSDataManagerImzML,
    imzml_path: Path,
    batch_size: int,
    writer_metadata: Mapping[str, Any],
    writer_mode: str,
    mask: np.ndarray | None,
    coordinate_base: str | None,
) -> None:
    with ImzMLWriter(
        str(imzml_path),
        mz_dtype=np.dtype(source.mz_dtype),
        intensity_dtype=np.dtype(source.intensity_dtype),
        mode=writer_mode,
        compute_stats=False,
        **writer_metadata,
    ) as writer:
        written = 0
        for mz_data, intensity_flat, lengths, coordinates in source.flat_generator(
            batch_size=batch_size,
            include_mz=True,
        ):
            lengths = np.asarray(lengths)
            if (
                lengths.ndim != 1
                or lengths.dtype.kind not in "iu"
                or np.any(lengths < 0)
            ):
                raise ValueError(
                    "lengths must be a one-dimensional array of nonnegative integers."
                )
            lengths = lengths.astype(np.int64, copy=False)
            coordinates = _prepare_imzml_coordinates(coordinates, lengths.size)
            intensity_flat = np.asarray(intensity_flat, dtype=source.intensity_dtype)
            total_points = int(lengths.sum(dtype=np.int64))
            if intensity_flat.ndim != 1 or intensity_flat.size != total_points:
                raise ValueError("intensity_flat size must equal sum(lengths).")
            mz_data = np.asarray(mz_data, dtype=source.mz_dtype)
            if writer_mode == "continuous":
                if mz_data.ndim != 1:
                    raise ValueError(
                        "continuous spectra require a shared 1-D m/z axis."
                    )
                if np.any(lengths != mz_data.size):
                    raise ValueError(
                        "continuous spectrum lengths must match the shared m/z axis."
                    )
            elif mz_data.ndim != 1 or mz_data.size != total_points:
                raise ValueError("processed m/z size must equal sum(lengths).")

            if mask is not None:
                selected = _select_coordinates(mask, coordinates, coordinate_base)
                if not bool(np.any(selected)):
                    continue
                if not bool(np.all(selected)):
                    point_selected = np.repeat(selected, lengths)
                    intensity_flat = intensity_flat[point_selected]
                    if writer_mode == "processed":
                        mz_data = mz_data[point_selected]
                    lengths = lengths[selected]
                    coordinates = coordinates[selected]

            writer.swap_flat_out2disk(
                mz_flat=mz_data,
                intensity_flat=intensity_flat,
                lengths=lengths,
                coords=coordinates,
            )
            written += lengths.size
        if written == 0:
            message = (
                "source selection contains no spectra."
                if mask is None
                else "mask does not retain any source pixels."
            )
            raise ValueError(message)


# --- Source opening and input validation -----------------------------------
def _positive_int(value: int, name: str) -> int:
    value = int(value)
    if value <= 0:
        raise ValueError(f"{name} must be positive.")
    return value


def _source_metadata(source: ExportSource):
    if isinstance(source, MSDataManagerImzML):
        return source.ms.meta
    if isinstance(source, (MSDataManagerZarr, IonImageDataManagerZarr)):
        return source._require_reader().metadata  # pylint: disable=protected-access
    return getattr(source, "meta", None)


def _writer_metadata_from_source(source: ExportSource) -> dict[str, Any]:
    metadata = _source_metadata(source)
    writer_metadata: dict[str, Any] = {}
    if metadata is not None:
        if metadata.get("profile_spectrum") is True:
            writer_metadata["spec_type"] = "profile"
        elif metadata.get("centroid_spectrum") is True:
            writer_metadata["spec_type"] = "centroid"
        scan_values = {
            "scan_direction": {"top_down", "bottom_up", "left_right", "right_left"},
            "line_scan_direction": {
                "line_left_right",
                "line_right_left",
                "line_top_down",
                "line_bottom_up",
            },
            "scan_pattern": {"one_way", "meandering", "random_access"},
            "scan_type": {"horizontal_line", "vertical_line"},
            "polarity": {"positive", "negative"},
        }
        for key, supported in scan_values.items():
            value = metadata.get(key)
            if isinstance(value, str) and value in supported:
                writer_metadata[key] = value
    if isinstance(source, MSDataManagerImzML):
        parser = source.parser
        if parser.spectrum_mode is not None:
            writer_metadata["spec_type"] = parser.spectrum_mode
        if parser.polarity in {"positive", "negative"}:
            writer_metadata["polarity"] = parser.polarity
    return writer_metadata


def _output_paths(output_path: str | os.PathLike[str]) -> tuple[Path, Path]:
    base = strip_known_suffix(output_path)
    return Path(f"{base}.imzML"), Path(f"{base}.ibd")


def _prepare_imzml_coordinates(
    coordinates: np.ndarray,
    n_pixels: int,
) -> np.ndarray:
    coords = np.asarray(coordinates)
    if coords.ndim != 2 or coords.shape[1] not in {2, 3}:
        raise ValueError("coordinates must have shape (n_pixels, 2) or (n_pixels, 3).")
    if coords.shape[0] != n_pixels:
        raise ValueError(
            "coordinates row count must equal the source pixel count: "
            f"coordinates={coords.shape[0]}, n_pixels={n_pixels}."
        )
    return coords


# ---------------------------------------------------------------------------
# Spatial mask helpers
# ---------------------------------------------------------------------------


def _spatial_dimensions(data_manager: ExportSource) -> tuple[int, int]:
    """Return the source image width and height."""
    if isinstance(data_manager, MSDataManagerImzML):
        meta = data_manager.ms.meta
        width = int(meta.max_count_of_pixels_x or 0)  # type: ignore
        height = int(meta.max_count_of_pixels_y or 0)  # type: ignore
        if width > 0 and height > 0:
            return width, height

        coordinates = np.asarray(data_manager.parser.coordinates)
        if coordinates.size == 0:
            raise ValueError("cannot infer image dimensions without coordinates.")
        xy = coordinates[:, :2].astype(np.int64, copy=False)
        zero_based = bool(np.any(xy == 0))
        if width <= 0:
            width = int(xy[:, 0].max()) + int(zero_based)
        if height <= 0:
            height = int(xy[:, 1].max()) + int(zero_based)
        return width, height
    if isinstance(data_manager, MSDataManagerZarr):
        reader = data_manager._require_reader()  # pylint: disable=protected-access
        height, width = reader.spatial_shape
        return width, height
    return int(data_manager.width), int(data_manager.height)


def _normalize_mask(mask: np.ndarray, *, width: int, height: int) -> np.ndarray:
    """Validate and normalize a spatial mask."""
    normalized = np.asarray(mask)
    if normalized.shape != (height, width):
        raise ValueError(
            f"mask must have shape {(height, width)}, got {normalized.shape}."
        )
    return normalized.astype(bool, copy=False)


def _select_coordinates(
    mask: np.ndarray,
    coordinates: np.ndarray,
    coordinate_base: str | None = None,
) -> np.ndarray:
    """Return which coordinates fall inside the mask."""
    xy = np.asarray(coordinates)[:, :2].astype(np.int64, copy=False)
    if coordinate_base == "zero":
        xy0 = xy
    elif coordinate_base == "one":
        xy0 = xy - 1
    else:
        xy0 = xy if bool(np.any(xy == 0)) else xy - 1
    return mask[xy0[:, 1], xy0[:, 0]]


def _filter_flat_batch(
    values: np.ndarray,
    lengths: np.ndarray,
    selected: np.ndarray,
) -> np.ndarray:
    """Filter a flattened variable-length spectrum batch by pixel."""
    flat_selected = np.repeat(selected, lengths.astype(np.int64, copy=False))
    return np.asarray(values)[flat_selected]


def _selected_indices(
    mask: np.ndarray,
    coordinates: np.ndarray,
    coordinate_base: str | None = None,
) -> np.ndarray:
    indices = np.flatnonzero(_select_coordinates(mask, coordinates, coordinate_base))
    if indices.size == 0:
        raise ValueError("mask does not retain any source pixels.")
    return indices


# ---------------------------------------------------------------------------
# Tiled temporary-file transpose
# ---------------------------------------------------------------------------


def _create_tiled_temp(
    source: IonImageSource,
    *,
    temp_dir: str | os.PathLike[str] | None,
    batch_size: int,
    selected_indices: np.ndarray | None = None,
) -> Path:
    temp_root = None if temp_dir is None else str(temp_dir)
    fd, temp_name = tempfile.mkstemp(
        prefix="massflow_export_imzml_",
        suffix=".bin",
        dir=temp_root,
    )
    temp_path = Path(temp_name)
    n_pixels = source.n_pixels if selected_indices is None else selected_indices.size
    total_bytes = n_pixels * source.n_ions * np.dtype(source.intensity_dtype).itemsize

    try:
        try:
            os.ftruncate(fd, total_bytes)
            block_rows = _target_block_rows(source, batch_size)

            for ion_start in range(0, source.n_ions, block_rows):
                ion_end = min(source.n_ions, ion_start + block_rows)
                ion_block = _read_ion_block(source, ion_start, ion_end)
                expected_shape = (ion_end - ion_start, source.n_pixels)
                if ion_block.shape != expected_shape:
                    raise ValueError(
                        "Ion-image source returned an unexpected block shape: "
                        f"expected={expected_shape}, actual={ion_block.shape}."
                    )

                for pixel_start, pixel_end in _pixel_batches(
                    n_pixels,
                    batch_size,
                ):
                    batch_pixels = pixel_end - pixel_start
                    columns = (
                        slice(pixel_start, pixel_end)
                        if selected_indices is None
                        else selected_indices[pixel_start:pixel_end]
                    )
                    tile = np.ascontiguousarray(ion_block[:, columns])
                    offset = _tile_byte_offset(
                        pixel_start=pixel_start,
                        batch_pixels=batch_pixels,
                        n_ions=source.n_ions,
                        ion_start=ion_start,
                        itemsize=np.dtype(source.intensity_dtype).itemsize,
                    )
                    _write_at(fd, tile, offset)
        finally:
            os.close(fd)
    except BaseException:
        temp_path.unlink(missing_ok=True)
        raise

    return temp_path


# --- Backend-aware block sizing and source reads ---------------------------
def _aligned_block_rows(target_rows: int, chunk_rows: int) -> int:
    target_rows = max(1, int(target_rows))
    chunk_rows = max(1, int(chunk_rows))
    if target_rows < chunk_rows:
        return chunk_rows
    return max(chunk_rows, (target_rows // chunk_rows) * chunk_rows)


def _target_block_rows(
    source: IonImageSource,
    batch_size: int,
) -> int:
    target_elements = batch_size * source.n_ions

    if isinstance(source, IonImageDataManagerH5):
        row_elements = source.height * source.width
        assert source._dset is not None  # pylint: disable=protected-access
        chunks = source._dset.chunks  # pylint: disable=protected-access
        chunk_rows = 1 if chunks is None else int(chunks[0])
    elif isinstance(source, IonImageDataManagerZarr):
        reader = source._require_reader()  # pylint: disable=protected-access
        row_elements = source.n_pixels
        chunk_points = int(reader.intensity.chunks[0])
        chunk_rows = max(1, chunk_points // source.n_pixels)
    else:
        row_elements = source.n_pixels
        chunk_rows = 1

    target_rows = max(1, target_elements // max(1, row_elements))
    return min(source.n_ions, _aligned_block_rows(target_rows, chunk_rows))


def _read_ion_block(source: IonImageSource, start: int, end: int) -> np.ndarray:
    return source.get_ion_block(start, end)
