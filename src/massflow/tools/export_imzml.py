"""Export MassFlow ion-image backends to imzML."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any, Iterator, Mapping

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerZarr,
)
from massflow.msi_zarr import (
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
)

IonImageSource = (
    IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr
)
ExportSource = IonImageSource | MSDataManagerZarr
# ---------------------------------------------------------------------------
# Positioned I/O and tile layout
# ---------------------------------------------------------------------------

_OS_PWRITE = getattr(os, "pwrite", None)


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
) -> tuple[Path, Path]:
    """Export a MassFlow dataset to ``.imzML``/``.ibd`` via official pyimzML.

    Ion-image backends (Binary/H5/Zarr) are transposed to a tiled temporary
    file and written in continuous mode. Pixel-major ``MSDataManagerZarr`` is
    streamed directly through ``flat_generator()`` in continuous or processed
    mode according to the Zarr encoding. The source manager is never opened,
    finalized, or closed by this function.
    """
    batch_size = _positive_int(batch_size, "batch_size")
    _validate_source(source)
    imzml_path, ibd_path = _output_paths(output_path)
    imzml_path.parent.mkdir(parents=True, exist_ok=True)

    writer_metadata = _writer_metadata_from_source(source)
    try:
        if isinstance(source, MSDataManagerZarr):
            _export_spectra_zarr(
                source=source,
                imzml_path=imzml_path,
                batch_size=batch_size,
                writer_metadata=writer_metadata,
            )
        else:
            _export_ion_image(
                source=source,
                imzml_path=imzml_path,
                temp_dir=temp_dir,
                batch_size=batch_size,
                writer_metadata=writer_metadata,
            )
    except BaseException:
        _remove_incomplete_output(imzml_path, ibd_path)
        raise

    return imzml_path, ibd_path


def _validate_source(source: ExportSource) -> None:
    if not isinstance(
        source,
        (
            IonImageDataManagerBinary,
            IonImageDataManagerH5,
            IonImageDataManagerZarr,
            MSDataManagerZarr,
        ),
    ):
        raise TypeError("source must be a supported MassFlow data manager.")

    if source.mode == "w" and not source._finalized:  # pylint: disable=protected-access
        raise ValueError("source writer must be finalized before export.")


def _remove_incomplete_output(
    imzml_path: Path,
    ibd_path: Path,
) -> None:
    imzml_path.unlink(missing_ok=True)
    ibd_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Official pyimzML writer output
# ---------------------------------------------------------------------------


def _export_ion_image(
    *,
    source: IonImageSource,
    imzml_path: Path,
    temp_dir: str | os.PathLike[str] | None,
    batch_size: int,
    writer_metadata: Mapping[str, Any],
) -> None:
    coordinates = _prepare_imzml_coordinates(
        source.coordinates,
        source.n_pixels,
    )
    mz_axis = np.asarray(source.mz_axis)
    if mz_axis.ndim != 1 or mz_axis.size != source.n_ions:
        raise ValueError(
            "mz_axis must be one-dimensional and match n_ions: "
            f"shape={mz_axis.shape}, n_ions={source.n_ions}."
        )
    intensity_dtype = np.dtype(source.intensity_dtype)
    temp_path = _create_tiled_temp(
        source,
        temp_dir=temp_dir,
        batch_size=batch_size,
    )
    try:
        _write_ion_image_imzml(
            imzml_path=imzml_path,
            temp_path=temp_path,
            mz_axis=mz_axis,
            intensity_dtype=intensity_dtype,
            coordinates=coordinates,
            n_pixels=source.n_pixels,
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
        mz_dtype=np.dtype(mz_axis.dtype).type,
        intensity_dtype=intensity_dtype.type,
        mode="continuous",
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
                for coord, intensity in zip(
                    coordinates[pixel_start:pixel_end],
                    spectra,
                ):
                    writer.addSpectrum(
                        mz_axis,
                        intensity,
                        tuple(coord),
                    )
            finally:
                del spectra
                del ion_major_batch


# ---------------------------------------------------------------------------
# Direct MSDataManagerZarr stream export
# ---------------------------------------------------------------------------


def _flat_offsets(lengths: np.ndarray) -> np.ndarray:
    lengths = np.asarray(lengths, dtype=np.int64)
    offsets = np.empty(lengths.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(lengths, dtype=np.int64, out=offsets[1:])
    return offsets


def _export_spectra_zarr(
    *,
    source: MSDataManagerZarr,
    imzml_path: Path,
    batch_size: int,
    writer_metadata: Mapping[str, Any],
) -> None:
    if source._target_indices().size == 0:  # pylint: disable=protected-access
        raise ValueError("source selection contains no spectra.")

    storage_mode = source.storage_mode
    writer_mode = "continuous" if storage_mode == ZARR_MODE_CONTINUOUS else "processed"

    with ImzMLWriter(
        str(imzml_path),
        mz_dtype=np.dtype(source.mz_dtype).type,
        intensity_dtype=np.dtype(source.intensity_dtype).type,
        mode=writer_mode,
        **writer_metadata,
    ) as writer:
        for mz_data, intensity_flat, lengths, coordinates in source.flat_generator(
            batch_size=batch_size,
            include_mz=True,
        ):
            lengths = np.asarray(lengths, dtype=np.int64)
            coordinates = _prepare_imzml_coordinates(
                coordinates,
                lengths.size,
            )
            intensity_flat = np.asarray(
                intensity_flat,
                dtype=source.intensity_dtype,
            )
            offsets = _flat_offsets(lengths)
            if intensity_flat.ndim != 1 or intensity_flat.size != offsets[-1]:
                raise ValueError("intensity_flat size must equal sum(lengths).")

            if storage_mode == ZARR_MODE_CONTINUOUS:
                shared_mz = np.asarray(mz_data, dtype=source.mz_dtype)
                if shared_mz.ndim != 1:
                    raise ValueError(
                        "continuous spectra Zarr requires a shared 1-D m/z axis."
                    )
                if np.any(lengths != shared_mz.size):
                    raise ValueError(
                        "continuous spectrum lengths must match the shared m/z axis."
                    )
                for index, coord in enumerate(coordinates):
                    start = int(offsets[index])
                    end = int(offsets[index + 1])
                    writer.addSpectrum(
                        shared_mz,
                        intensity_flat[start:end],
                        tuple(coord),
                    )
                continue

            point_mz = np.asarray(mz_data, dtype=source.mz_dtype)
            if point_mz.ndim != 1 or point_mz.size != offsets[-1]:
                raise ValueError("processed m/z size must equal sum(lengths).")
            for index, coord in enumerate(coordinates):
                start = int(offsets[index])
                end = int(offsets[index + 1])
                writer.addSpectrum(
                    point_mz[start:end],
                    intensity_flat[start:end],
                    tuple(coord),
                )


# --- Source opening and input validation -----------------------------------
def _positive_int(value: int, name: str) -> int:
    value = int(value)
    if value <= 0:
        raise ValueError(f"{name} must be positive.")
    return value


def _source_metadata(source: ExportSource):
    if isinstance(source, MSDataManagerZarr):
        return source._require_reader().metadata  # pylint: disable=protected-access
    return getattr(source, "meta", None)


def _writer_metadata_from_source(source: ExportSource) -> dict[str, Any]:
    metadata = _source_metadata(source)
    if metadata is None:
        return {}

    writer_metadata: dict[str, Any] = {}
    if metadata.get("profile_spectrum") is True:
        writer_metadata["spec_type"] = "profile"
    elif metadata.get("centroid_spectrum") is True:
        writer_metadata["spec_type"] = "centroid"

    for key in (
        "scan_direction",
        "line_scan_direction",
        "scan_pattern",
        "scan_type",
        "polarity",
    ):
        value = metadata.get(key)
        if value is not None:
            writer_metadata[key] = value

    return writer_metadata


def _output_paths(output_path: str | os.PathLike[str]) -> tuple[Path, Path]:
    output_path = Path(output_path)
    base_path = output_path.with_suffix("") if output_path.suffix else output_path
    return base_path.with_suffix(".imzML"), base_path.with_suffix(".ibd")


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
# Tiled temporary-file transpose
# ---------------------------------------------------------------------------


def _create_tiled_temp(
    source: IonImageSource,
    *,
    temp_dir: str | os.PathLike[str] | None,
    batch_size: int,
) -> Path:
    temp_root = None if temp_dir is None else str(temp_dir)
    fd, temp_name = tempfile.mkstemp(
        prefix="massflow_export_imzml_",
        suffix=".bin",
        dir=temp_root,
    )
    temp_path = Path(temp_name)
    total_bytes = (
        source.n_pixels * source.n_ions * np.dtype(source.intensity_dtype).itemsize
    )

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
                    source.n_pixels,
                    batch_size,
                ):
                    batch_pixels = pixel_end - pixel_start
                    tile = np.ascontiguousarray(ion_block[:, pixel_start:pixel_end])
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
