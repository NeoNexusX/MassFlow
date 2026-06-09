"""Export MassFlow ion-image backends to imzML."""

from __future__ import annotations

import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Mapping

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
)

IonImageSource = (
    IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr
)

_SOURCE_CLASSES = {
    "binary": IonImageDataManagerBinary,
    "h5": IonImageDataManagerH5,
    "zarr": IonImageDataManagerZarr,
}
_WRITER_METADATA_KEYS = {
    "spec_type",
    "scan_direction",
    "line_scan_direction",
    "scan_pattern",
    "scan_type",
    "polarity",
}


def export_imzml(
    source: str | os.PathLike[str] | IonImageSource,
    output_filepath: str | os.PathLike[str],
    *,
    source_format: str = "auto",
    temp_dir: str | os.PathLike[str] | None = None,
    ion_block_size: int = 256,
    spectrum_batch_size: int = 256,
    metadata: Mapping[str, Any] | None = None,
) -> tuple[Path, Path]:
    """Export a continuous MassFlow ion-image dataset to ``.imzML``/``.ibd``.

    The source backends store intensities in ion-major order. This function
    transposes them to a temporary pixel-major memmap, then writes spectra with
    the official :class:`pyimzml.ImzMLWriter`.

    ``metadata`` accepts the writer fields ``spec_type``, ``scan_direction``,
    ``line_scan_direction``, ``scan_pattern``, ``scan_type``, and ``polarity``.
    Fields not stored by MassFlow ion-image backends must be supplied explicitly.
    """
    ion_block_size = _positive_int(ion_block_size, "ion_block_size")
    spectrum_batch_size = _positive_int(spectrum_batch_size, "spectrum_batch_size")
    writer_metadata = _writer_metadata(metadata)
    imzml_path, ibd_path = _output_paths(output_filepath)
    imzml_path.parent.mkdir(parents=True, exist_ok=True)

    with _open_source(source, source_format) as data_manager:
        coordinates = _prepare_imzml_coordinates(data_manager.coordinates)
        if coordinates.shape[0] != data_manager.n_pixels:
            raise ValueError(
                "coordinates row count must equal the source pixel count: "
                f"coordinates={coordinates.shape[0]}, n_pixels={data_manager.n_pixels}."
            )
        mz_axis = np.asarray(data_manager.mz_axis)
        intensity_dtype = np.dtype(data_manager.intensity_dtype)
        temp_path = _create_pixel_major_temp(
            data_manager,
            temp_dir=temp_dir,
            ion_block_size=ion_block_size,
        )
        try:
            pixel_major = np.memmap(
                temp_path,
                mode="r",
                dtype=intensity_dtype,
                shape=(data_manager.n_pixels, data_manager.n_ions),
            )
            try:
                with ImzMLWriter(
                    str(imzml_path),
                    mz_dtype=np.dtype(mz_axis.dtype).type,
                    intensity_dtype=intensity_dtype.type,
                    mode="continuous",
                    **writer_metadata,
                ) as writer:
                    for start in range(0, data_manager.n_pixels, spectrum_batch_size):
                        end = min(data_manager.n_pixels, start + spectrum_batch_size)
                        spectra = np.asarray(pixel_major[start:end])
                        for coord, intensity in zip(coordinates[start:end], spectra):
                            writer.addSpectrum(mz_axis, intensity, tuple(coord))
            finally:
                del pixel_major
        finally:
            temp_path.unlink(missing_ok=True)

    return imzml_path, ibd_path


def _positive_int(value: int, name: str) -> int:
    value = int(value)
    if value <= 0:
        raise ValueError(f"{name} must be positive.")
    return value


def _writer_metadata(metadata: Mapping[str, Any] | None) -> dict[str, Any]:
    if metadata is None:
        return {}
    unknown = set(metadata) - _WRITER_METADATA_KEYS
    if unknown:
        raise ValueError(f"Unsupported imzML metadata fields: {sorted(unknown)}")
    return dict(metadata)


def _output_paths(output_filepath: str | os.PathLike[str]) -> tuple[Path, Path]:
    output_path = Path(output_filepath)
    base_path = output_path.with_suffix("") if output_path.suffix else output_path
    return base_path.with_suffix(".imzML"), base_path.with_suffix(".ibd")


@contextmanager
def _open_source(
    source: str | os.PathLike[str] | IonImageSource,
    source_format: str,
) -> Iterator[IonImageSource]:
    if isinstance(
        source,
        (IonImageDataManagerZarr, IonImageDataManagerH5, IonImageDataManagerBinary),
    ):
        yield source
        return

    source_path = Path(source)
    resolved_format = _resolve_source_format(source_path, source_format)
    data_manager = _SOURCE_CLASSES[resolved_format](source_path, mode="r")
    try:
        yield data_manager
    finally:
        data_manager.close()


def _resolve_source_format(source_path: Path, source_format: str) -> str:
    normalized = source_format.lower()
    if normalized != "auto":
        if normalized not in _SOURCE_CLASSES:
            raise ValueError(f"Unsupported source_format: {source_format!r}")
        return normalized

    suffix = source_path.suffix.lower()
    if suffix in {".h5", ".hdf5"}:
        return "h5"
    if suffix == ".zarr":
        return "zarr"
    if suffix in {".iib", ".bin", ".binary"}:
        return "binary"
    raise ValueError(
        f"Cannot infer source format from {source_path}. "
        "Use source_format='binary', 'h5', or 'zarr'."
    )


def _prepare_imzml_coordinates(coordinates: np.ndarray) -> np.ndarray:
    coords = np.asarray(coordinates)
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("coordinates must have shape (n_pixels, >=2).")

    xy = coords[:, :2].astype(np.int64, copy=True)
    if np.any(xy == 0):
        xy += 1
    z = np.ones((xy.shape[0], 1), dtype=np.int64)
    return np.hstack((xy, z))


def _create_pixel_major_temp(
    source: IonImageSource,
    *,
    temp_dir: str | os.PathLike[str] | None,
    ion_block_size: int,
) -> Path:
    temp_root = None if temp_dir is None else str(temp_dir)
    fd, temp_name = tempfile.mkstemp(
        prefix="massflow_export_imzml_", suffix=".bin", dir=temp_root
    )
    os.close(fd)
    temp_path = Path(temp_name)

    pixel_major = np.memmap(
        temp_path,
        mode="w+",
        dtype=source.intensity_dtype,
        shape=(source.n_pixels, source.n_ions),
    )
    try:
        for start in range(0, source.n_ions, ion_block_size):
            end = min(source.n_ions, start + ion_block_size)
            ion_block = _read_ion_block(source, start, end)
            if ion_block.shape != (end - start, source.n_pixels):
                raise ValueError(
                    "Ion-image source returned an unexpected block shape: "
                    f"expected={(end - start, source.n_pixels)}, actual={ion_block.shape}."
                )
            pixel_major[:, start:end] = ion_block.T
        pixel_major.flush()
    except BaseException:
        temp_path.unlink(missing_ok=True)
        raise
    finally:
        del pixel_major

    return temp_path


def _read_ion_block(source: IonImageSource, start: int, end: int) -> np.ndarray:
    if isinstance(source, (IonImageDataManagerZarr, IonImageDataManagerH5)):
        dense = np.asarray(source.dset[start:end])
        return dense.reshape(end - start, -1)[:, source.flat_idx]
    return np.asarray(source.intensity_matrix[start:end])
