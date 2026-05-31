"""Unified ion-image and TIC-image plotting helpers."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator, Literal

import matplotlib.pyplot as plt
import numpy as np

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
)
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tools")

InputFormat = Literal["imzml", "ion_image_binary", "ion_image_h5", "ion_image_zarr"]
_ION_IMAGE_FORMATS = {"ion_image_binary", "ion_image_h5", "ion_image_zarr"}


def _infer_format_from_path(path: Path) -> InputFormat:
    suffix = path.suffix.lower()
    if suffix == ".imzml":
        return "imzml"
    if suffix in {".iib", ".bin"}:
        return "ion_image_binary"
    if suffix in {".h5", ".hdf5"}:
        return "ion_image_h5"
    if suffix == ".zarr":
        return "ion_image_zarr"
    raise ValueError(f"Cannot infer input format from path suffix: {path}")


def _infer_format_from_manager(source: Any) -> InputFormat:
    if isinstance(source, MSDataManagerImzML):
        return "imzml"
    if isinstance(source, IonImageDataManagerZarr):
        return "ion_image_zarr"
    if isinstance(source, IonImageDataManagerH5):
        return "ion_image_h5"
    if isinstance(source, IonImageDataManagerBinary):
        return "ion_image_binary"
    if hasattr(source, "flat_generator") and hasattr(source, "load_head_data"):
        return "imzml"
    if hasattr(source, "get_ion_image"):
        return "ion_image_binary"
    raise TypeError("source must be a supported path or data manager.")


@contextmanager
def _resolve_source(
    source: str | Path | Any,
) -> Generator[tuple[Any, InputFormat], None, None]:
    if isinstance(source, (str, Path)):
        path = Path(source)
        source_format = _infer_format_from_path(path)
        if source_format == "imzml":
            manager = MSDataManagerImzML(filepath=str(path))
            _ensure_imzml_loaded(manager)
        elif source_format == "ion_image_binary":
            manager = IonImageDataManagerBinary(filepath=path, mode="r")
        elif source_format == "ion_image_h5":
            manager = IonImageDataManagerH5(filepath=path, mode="r")
        elif source_format == "ion_image_zarr":
            manager = IonImageDataManagerZarr(filepath=path, mode="r")
        else:
            raise ValueError(f"Unsupported input format: {source_format}")

        try:
            yield manager, source_format
        finally:
            close = getattr(manager, "close", None)
            if callable(close):
                close()
        return

    source_format = _infer_format_from_manager(source)
    if source_format == "imzml":
        _ensure_imzml_loaded(source)
    yield source, source_format


def _ensure_imzml_loaded(manager: Any) -> None:
    if len(manager.ms) == 0:
        manager.load_head_data()


def _imzml_is_continuous(manager: Any) -> bool | None:
    meta = getattr(getattr(manager, "ms", None), "meta", None)
    continuous = getattr(meta, "continuous", None)
    return None if continuous is None else bool(continuous)


def _batch_uses_shared_mz(
    manager: Any,
    mz_data: np.ndarray | None,
    lengths: np.ndarray,
) -> bool:
    continuous = _imzml_is_continuous(manager)
    if continuous is not None:
        return continuous
    if mz_data is None or mz_data.ndim != 1 or lengths.size == 0:
        return False
    if mz_data.size != int(np.sum(lengths, dtype=np.int64)):
        return True
    return bool(lengths.size > 1 and np.all(lengths == mz_data.size))


def _nearest_mz_index(
    mz_axis: np.ndarray,
    mz: float,
) -> int:
    if mz_axis.ndim != 1 or mz_axis.size == 0:
        raise ValueError("mz_axis must be a non-empty 1-D array.")
    return int(np.nanargmin(np.abs(mz_axis.astype(float, copy=False) - float(mz))))


def _mz_window_indices(
    mz_axis: np.ndarray,
    mz: float,
    mz_tolerance: float,
) -> np.ndarray:
    low = float(mz) - float(mz_tolerance)
    high = float(mz) + float(mz_tolerance)
    return np.flatnonzero((mz_axis >= low) & (mz_axis <= high))


def _sum_by_lengths(intensity_flat: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    dtype = np.result_type(intensity_flat.dtype, np.float64)
    if lengths.size == 0:
        return np.zeros(0, dtype=dtype)

    offsets = np.empty(lengths.size, dtype=np.int64)
    offsets[0] = 0
    if lengths.size > 1:
        np.cumsum(lengths[:-1], dtype=np.int64, out=offsets[1:])

    values = np.zeros(lengths.size, dtype=dtype)
    nonzero = lengths > 0
    if np.any(nonzero):
        values[nonzero] = np.add.reduceat(intensity_flat, offsets[nonzero])
    return values


def _sum_processed_mz_values(
    mz_flat: np.ndarray | None,
    intensity_flat: np.ndarray,
    lengths: np.ndarray,
    mz: float,
    mz_tolerance: float | None,
) -> np.ndarray:
    if mz_flat is None:
        raise ValueError("Processed imzML ion image requires m/z data.")
    if mz_flat.size != int(np.sum(lengths, dtype=np.int64)):
        raise ValueError("Processed imzML ion image requires flattened per-spectrum m/z data.")

    values = np.zeros(lengths.size, dtype=np.result_type(intensity_flat.dtype, np.float64))
    if mz_tolerance is not None:
        low = float(mz) - float(mz_tolerance)
        high = float(mz) + float(mz_tolerance)

    offset = 0
    for index, length_value in enumerate(lengths):
        length = int(length_value)
        end = offset + length
        mz_slice = mz_flat[offset:end]
        intensity_slice = intensity_flat[offset:end]
        if length > 0:
            if mz_tolerance is None:
                nearest = int(np.nanargmin(np.abs(mz_slice.astype(float, copy=False) - float(mz))))
                values[index] = intensity_slice[nearest]
            else:
                mask = (mz_slice >= low) & (mz_slice <= high) # pylint: disable=E0606
                values[index] = np.sum(intensity_slice[mask])
        offset = end

    return values


def _infer_image_geometry(
    coordinates: np.ndarray,
    *,
    width: int | None = None,
    height: int | None = None,
) -> tuple[bool, int, int]:
    coords = np.asarray(coordinates)
    xy = coords[:, :2].astype(np.int64, copy=False)
    xs = xy[:, 0]
    ys = xy[:, 1]
    zero_based = bool((xs.size > 0 and int(xs.min()) == 0) or (ys.size > 0 and int(ys.min()) == 0))

    if width is None:
        width = int(xs.max()) + 1 if zero_based else int(xs.max())
    if height is None:
        height = int(ys.max()) + 1 if zero_based else int(ys.max())
    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid image dimensions: width={width}, height={height}")

    return zero_based, width, height


def _scatter_values_into_image(
    image: np.ndarray,
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    zero_based: bool,
) -> None:
    xy = coordinates[:, :2].astype(np.int64, copy=False)
    xs = xy[:, 0]
    ys = xy[:, 1]
    scatter_x = xs if zero_based else xs - 1
    scatter_y = ys if zero_based else ys - 1
    height, width = image.shape[:2]
    valid = (scatter_x >= 0) & (scatter_x < width) & (scatter_y >= 0) & (scatter_y < height)
    image[scatter_y[valid], scatter_x[valid]] = values[valid]


def _imzml_image_geometry(manager: Any) -> tuple[bool, int | None, int | None]:
    meta = getattr(getattr(manager, "ms", None), "meta", None)
    width = getattr(meta, "max_count_of_pixels_x", None)
    height = getattr(meta, "max_count_of_pixels_y", None)
    width = int(width) if width else None
    height = int(height) if height else None
    min_x = min_y = max_x = max_y = None

    for spectrum in manager.ms:
        x = int(spectrum.x)
        y = int(spectrum.y)
        min_x = x if min_x is None else min(min_x, x)
        min_y = y if min_y is None else min(min_y, y)
        max_x = x if max_x is None else max(max_x, x)
        max_y = y if max_y is None else max(max_y, y)

    if min_x is None or min_y is None or max_x is None or max_y is None:
        return False, width, height

    zero_based = bool(min_x == 0 or min_y == 0)
    if width is None:
        width = max_x + 1 if zero_based else max_x
    if height is None:
        height = max_y + 1 if zero_based else max_y
    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid image dimensions: width={width}, height={height}")

    return zero_based, width, height


def _read_imzml_tic_image(
    manager: Any,
    *,
    batch_size: int,
    max_threads: int,
) -> np.ndarray:
    zero_based, width, height = _imzml_image_geometry(manager)
    image: np.ndarray | None = None

    for _mz_data, intensity_flat, lengths, coordinates in manager.flat_generator(
        batch_size=batch_size,
        include_mz=False,
        max_threads=max_threads,
    ):
        intensity_flat = np.asarray(intensity_flat)
        lengths = np.asarray(lengths, dtype=np.int64)
        coordinates = np.asarray(coordinates)
        values = _sum_by_lengths(intensity_flat, lengths)
        if image is None:
            if width is None or height is None:
                zero_based, width, height = _infer_image_geometry(coordinates, width=width, height=height)
            image = np.zeros((height, width), dtype=values.dtype)
        _scatter_values_into_image(image, coordinates, values, zero_based=zero_based)

    if image is None:
        raise ValueError("No pixels were read from imzML source.")

    return image


def _read_imzml_ion_image(
    manager: Any,
    *,
    mz: float | None,
    batch_size: int,
    max_threads: int,
    mz_tolerance: float | None,
) -> np.ndarray:
    if mz is None:
        raise ValueError("'mz' is required.")

    selected_idx: int | None = None
    selected_indices: np.ndarray | None = None
    expected_n_ions: int | None = None
    zero_based, width, height = _imzml_image_geometry(manager)
    image: np.ndarray | None = None

    for mz_data, intensity_flat, lengths, coordinates in manager.flat_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    ):
        mz_arr = None if mz_data is None else np.asarray(mz_data)
        intensity_flat = np.asarray(intensity_flat)
        lengths = np.asarray(lengths, dtype=np.int64)
        coordinates = np.asarray(coordinates)

        if _batch_uses_shared_mz(manager, mz_arr, lengths):
            if mz_arr is None:
                raise ValueError("Shared m/z imzML batch did not include mz_axis.")
            n_ions = int(mz_arr.size)
            if expected_n_ions is None:
                expected_n_ions = n_ions
            elif n_ions != expected_n_ions:
                raise ValueError("All imzML batches must use the same shared mz_axis.")

            if not np.all(lengths == n_ions):
                raise ValueError("Shared m/z imzML batches must have uniform spectrum lengths.")
            intensity_matrix = intensity_flat.reshape(lengths.size, n_ions)
            if mz_tolerance is None:
                if selected_idx is None:
                    selected_idx = _nearest_mz_index(mz_arr, float(mz))
                values = intensity_matrix[:, selected_idx]
            else:
                if selected_indices is None:
                    selected_indices = _mz_window_indices(mz_arr, float(mz), float(mz_tolerance))
                values = np.sum(
                    intensity_matrix[:, selected_indices],
                    axis=1,
                    dtype=np.result_type(intensity_matrix.dtype, np.float64),
                )
        else:
            values = _sum_processed_mz_values(
                mz_arr,
                intensity_flat,
                lengths,
                float(mz),
                None if mz_tolerance is None else float(mz_tolerance),
            )

        values = np.asarray(values)
        if image is None:
            if width is None or height is None:
                zero_based, width, height = _infer_image_geometry(coordinates, width=width, height=height)
            image = np.zeros((height, width), dtype=values.dtype)
        _scatter_values_into_image(image, coordinates, values, zero_based=zero_based)

    if image is None:
        raise ValueError("No pixels were read from imzML source.")

    return image


def _read_ion_manager_ion_image(
    manager: Any,
    *,
    mz: float | None,
    mz_tolerance: float | None,
) -> np.ndarray:
    if mz is None:
        raise ValueError("'mz' is required.")

    mz_axis = np.asarray(manager.mz_axis)
    if mz_tolerance is not None:
        ion_indices = _mz_window_indices(mz_axis, float(mz), float(mz_tolerance))
        return _sum_ion_manager_images(manager, ion_indices)
    ion_idx = _nearest_mz_index(mz_axis, float(mz))
    return np.asarray(manager.get_ion_image(ion_idx))


def _sum_ion_manager_images(manager: Any, ion_indices: np.ndarray) -> np.ndarray:
    height = int(manager.height)
    width = int(manager.width)
    dtype = np.result_type(getattr(manager, "intensity_dtype", np.float32), np.float64)
    image = np.zeros((height, width), dtype=dtype)

    for ion_idx in ion_indices:
        image += np.asarray(manager.get_ion_image(int(ion_idx)), dtype=dtype)

    return image


def _read_ion_manager_tic_image(manager: Any) -> np.ndarray:
    n_ions = int(manager.n_ions)
    if n_ions <= 0:
        raise ValueError("Ion image source contains no ions.")

    if isinstance(manager, (IonImageDataManagerZarr, IonImageDataManagerH5)):
        return _read_dense_ion_manager_tic_image(manager)
    if isinstance(manager, IonImageDataManagerBinary):
        return _read_binary_ion_manager_tic_image(manager)

    ion_image = np.asarray(manager.get_ion_image(0))
    tic_image = np.zeros(ion_image.shape, dtype=np.result_type(ion_image.dtype, np.float64))
    tic_image += ion_image
    for ion_idx in range(1, n_ions):
        tic_image += np.asarray(manager.get_ion_image(ion_idx))
    return tic_image


def _read_dense_ion_manager_tic_image(
    manager: IonImageDataManagerZarr | IonImageDataManagerH5,
) -> np.ndarray:
    n_ions = int(manager.n_ions)
    height = int(manager.height)
    width = int(manager.width)
    dtype = np.result_type(manager.dset.dtype, np.float64)
    tic_image = np.zeros((height, width), dtype=dtype)
    ion_block = max(1, int(getattr(manager, "ion_block", 256)))

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        dense_block = np.asarray(manager.dset[ion_start:ion_end])
        tic_image += np.sum(dense_block, axis=0, dtype=dtype)

    return tic_image


def _read_binary_ion_manager_tic_image(manager: IonImageDataManagerBinary) -> np.ndarray:
    values = np.sum(
        manager.intensity_matrix,
        axis=0,
        dtype=np.result_type(manager.intensity_matrix.dtype, np.float64),
    )
    coords = np.asarray(manager.coordinates)
    zero_based, width, height = _infer_image_geometry(
        coords,
        width=int(manager.width),
        height=int(manager.height),
    )
    image = np.zeros((height, width), dtype=values.dtype)
    _scatter_values_into_image(image, coords, values, zero_based=zero_based)
    return image


def read_ion_image(
    source: str | Path | Any,
    *,
    mz: float | None = None,
    batch_size: int = 256,
    max_threads: int = 2,
    mz_tolerance: float | None = None,
) -> np.ndarray:
    """Read a single ion image from imzML or an ion-image backend."""
    with _resolve_source(source) as (manager, source_format):
        logger.info(f"Reading ion image from {source_format}.")
        if source_format == "imzml":
            return _read_imzml_ion_image(
                manager,
                mz=mz,
                batch_size=int(batch_size),
                max_threads=int(max_threads),
                mz_tolerance=mz_tolerance,
            )
        if source_format in _ION_IMAGE_FORMATS:
            return _read_ion_manager_ion_image(
                manager,
                mz=mz,
                mz_tolerance=mz_tolerance,
            )
    raise ValueError(f"Unsupported source format: {source_format}")


def read_tic_image(
    source: str | Path | Any,
    *,
    batch_size: int = 256,
    max_threads: int = 2,
) -> np.ndarray:
    """Read a TIC image from imzML or an ion-image backend."""
    with _resolve_source(source) as (manager, source_format):
        logger.info(f"Reading TIC image from {source_format}.")
        if source_format == "imzml":
            return _read_imzml_tic_image(
                manager,
                batch_size=int(batch_size),
                max_threads=int(max_threads),
            )
        if source_format in _ION_IMAGE_FORMATS:
            return _read_ion_manager_tic_image(manager)
    raise ValueError(f"Unsupported source format: {source_format}")


def _plot_image(
    image: np.ndarray,
    *,
    save_path: str | Path | None,
    show: bool,
    dpi: int,
    figsize: tuple[float, float],
    title: str | None,
    structure_only: bool = False,
) -> None:
    render_figsize = figsize
    if structure_only and image.ndim >= 2 and image.shape[0] > 0:
        height, width = image.shape[:2]
        render_figsize = (float(figsize[1]) * float(width) / float(height), float(figsize[1]))

    fig, ax = plt.subplots(1, 1, figsize=render_figsize, dpi=int(dpi))
    vmax = float(np.nanpercentile(image, 95)) if image.size else None
    if vmax is not None and not np.isfinite(vmax):
        vmax = None

    if structure_only:
        ax.imshow(image, cmap="inferno", vmax=vmax)
        ax.set_axis_off()
        ax.set_position((0, 0, 1, 1))
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    else:
        mappable = ax.imshow(image, cmap="inferno", vmax=vmax)
        if title:
            ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(mappable, ax=ax)
        fig.tight_layout()

    if save_path:
        output_path = Path(save_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if structure_only:
            fig.savefig(output_path, dpi=int(dpi), bbox_inches=None, pad_inches=0)
        else:
            fig.savefig(output_path, dpi=int(dpi), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def plot_ion_image(
    source: str | Path | Any,
    image_matrix: np.ndarray | None = None,
    *,
    mz: float | None = None,
    batch_size: int = 256,
    max_threads: int = 2,
    mz_tolerance: float | None = None,
    save_path: str | Path | None = None,
    show: bool = True,
    dpi: int = 300,
    figsize: tuple[float, float] = (20, 5),
    title: str | None = None,
) -> None:
    """Read and plot one ion image."""
    if image_matrix is None:
        image_matrix = read_ion_image(
            source,
            mz=mz,
            batch_size=batch_size,
            max_threads=max_threads,
            mz_tolerance=mz_tolerance,
        )
    _plot_image(
        image_matrix,
        save_path=save_path,
        show=bool(show),
        dpi=int(dpi),
        figsize=(float(figsize[0]), float(figsize[1])),
        title=title,
    )


def plot_tic_image(
    source: str | Path | Any,
    image_matrix: np.ndarray | None = None,
    *,
    batch_size: int = 256,
    max_threads: int = 2,
    save_path: str | Path | None = None,
    show: bool = True,
    dpi: int = 300,
    figsize: tuple[float, float] = (20, 5),
    title: str | None = None,
    structure_only: bool = False,
) -> None:
    """Read and plot a TIC image."""
    if image_matrix is None:
        image_matrix = read_tic_image(
            source,
            batch_size=batch_size,
            max_threads=max_threads,
        )
    _plot_image(
        image_matrix,
        save_path=save_path,
        show=bool(show),
        dpi=int(dpi),
        figsize=(float(figsize[0]), float(figsize[1])),
        title=title,
        structure_only=structure_only,
    )
