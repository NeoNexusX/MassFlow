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
    is_loaded = len(manager.ms) > 0

    if not is_loaded:
        manager.load_head_data()


def _validate_ion_selector(mz: float | None, ion_idx: int | None) -> None:
    if (mz is None) == (ion_idx is None):
        raise ValueError("Exactly one of 'mz' or 'ion_idx' must be provided.")


def _metadata_dimensions(manager: Any) -> tuple[int | None, int | None]:
    meta = getattr(getattr(manager, "ms", None), "meta", None)
    width = getattr(meta, "max_count_of_pixels_x", None)
    height = getattr(meta, "max_count_of_pixels_y", None)
    width = int(width) if width else None
    height = int(height) if height else None
    return width, height


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
    mz_tolerance: float | None = None,
) -> int:
    if mz_axis.ndim != 1 or mz_axis.size == 0:
        raise ValueError("mz_axis must be a non-empty 1-D array.")
    ion_idx = int(np.nanargmin(np.abs(mz_axis.astype(float, copy=False) - float(mz))))
    if mz_tolerance is not None:
        delta = abs(float(mz_axis[ion_idx]) - float(mz))
        if delta > float(mz_tolerance):
            raise ValueError(
                f"No m/z value within mz_tolerance={mz_tolerance}; "
                f"nearest={float(mz_axis[ion_idx])}, requested={mz}."
            )
    return ion_idx


def _sum_by_lengths(intensity_flat: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    values = np.zeros(lengths.size, dtype=np.result_type(intensity_flat.dtype, np.float64))
    offset = 0
    for index, length_value in enumerate(lengths):
        length = int(length_value)
        end = offset + length
        values[index] = np.sum(intensity_flat[offset:end])
        offset = end
    return values


def _sum_processed_mz_window(
    mz_flat: np.ndarray | None,
    intensity_flat: np.ndarray,
    lengths: np.ndarray,
    mz: float,
    mz_tolerance: float,
) -> np.ndarray:
    if mz_flat is None:
        raise ValueError("Processed imzML ion image requires m/z data.")
    if mz_flat.size != int(np.sum(lengths, dtype=np.int64)):
        raise ValueError("Processed imzML ion image requires flattened per-spectrum m/z data.")

    low = float(mz) - float(mz_tolerance)
    high = float(mz) + float(mz_tolerance)
    values = np.zeros(lengths.size, dtype=np.result_type(intensity_flat.dtype, np.float64))

    offset = 0
    for index, length_value in enumerate(lengths):
        length = int(length_value)
        end = offset + length
        mz_slice = mz_flat[offset:end]
        intensity_slice = intensity_flat[offset:end]
        mask = (mz_slice >= low) & (mz_slice <= high)
        values[index] = np.sum(intensity_slice[mask])
        offset = end

    return values


def _dense_image_from_points(
    coordinates: np.ndarray,
    values: np.ndarray,
    *,
    width: int | None = None,
    height: int | None = None,
) -> np.ndarray:
    coords = np.asarray(coordinates)
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("coordinates must have shape (n_pixels, >=2).")
    values = np.asarray(values)
    if coords.shape[0] != values.shape[0]:
        raise ValueError("coordinates and values must have the same row count.")
    if coords.shape[0] == 0:
        raise ValueError("No pixels were read from the source.")

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

    scatter_x = xs if zero_based else xs - 1
    scatter_y = ys if zero_based else ys - 1
    image = np.zeros((height, width), dtype=values.dtype)
    valid = (scatter_x >= 0) & (scatter_x < width) & (scatter_y >= 0) & (scatter_y < height)
    image[scatter_y[valid], scatter_x[valid]] = values[valid]
    return image


def _read_imzml_tic_image(
    manager: Any,
    *,
    batch_size: int,
    max_threads: int,
) -> np.ndarray:
    coord_parts: list[np.ndarray] = []
    value_parts: list[np.ndarray] = []

    for _mz_data, intensity_flat, lengths, coordinates in manager.flat_generator(
        batch_size=batch_size,
        include_mz=False,
        max_threads=max_threads,
    ):
        intensity_flat = np.asarray(intensity_flat)
        lengths = np.asarray(lengths, dtype=np.int64)
        coord_parts.append(np.asarray(coordinates))
        value_parts.append(_sum_by_lengths(intensity_flat, lengths))

    if not coord_parts:
        raise ValueError("No pixels were read from imzML source.")

    width, height = _metadata_dimensions(manager)
    return _dense_image_from_points(
        np.vstack(coord_parts),
        np.concatenate(value_parts),
        width=width,
        height=height,
    )


def _read_imzml_ion_image(
    manager: Any,
    *,
    mz: float | None,
    ion_idx: int | None,
    batch_size: int,
    max_threads: int,
    mz_tolerance: float | None,
) -> np.ndarray:
    _validate_ion_selector(mz, ion_idx)
    if mz_tolerance is not None and mz is None:
        raise ValueError("'mz_tolerance' can only be used with 'mz'.")

    selected_idx = None if ion_idx is None else int(ion_idx)
    expected_n_ions: int | None = None
    coord_parts: list[np.ndarray] = []
    value_parts: list[np.ndarray] = []

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

            if selected_idx is None:
                if mz is None:
                    raise ValueError("'mz' is required when 'ion_idx' is not provided.")
                selected_idx = _nearest_mz_index(mz_arr, float(mz), mz_tolerance)
            if selected_idx < 0 or selected_idx >= n_ions:
                raise IndexError(f"ion_idx out of range: {selected_idx}")
            if not np.all(lengths == n_ions):
                raise ValueError("Shared m/z imzML batches must have uniform spectrum lengths.")

            values = intensity_flat.reshape(lengths.size, n_ions)[:, selected_idx]
        else:
            if mz is None:
                raise ValueError(
                    "Variable-length imzML ion image requires 'mz' and 'mz_tolerance'; "
                    "'ion_idx' only works for shared m/z imzML data."
                )
            if mz_tolerance is None:
                raise ValueError("Processed imzML ion image requires 'mz_tolerance'.")
            values = _sum_processed_mz_window(
                mz_arr,
                intensity_flat,
                lengths,
                float(mz),
                float(mz_tolerance),
            )

        coord_parts.append(coordinates)
        value_parts.append(np.asarray(values))

    if not coord_parts:
        raise ValueError("No pixels were read from imzML source.")

    width, height = _metadata_dimensions(manager)
    return _dense_image_from_points(
        np.vstack(coord_parts),
        np.concatenate(value_parts),
        width=width,
        height=height,
    )


def _read_ion_manager_ion_image(
    manager: Any,
    *,
    mz: float | None,
    ion_idx: int | None,
    mz_tolerance: float | None,
) -> np.ndarray:
    _validate_ion_selector(mz, ion_idx)
    if ion_idx is None:
        mz_axis = np.asarray(manager.mz_axis)
        ion_idx = _nearest_mz_index(mz_axis, float(mz), mz_tolerance)  # type: ignore[arg-type]
    else:
        ion_idx = int(ion_idx)
        n_ions = int(manager.n_ions)
        if ion_idx < 0 or ion_idx >= n_ions:
            raise IndexError(f"ion_idx out of range: {ion_idx}")

    return np.asarray(manager.get_ion_image(ion_idx))


def _read_ion_manager_tic_image(manager: Any) -> np.ndarray:
    n_ions = int(manager.n_ions)
    if n_ions <= 0:
        raise ValueError("Ion image source contains no ions.")

    tic_image: np.ndarray | None = None
    for ion_idx in range(n_ions):
        ion_image = np.asarray(manager.get_ion_image(ion_idx))
        if tic_image is None:
            tic_image = np.zeros(ion_image.shape, dtype=np.result_type(ion_image.dtype, np.float64))
        tic_image += ion_image

    if tic_image is None:
        raise ValueError("Ion image source contains no ions.")
    return tic_image


def read_ion_image(
    source: str | Path | Any,
    *,
    mz: float | None = None,
    ion_idx: int | None = None,
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
                ion_idx=ion_idx,
                batch_size=int(batch_size),
                max_threads=int(max_threads),
                mz_tolerance=mz_tolerance,
            )
        if source_format in _ION_IMAGE_FORMATS:
            return _read_ion_manager_ion_image(
                manager,
                mz=mz,
                ion_idx=ion_idx,
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
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=int(dpi))
    vmax = float(np.nanpercentile(image, 95)) if image.size else None
    if vmax is not None and not np.isfinite(vmax):
        vmax = None

    if structure_only:
        ax.imshow(image, cmap="inferno", vmax=vmax, aspect="auto")
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
    *,
    mz: float | None = None,
    ion_idx: int | None = None,
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
    image = read_ion_image(
        source,
        mz=mz,
        ion_idx=ion_idx,
        batch_size=batch_size,
        max_threads=max_threads,
        mz_tolerance=mz_tolerance,
    )
    _plot_image(
        image,
        save_path=save_path,
        show=bool(show),
        dpi=int(dpi),
        figsize=(float(figsize[0]), float(figsize[1])),
        title=title,
    )


def plot_tic_image(
    source: str | Path | Any,
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
    image = read_tic_image(
        source,
        batch_size=batch_size,
        max_threads=max_threads,
    )
    _plot_image(
        image,
        save_path=save_path,
        show=bool(show),
        dpi=int(dpi),
        figsize=(float(figsize[0]), float(figsize[1])),
        title=title,
        structure_only=structure_only,
    )
