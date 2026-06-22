"""Apply a spatial region mask to MSI data managers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
)

DataManager = (
    MSDataManagerImzML
    | IonImageDataManagerBinary
    | IonImageDataManagerH5
    | IonImageDataManagerZarr
)
IonImageDataManager = (
    IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr
)


# ---------------------------------------------------------------------------
# Spatial mask helpers
# ---------------------------------------------------------------------------


def _spatial_dimensions(data_manager: DataManager) -> tuple[int, int]:
    """Return the source image width and height."""
    if isinstance(data_manager, MSDataManagerImzML):
        meta = data_manager.ms.meta
        return int(meta.max_count_of_pixels_x), int(meta.max_count_of_pixels_y)  # type: ignore
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


# ---------------------------------------------------------------------------
# imzML backend
# ---------------------------------------------------------------------------


def _apply_imzml_mask(
    data_manager: MSDataManagerImzML,
    mask: np.ndarray,
    output_path: str | Path,
    *,
    batch_size: int,
    temp_dir: str | Path | None,
) -> MSDataManagerImzML:
    """Write selected imzML spectra to a new imzML manager."""
    output = MSDataManagerImzML(
        filepath=str(output_path),
        max_threads=data_manager.max_threads,
        temp_dir=None if temp_dir is None else str(temp_dir),
        mz_dtype=data_manager.mz_dtype,
        intensity_dtype=data_manager.intensity_dtype,
    )
    output.copy_meta(data_manager)
    output.ms.meta.continuous = data_manager.ms.meta.continuous
    output.ms.meta.processed = data_manager.ms.meta.processed
    continuous = bool(data_manager.ms.meta.continuous)
    selected_pixels = 0

    for mz_data, intensity_flat, lengths, coordinates in data_manager.flat_generator(
        batch_size=batch_size,
        include_mz=True,
    ):
        selected = _select_coordinates(mask, coordinates)
        if not bool(np.any(selected)):
            continue

        selected_lengths = np.asarray(lengths)[selected]
        selected_intensity = _filter_flat_batch(intensity_flat, lengths, selected)
        selected_mz = (
            mz_data
            if continuous or mz_data is None
            else _filter_flat_batch(mz_data, lengths, selected)
        )
        selected_coordinates = np.asarray(coordinates)[selected]
        output.write_flat_batch(
            intensity_flat=selected_intensity,
            lengths=selected_lengths,
            coordinates=selected_coordinates,
            mz_axis=selected_mz,
        )
        selected_pixels += int(selected_lengths.size)

    if selected_pixels == 0:
        output.close()
        raise ValueError("mask does not retain any source pixels.")

    output.close_writer()
    output.load_head_data()
    return output


# ---------------------------------------------------------------------------
# Ion-image backend helpers
# ---------------------------------------------------------------------------


def _selected_coordinates(
    data_manager: IonImageDataManager,
    mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return selected coordinate rows and their source indices."""
    selected = _select_coordinates(
        mask,
        data_manager.coordinates,
        (
            "zero"
            if getattr(data_manager, "coordinate_base", None) == 0
            else None
        ),
    )
    selected_indices = np.flatnonzero(selected)
    if selected_indices.size == 0:
        raise ValueError("mask does not retain any source pixels.")
    return np.asarray(data_manager.coordinates)[selected_indices], selected_indices


def _apply_ion_image_mask(
    data_manager: IonImageDataManager,
    mask: np.ndarray,
    output_path: str | Path,
    *,
    batch_size: int,
    temp_dir: str | Path | None,
) -> IonImageDataManager:
    """Copy selected pixels through the common ion-image writer contract."""
    coordinates, selected_indices = _selected_coordinates(data_manager, mask)
    mz_axis = np.asarray(data_manager.mz_axis)
    if isinstance(data_manager, IonImageDataManagerZarr):
        output: IonImageDataManager = IonImageDataManagerZarr(
            output_path,
            mode="w",
            width=data_manager.width,
            height=data_manager.height,
            temp_dir=temp_dir,
            intensity_dtype=data_manager.intensity_dtype,
            metadata=data_manager.meta,
            shared_mz_axis=mz_axis,
        )
    elif isinstance(data_manager, IonImageDataManagerH5):
        output = IonImageDataManagerH5(
            output_path,
            mode="w",
            width=data_manager.width,
            height=data_manager.height,
            temp_dir=temp_dir,
            mz_dtype=data_manager.mz_dtype,
            intensity_dtype=data_manager.intensity_dtype,
            shared_mz_axis=mz_axis,
        )
    else:
        output = IonImageDataManagerBinary(
            output_path,
            mode="w",
            width=data_manager.width,
            height=data_manager.height,
            temp_dir=temp_dir,
            mz_dtype=data_manager.mz_dtype,
            intensity_dtype=data_manager.intensity_dtype,
            shared_mz_axis=mz_axis,
        )

    n_ions = int(data_manager.n_ions)
    ion_block = max(1, int(data_manager.ion_block))
    try:
        for pixel_start in range(0, selected_indices.size, batch_size):
            pixel_end = min(selected_indices.size, pixel_start + batch_size)
            pixel_indices = selected_indices[pixel_start:pixel_end]
            matrix = np.empty(
                (pixel_indices.size, n_ions),
                dtype=data_manager.intensity_dtype,
            )
            for ion_start in range(0, n_ions, ion_block):
                ion_end = min(n_ions, ion_start + ion_block)
                matrix[:, ion_start:ion_end] = data_manager.get_ion_block(
                    ion_start,
                    ion_end,
                )[:, pixel_indices].T
            output.write_flat_batch(
                intensity_flat=matrix.reshape(-1),
                lengths=np.full(pixel_indices.size, n_ions, dtype="<u4"),
                coordinates=coordinates[pixel_start:pixel_end],
                mz_axis=mz_axis,
            )
        output.finalize()
        return output
    except BaseException:
        output.close()
        raise


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def apply_region_mask(
    data_manager: DataManager,
    mask: np.ndarray,
    output_path: str | Path,
    *,
    batch_size: int = 256,
    temp_dir: str | Path | None = None,
) -> DataManager:
    """Return a same-type manager containing only mask-valid pixels.

    ``batch_size`` controls spectrum batches. ``temp_dir`` is used by output
    writers that stage intermediate data.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    width, height = _spatial_dimensions(data_manager)
    normalized_mask = _normalize_mask(mask, width=width, height=height)

    if isinstance(data_manager, MSDataManagerImzML):
        return _apply_imzml_mask(
            data_manager,
            normalized_mask,
            output_path,
            batch_size=batch_size,
            temp_dir=temp_dir,
        )

    return _apply_ion_image_mask(
        data_manager,
        normalized_mask,
        output_path,
        batch_size=batch_size,
        temp_dir=temp_dir,
    )
