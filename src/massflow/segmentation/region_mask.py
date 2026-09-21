"""Apply a spatial region mask to MSI data managers."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import h5py
import hdf5plugin
import numpy as np
import zarr

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
)

from massflow.data_manager.ion_image_data_manager_binary import (
    HEADER_SIZE,
    UINT32_SIZE,
    IonImageBinaryHeader,
)

from massflow.msi_zarr.schema import (
    ZARR_FORMAT_NAME,
    ZARR_FORMAT_VERSION,
    ZARR_GROUP_ION_IMAGE,
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
        getattr(data_manager, "_coordinate_base", None),
    )
    selected_indices = np.flatnonzero(selected)
    if selected_indices.size == 0:
        raise ValueError("mask does not retain any source pixels.")
    return np.asarray(data_manager.coordinates)[selected_indices], selected_indices


# --- Binary ion-image backend ----------------------------------------------


def _write_binary_masked(
    data_manager: IonImageDataManagerBinary,
    mask: np.ndarray,
    output_path: str | Path,
    *,
    batch_size: int,
) -> IonImageDataManagerBinary:
    """Write selected Binary pixel columns directly to a new Binary file."""
    coordinates, selected_indices = _selected_coordinates(data_manager, mask)
    mz_axis = np.asarray(data_manager.mz_axis, dtype=data_manager.mz_dtype)
    coordinates = coordinates[:, :2].astype("<u4", copy=False)
    n_pixels = int(selected_indices.size)
    mz_offset = HEADER_SIZE
    coordinate_offset = mz_offset + data_manager.n_ions * data_manager.mz_dtype.itemsize
    intensity_offset = coordinate_offset + n_pixels * 2 * UINT32_SIZE
    header = IonImageBinaryHeader(
        width=data_manager.width,
        height=data_manager.height,
        n_ions=data_manager.n_ions,
        n_pixels=n_pixels,
        mz_dtype_bits=data_manager.mz_dtype.itemsize * 8,
        intensity_dtype_bits=data_manager.intensity_dtype.itemsize * 8,
        mz_offset=mz_offset,
        coordinate_offset=coordinate_offset,
        intensity_offset=intensity_offset,
    )

    with Path(output_path).open("wb") as output_file:
        output_file.write(header.to_bytes())
        output_file.write(memoryview(np.ascontiguousarray(mz_axis)))
        output_file.write(memoryview(np.ascontiguousarray(coordinates)))
        for start in range(0, data_manager.n_ions, batch_size):
            end = min(data_manager.n_ions, start + batch_size)
            block = np.ascontiguousarray(
                data_manager.intensity_matrix[start:end, selected_indices]
            )
            output_file.write(memoryview(block))

    return IonImageDataManagerBinary(output_path, mode="r")


# --- Dense ion-image block masking (H5 / Zarr) -----------------------------


def _masked_dense_block(
    data_manager: IonImageDataManagerH5,
    mask: np.ndarray,
    start: int,
    end: int,
) -> np.ndarray:
    """Read a dense ion-image block and zero pixels outside the mask."""
    block = np.asarray(data_manager._dset[start:end])  # type: ignore  # pylint: disable=protected-access
    block *= mask
    return block


# --- H5 ion-image backend ---------------------------------------------------


def _write_h5_masked(
    data_manager: IonImageDataManagerH5,
    mask: np.ndarray,
    output_path: str | Path,
    *,
    batch_size: int,
) -> IonImageDataManagerH5:
    """Write masked dense ion-image blocks directly to a new H5 file."""
    coordinates, _ = _selected_coordinates(data_manager, mask)
    chunk_size = int(data_manager._dset.chunks[0])  # type: ignore  # pylint: disable=protected-access
    compression = hdf5plugin.Bitshuffle(nelems=0, cname="lz4", clevel=1)  # type: ignore

    with h5py.File(output_path, "w") as output_file:
        output_file.attrs["width"] = data_manager.width
        output_file.attrs["height"] = data_manager.height
        output_file.attrs["n_ions"] = data_manager.n_ions
        output_file.attrs["n_pixels"] = coordinates.shape[0]
        output_file.create_dataset("mz_axis", data=np.asarray(data_manager.mz_axis))
        output_file.create_dataset("coordinates", data=coordinates, dtype="<u4")
        output_images = output_file.create_dataset(
            "ion_images",
            shape=(data_manager.n_ions, data_manager.height, data_manager.width),
            chunks=(chunk_size, data_manager.height, data_manager.width),
            dtype=data_manager.intensity_dtype,
            **compression,
        )
        for start in range(0, data_manager.n_ions, batch_size):
            end = min(data_manager.n_ions, start + batch_size)
            output_images[start:end] = _masked_dense_block(
                data_manager,
                mask,
                start,
                end,
            )

    return IonImageDataManagerH5(output_path, mode="r")


# --- Zarr ion-image backend -------------------------------------------------


def _write_zarr_masked(
    data_manager: IonImageDataManagerZarr,
    mask: np.ndarray,
    output_path: str | Path,
    *,
    batch_size: int,
) -> IonImageDataManagerZarr:
    """Write selected compact ion vectors to a new Zarr group."""
    output_path = Path(output_path)
    if output_path.suffix != ".zarr":
        output_path = Path(f"{output_path}.zarr")
    coordinates, selected_indices = _selected_coordinates(data_manager, mask)
    n_pixels = int(selected_indices.size)
    reader = data_manager._reader  # pylint: disable=protected-access
    assert reader is not None
    chunk_points = int(reader.intensity.chunks[0])
    output_group = zarr.open_group(str(output_path), mode="w")
    output_group.attrs.update(
        {
            "format": ZARR_FORMAT_NAME,
            "format_version": ZARR_FORMAT_VERSION,
            "write_state": "writing",
            "coordinate_order": ["x", "y", "z"],
            "coordinate_base": 1,
            "spatial_shape": [data_manager.height, data_manager.width],
        }
    )
    axes_group = output_group.create_group("axes")
    axes_group.create_array(
        "coordinates",
        data=coordinates.astype("<u4", copy=False),
    )
    mz_axis = np.asarray(data_manager.mz_axis, dtype=data_manager.mz_dtype)
    mz_array = axes_group.create_array("mz", data=mz_axis)
    mz_array.attrs["unit"] = "m/z"
    data_group = output_group.create_group(ZARR_GROUP_ION_IMAGE)
    data_group.attrs.update({"row_axis": "ion", "encoding": "continuous"})
    output_intensity = data_group.create_array(
        "intensity",
        shape=(data_manager.n_ions * n_pixels,),
        chunks=(max(n_pixels, chunk_points),),
        dtype=data_manager.intensity_dtype,
    )
    data_group.create_array(
        "offsets",
        data=(np.arange(data_manager.n_ions + 1, dtype=np.int64) * n_pixels),
    )
    # ion-only 输出的 mean_spectrum 基于筛选后的有效像素重新计算。
    mean_spectrum = np.empty(data_manager.n_ions, dtype=np.float64)
    for start in range(0, data_manager.n_ions, batch_size):
        end = min(data_manager.n_ions, start + batch_size)
        block = data_manager.get_ion_block(start, end)[:, selected_indices]
        output_intensity[start * n_pixels : end * n_pixels] = block.reshape(-1)
        mean_spectrum[start:end] = block.sum(axis=1, dtype=np.float64) / n_pixels

    stats = output_group.create_group("stats")
    mean_array = stats.create_array("mean_spectrum", data=mean_spectrum)
    mean_array.attrs.update(
        {
            "source": f"/{ZARR_GROUP_ION_IMAGE}",
            "reduced_axis": "pixel",
            "denominator": "sampled_pixels",
        }
    )

    meta = deepcopy(data_manager.meta)
    meta.spectrum_count_num = n_pixels
    meta.mask = np.asarray(mask, dtype=bool)
    meta.continuous = True
    meta.processed = None
    meta.to_zarr(output_group)
    output_group.attrs["write_state"] = "complete"
    return IonImageDataManagerZarr(output_path, mode="r")


# --- Ion-image backend dispatch --------------------------------------------


def _apply_ion_image_mask(
    data_manager: IonImageDataManager,
    mask: np.ndarray,
    output_path: str | Path,
    *,
    batch_size: int,
) -> IonImageDataManager:
    """Apply a mask directly to ion-major blocks and return the same backend type."""
    if isinstance(data_manager, IonImageDataManagerZarr):
        return _write_zarr_masked(
            data_manager,
            mask,
            output_path,
            batch_size=batch_size,
        )
    if isinstance(data_manager, IonImageDataManagerH5):
        return _write_h5_masked(
            data_manager,
            mask,
            output_path,
            batch_size=batch_size,
        )
    return _write_binary_masked(
        data_manager,
        mask,
        output_path,
        batch_size=batch_size,
    )


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

    ``batch_size`` controls spectrum batches for imzML and ion-image blocks for
    Binary/H5/Zarr. ``temp_dir`` is only used by the imzML writer.
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
    )
