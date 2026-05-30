"""Data loading helpers for feature analysis."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

from massflow.data_manager import IonImageDataManagerBinary
from massflow.data_manager import IonImageDataManagerH5
from massflow.data_manager import IonImageDataManagerZarr


def _read_binary_pixel_columns(
    manager: IonImageDataManagerBinary,
    pixel_indexer: slice | np.ndarray,
    out: np.ndarray,
    ion_block: int,
) -> None:
    """Fill *out* ``(n_pixels_subset, n_ions)`` from a binary ion-image backend."""
    n_ions = int(manager.n_ions)
    intensity_matrix = manager.intensity_matrix

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        ion_slice = intensity_matrix[ion_start:ion_end, pixel_indexer]
        out[:, ion_start:ion_end] = np.asarray(ion_slice).T


def read_binary_sample_matrix(
    manager: IonImageDataManagerBinary,
    sample_ratio: float,
    *,
    random_state: int | None = None,
    ion_block: int = 256,
) -> np.ndarray:
    """Randomly sample nonzero spectra from a binary ion-image backend."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    sample_count = max(1, int(np.ceil(n_pixels * sample_ratio)))

    rng = np.random.default_rng(random_state)
    pixel_indices = np.sort(rng.choice(n_pixels, size=sample_count, replace=False))

    sample_matrix = np.empty((sample_count, n_ions), dtype=np.float32)
    _read_binary_pixel_columns(manager, pixel_indices, sample_matrix, ion_block)

    nonzero_mask = np.any(sample_matrix != 0, axis=1)
    sample_matrix = sample_matrix[nonzero_mask]

    if sample_matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")

    return sample_matrix


def iter_binary_matrix_batches(
    manager: IonImageDataManagerBinary,
    *,
    pixel_batch_size: int = 20000,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield ``(batch_matrix, batch_coordinates)`` by pixel-batch from binary data."""
    if pixel_batch_size <= 0:
        raise ValueError("pixel_batch_size must be positive.")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    coordinates = np.asarray(manager.coordinates)

    for pixel_start in range(0, n_pixels, pixel_batch_size):
        pixel_end = min(n_pixels, pixel_start + pixel_batch_size)
        pixel_slice = slice(pixel_start, pixel_end)
        batch_coords = coordinates[pixel_slice, :2]

        batch_matrix = np.empty((pixel_end - pixel_start, n_ions), dtype=np.float32)
        _read_binary_pixel_columns(manager, pixel_slice, batch_matrix, ion_block)

        nonzero_mask = np.any(batch_matrix != 0, axis=1)
        if not np.any(nonzero_mask):
            continue

        if np.all(nonzero_mask):
            yield batch_matrix, batch_coords
        else:
            yield batch_matrix[nonzero_mask], batch_coords[nonzero_mask]


def _read_zarr_pixel_columns(
    manager: IonImageDataManagerZarr,
    pixel_flat_idx: np.ndarray,
    out: np.ndarray,
    ion_block: int,
) -> None:
    """Fill *out* ``(n_pixels_subset, n_ions)`` from a Zarr ion-image backend."""
    n_ions = int(manager.n_ions)
    height = int(manager.height)
    width = int(manager.width)
    image_size = height * width

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        dense_block = np.asarray(manager.dset[ion_start:ion_end])
        flat_block = dense_block.reshape(ion_end - ion_start, image_size)
        out[:, ion_start:ion_end] = flat_block[:, pixel_flat_idx].T


def read_zarr_sample_matrix(
    manager: IonImageDataManagerZarr,
    sample_ratio: float,
    *,
    random_state: int | None = None,
    ion_block: int = 256,
) -> np.ndarray:
    """Randomly sample nonzero spectra from a Zarr ion-image backend."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    flat_idx = manager.flat_idx
    sample_count = max(1, int(np.ceil(n_pixels * sample_ratio)))

    rng = np.random.default_rng(random_state)
    pixel_indices = np.sort(rng.choice(n_pixels, size=sample_count, replace=False))
    sample_flat_idx = flat_idx[pixel_indices]

    sample_matrix = np.empty((sample_count, n_ions), dtype=np.float32)
    _read_zarr_pixel_columns(manager, sample_flat_idx, sample_matrix, ion_block)

    nonzero_mask = np.any(sample_matrix != 0, axis=1)
    sample_matrix = sample_matrix[nonzero_mask]

    if sample_matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")

    return sample_matrix


def iter_zarr_matrix_batches(
    manager: IonImageDataManagerZarr,
    *,
    pixel_batch_size: int = 20000,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield ``(batch_matrix, batch_coordinates)`` by pixel-batch from Zarr data."""
    if pixel_batch_size <= 0:
        raise ValueError("pixel_batch_size must be positive.")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    flat_idx = manager.flat_idx
    coordinates = np.asarray(manager.coordinates)

    for pixel_start in range(0, n_pixels, pixel_batch_size):
        pixel_end = min(n_pixels, pixel_start + pixel_batch_size)
        batch_flat_idx = flat_idx[pixel_start:pixel_end]
        batch_coords = coordinates[pixel_start:pixel_end, :2]

        batch_matrix = np.empty((pixel_end - pixel_start, n_ions), dtype=np.float32)
        _read_zarr_pixel_columns(manager, batch_flat_idx, batch_matrix, ion_block)

        nonzero_mask = np.any(batch_matrix != 0, axis=1)
        if not np.any(nonzero_mask):
            continue

        if np.all(nonzero_mask):
            yield batch_matrix, batch_coords
        else:
            yield batch_matrix[nonzero_mask], batch_coords[nonzero_mask]


def _read_h5_pixel_columns(
    manager: IonImageDataManagerH5,
    pixel_flat_idx: np.ndarray,
    out: np.ndarray,
    ion_block: int,
) -> None:
    """Fill *out* ``(n_pixels_subset, n_ions)`` from an H5 ion-image backend."""
    n_ions = int(manager.n_ions)
    height = int(manager.height)
    width = int(manager.width)
    image_size = height * width

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        dense_block = np.asarray(manager.dset[ion_start:ion_end])
        flat_block = dense_block.reshape(ion_end - ion_start, image_size)
        out[:, ion_start:ion_end] = flat_block[:, pixel_flat_idx].T


def read_h5_sample_matrix(
    manager: IonImageDataManagerH5,
    sample_ratio: float,
    *,
    random_state: int | None = None,
    ion_block: int = 256,
) -> np.ndarray:
    """Randomly sample nonzero spectra from an H5 ion-image backend."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    flat_idx = manager.flat_idx
    sample_count = max(1, int(np.ceil(n_pixels * sample_ratio)))

    rng = np.random.default_rng(random_state)
    pixel_indices = np.sort(rng.choice(n_pixels, size=sample_count, replace=False))
    sample_flat_idx = flat_idx[pixel_indices]

    sample_matrix = np.empty((sample_count, n_ions), dtype=np.float32)
    _read_h5_pixel_columns(manager, sample_flat_idx, sample_matrix, ion_block)

    nonzero_mask = np.any(sample_matrix != 0, axis=1)
    sample_matrix = sample_matrix[nonzero_mask]

    if sample_matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")

    return sample_matrix


def iter_h5_matrix_batches(
    manager: IonImageDataManagerH5,
    *,
    pixel_batch_size: int = 20000,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield ``(batch_matrix, batch_coordinates)`` by pixel-batch from H5 data."""
    if pixel_batch_size <= 0:
        raise ValueError("pixel_batch_size must be positive.")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    flat_idx = manager.flat_idx
    coordinates = np.asarray(manager.coordinates)

    for pixel_start in range(0, n_pixels, pixel_batch_size):
        pixel_end = min(n_pixels, pixel_start + pixel_batch_size)
        batch_flat_idx = flat_idx[pixel_start:pixel_end]
        batch_coords = coordinates[pixel_start:pixel_end, :2]

        batch_matrix = np.empty((pixel_end - pixel_start, n_ions), dtype=np.float32)
        _read_h5_pixel_columns(manager, batch_flat_idx, batch_matrix, ion_block)

        nonzero_mask = np.any(batch_matrix != 0, axis=1)
        if not np.any(nonzero_mask):
            continue

        if np.all(nonzero_mask):
            yield batch_matrix, batch_coords
        else:
            yield batch_matrix[nonzero_mask], batch_coords[nonzero_mask]


__all__ = [
    "iter_binary_matrix_batches",
    "iter_h5_matrix_batches",
    "iter_zarr_matrix_batches",
    "read_binary_sample_matrix",
    "read_h5_sample_matrix",
    "read_zarr_sample_matrix",
]
