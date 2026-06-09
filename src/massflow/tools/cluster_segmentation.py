"""Segmentation tools for clustering data."""

from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import zarr
from scipy.ndimage import binary_propagation


def _read_array(result_group: zarr.Group, name: str) -> np.ndarray:
    """Read a required Zarr array as a NumPy array."""
    if name not in result_group:
        raise ValueError(f"Zarr group must contain '{name}'.")
    array = result_group[name]
    if not isinstance(array, zarr.Array):
        raise ValueError(f"Zarr group item '{name}' must be an array.")
    return np.asarray(array)


def _save_transparent_image(
    image: np.ndarray,
    mask: np.ndarray,
    output_path: Path,
) -> None:
    """Apply mask as alpha and save an RGB image as transparent PNG."""
    if image.shape != (*mask.shape, 3):
        raise ValueError("Feature-analysis image must have shape (height, width, 3).")

    rgba = np.empty((*mask.shape, 4), dtype=np.uint8)
    rgba[..., :3] = image
    rgba[..., 3] = mask * 255
    plt.imsave(str(output_path), rgba)


def _read_feature_analysis(
    zarr_path: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read label, KMeans, and UMAP images from a feature-analysis Zarr group."""
    result_group = zarr.open_group(str(zarr_path), mode="r")
    label_image = _read_array(result_group, "kmeans_label_image")
    kmeans_image = _read_array(result_group, "kmeans_image")
    umap_image = _read_array(result_group, "umap_image")
    if label_image.ndim != 2:
        raise ValueError("kmeans_label_image must be a two-dimensional array.")
    return label_image, kmeans_image, umap_image


def _save_outputs(
    mask: np.ndarray,
    kmeans_image: np.ndarray,
    umap_image: np.ndarray,
    output_dir: str | Path,
) -> None:
    """Save transparent UMAP and KMeans images."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    _save_transparent_image(umap_image, mask, out / "umap_image.png")
    _save_transparent_image(kmeans_image, mask, out / "kmeans_image.png")


def _normalize_labels(labels: Sequence[int] | None, *, parameter_name: str) -> np.ndarray | None:
    """Validate optional KMeans labels."""
    if labels is None:
        return None
    normalized = np.asarray(labels, dtype=np.int64)
    if normalized.ndim != 1:
        raise ValueError(f"{parameter_name} must be a one-dimensional sequence.")
    if np.any(normalized < 0):
        raise ValueError(f"{parameter_name} must contain non-negative KMeans labels.")
    return normalized if normalized.size else None


def remove_clusters(
    zarr_path: str | Path,
    labels_to_remove: Sequence[int],
    *,
    output_dir: str | Path = "./cluster_segmentation_output",
) -> np.ndarray:
    """Remove every pixel belonging to the specified labels."""
    labels = _normalize_labels(labels_to_remove, parameter_name="labels_to_remove")
    if labels is None:
        raise ValueError("labels_to_remove must contain at least one label.")

    label_image, kmeans_image, umap_image = _read_feature_analysis(zarr_path)
    mask = ((label_image >= 0) & ~np.isin(label_image, labels)).astype(np.uint8)
    _save_outputs(mask, kmeans_image, umap_image, output_dir)
    return mask


def _exterior_background(remaining: np.ndarray) -> np.ndarray:
    """Return non-remaining pixels connected to the image boundary."""
    background = ~remaining
    seeds = np.zeros(remaining.shape, dtype=bool)
    seeds[0, :] = background[0, :]
    seeds[-1, :] = background[-1, :]
    seeds[:, 0] = background[:, 0]
    seeds[:, -1] = background[:, -1]
    return np.asarray(binary_propagation(seeds, mask=background), dtype=bool)


def _exterior_boundary(remaining: np.ndarray, exterior: np.ndarray) -> np.ndarray:
    """Return remaining pixels touching the exterior through 4-neighbour connectivity."""
    padded = np.pad(exterior, 1, mode="constant", constant_values=False)
    boundary = remaining & (
        padded[:-2, 1:-1]
        | padded[2:, 1:-1]
        | padded[1:-1, :-2]
        | padded[1:-1, 2:]
    )
    boundary[0, :] |= remaining[0, :]
    boundary[-1, :] |= remaining[-1, :]
    boundary[:, 0] |= remaining[:, 0]
    boundary[:, -1] |= remaining[:, -1]
    return boundary


def _select_outermost_label(
    label_image: np.ndarray,
    boundary: np.ndarray,
    labels: np.ndarray | None,
) -> int | None:
    """Select the eligible label with the most pixels on the exterior boundary."""
    if labels is not None:
        boundary = boundary & np.isin(label_image, labels)
    boundary_labels, counts = np.unique(label_image[boundary], return_counts=True)
    if boundary_labels.size == 0:
        return None
    return int(boundary_labels[np.argmax(counts)])


def _adjacent_remaining(remaining: np.ndarray, removed: np.ndarray) -> np.ndarray:
    """Return remaining pixels that touch the previously removed layer."""
    padded = np.pad(removed, 1, mode="constant", constant_values=False)
    return remaining & (
        padded[:-2, 1:-1]
        | padded[2:, 1:-1]
        | padded[1:-1, :-2]
        | padded[1:-1, 2:]
    )


def remove_background(
    zarr_path: str | Path,
    labels_to_remove: Sequence[int] | None = None,
    *,
    remove_num: int = 1,
    output_dir: str | Path = "./background_removal_output",
) -> np.ndarray:
    """Remove exterior label-component layers and save transparent result images."""
    if remove_num < 1:
        raise ValueError("remove_num must be a positive integer.")

    labels = _normalize_labels(labels_to_remove, parameter_name="labels_to_remove")
    label_image, kmeans_image, umap_image = _read_feature_analysis(zarr_path)
    remaining = label_image >= 0
    previous_removed: np.ndarray | None = None

    for _ in range(remove_num):
        if previous_removed is None:
            frontier = _exterior_boundary(remaining, _exterior_background(remaining))
        else:
            frontier = _adjacent_remaining(remaining, previous_removed)

        outermost_label = _select_outermost_label(label_image, frontier, labels)
        if outermost_label is None:
            break

        label_region = remaining & (label_image == outermost_label)
        seeds = frontier & label_region
        removed = np.asarray(binary_propagation(seeds, mask=label_region), dtype=bool)
        remaining[removed] = False
        previous_removed = removed

    mask = remaining.astype(np.uint8)
    _save_outputs(mask, kmeans_image, umap_image, output_dir)
    return mask


def keep_clusters(
    zarr_path: str | Path,
    labels_to_keep: Sequence[int],
    *,
    output_dir: str | Path = "./cluster_segmentation_output",
) -> np.ndarray:
    """Keep specified labels and all regions topologically enclosed by them."""
    labels = _normalize_labels(labels_to_keep, parameter_name="labels_to_keep")
    if labels is None:
        raise ValueError("labels_to_keep must contain at least one label.")

    label_image, kmeans_image, umap_image = _read_feature_analysis(zarr_path)
    valid = label_image >= 0
    keep_labels = valid & np.isin(label_image, labels)
    traversable = ~keep_labels

    exterior_seeds = np.zeros(label_image.shape, dtype=bool)
    exterior_seeds[0, :] = traversable[0, :]
    exterior_seeds[-1, :] = traversable[-1, :]
    exterior_seeds[:, 0] = traversable[:, 0]
    exterior_seeds[:, -1] = traversable[:, -1]
    exterior_reachable = np.asarray(
        binary_propagation(exterior_seeds, mask=traversable),
        dtype=bool,
    )

    mask = (valid & ~exterior_reachable).astype(np.uint8)
    _save_outputs(mask, kmeans_image, umap_image, output_dir)
    return mask
