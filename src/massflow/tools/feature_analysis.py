"""UMAP/KMeans feature analysis for MSI data."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import umap
import zarr
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler

from massflow.data_manager import MSDataManagerImzML
from massflow.data_manager import IonImageDataManagerZarr
from massflow.tools.plot_ion_image import read_tic_image
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tools.feature_analysis")


# ---------------------------------------------------------------------------
# Low-level algorithm functions
# ---------------------------------------------------------------------------

def compute_umap_embeddings(
    matrix: np.ndarray,
    *,
    n_components: int = 3,
    metric: str = "cosine",
    random_state: int = 0,
    umap_sample: int | None = None,
    transform_batch_size: int | None = None,
) -> np.ndarray:
    """Compute UMAP embeddings and MinMax-scale to [0, 1]."""
    n_samples = int(matrix.shape[0])
    reducer = umap.UMAP(
        metric=metric,
        n_components=n_components,
        random_state=random_state,
        force_approximation_algorithm=True,
    )

    logger.info(f"Computing UMAP: n_components={n_components}, metric={metric}")

    # method1: fit on all data at once (fast for small datasets, but high memory)
    if umap_sample is None or int(umap_sample) >= n_samples:
        matrix = np.nan_to_num(matrix.astype(np.float32, copy=False), copy=False)
        embeddings = reducer.fit_transform(matrix)
        return MinMaxScaler().fit_transform(embeddings).astype(np.float32)  # type: ignore

    # method2: fit on a random sample, then transform in batches (lower memory, more I/O)
    sample_size = int(umap_sample)

    batch_size = 256 if transform_batch_size is None else int(transform_batch_size)

    rng = np.random.default_rng(random_state)
    sample_indices = np.sort(rng.choice(n_samples, size=sample_size, replace=False))
    sample_matrix = np.array(matrix[sample_indices], dtype=np.float32, copy=True)
    sample_matrix = np.nan_to_num(sample_matrix, copy=False)
    reducer.fit(sample_matrix)

    embeddings = np.empty((n_samples, n_components), dtype=np.float32)
    for start in range(0, n_samples, batch_size):
        end = min(start + batch_size, n_samples)
        batch = np.array(matrix[start:end], dtype=np.float32, copy=True)
        batch = np.nan_to_num(batch, copy=False)
        embeddings[start:end] = reducer.transform(batch).astype(np.float32, copy=False) # type: ignore

    return MinMaxScaler().fit_transform(embeddings).astype(np.float32)  # type: ignore


def kmeans_clustering(
    embeddings: np.ndarray,
    *,
    n_clusters: int = 10,
    random_state: int = 0,
) -> np.ndarray:
    """Perform KMeans clustering on embeddings."""
    logger.info(f"Performing KMeans clustering with {n_clusters} clusters.")

    clusters = KMeans(
        n_clusters=n_clusters,
        random_state=random_state,
        n_init="auto",
    ).fit_predict(embeddings)

    return clusters


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def imzml_build_spectra_matrix(
    manager: Any,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a cleaned spectra matrix from ``matrix_generator`` (two-pass).

    * **Pass 1** – lightweight scan: record nonzero flags and xy
      coordinates per row.  Compute keep mask (nonzero AND last
      occurrence of each xy).
    * **Pass 2** – re-read ``matrix_generator``, pre-allocate a plain
      float32 array, and copy only kept rows into it.
    """
    # ------------------------------------------------------------------
    # Pass 1 – metadata only, no intensity storage
    # ------------------------------------------------------------------
    nonzero_flags: list[np.ndarray] = []
    xy_coords: list[np.ndarray] = []
    n_features: int | None = None
    total_rows = 0

    for _, intensity_matrix, _lengths, coordinates in manager.matrix_generator(include_mz=False):
        intensity = np.asarray(intensity_matrix)
        coordinates = np.asarray(coordinates)
        if n_features is None:
            n_features = int(intensity.shape[1])

        nonzero_flags.append(np.any(intensity[:, :n_features] != 0, axis=1))
        xy_coords.append(coordinates[:, :2].astype(np.int64, copy=False))
        total_rows += int(intensity.shape[0])

    if total_rows == 0:
        raise ValueError("No spectra were read from the source.")
    assert n_features is not None

    # Build keep mask
    all_nonzero = np.concatenate(nonzero_flags)
    all_xy = np.vstack(xy_coords)
    del nonzero_flags, xy_coords

    last_occurrence: dict[tuple[int, int], int] = {}
    for i in range(total_rows):
        last_occurrence[(int(all_xy[i, 0]), int(all_xy[i, 1]))] = i

    is_last = np.zeros(total_rows, dtype=bool)
    for idx in last_occurrence.values():
        is_last[idx] = True
    del last_occurrence

    keep_mask = is_last & all_nonzero
    n_keep = int(keep_mask.sum())
    del all_nonzero, is_last

    if n_keep == 0:
        raise ValueError("No nonzero spectra remain after removing duplicates and all-zero rows.")

    # ------------------------------------------------------------------
    # Pass 2 – re-read and fill pre-allocated output
    # ------------------------------------------------------------------
    spectra_matrix = np.empty((n_keep, n_features), dtype=np.float32)
    out_coordinates = np.empty((n_keep, 2), dtype=np.int64)

    write_cursor = 0
    batch_offset = 0

    for _, intensity_matrix, _lengths, coordinates in manager.matrix_generator(include_mz=False):
        intensity = np.asarray(intensity_matrix)
        coordinates = np.asarray(coordinates)
        batch_size = int(intensity.shape[0])

        batch_keep = keep_mask[batch_offset : batch_offset + batch_size]
        n_batch_keep = int(batch_keep.sum())
        if n_batch_keep > 0:
            spectra_matrix[write_cursor : write_cursor + n_batch_keep] = (
                intensity[batch_keep, :n_features].astype(np.float32, copy=False)
            )
            out_coordinates[write_cursor : write_cursor + n_batch_keep] = (
                coordinates[batch_keep, :2].astype(np.int64, copy=False)
            )
            write_cursor += n_batch_keep
        batch_offset += batch_size

    return spectra_matrix, out_coordinates


def zarr_build_spectra_matrix(
    manager: IonImageDataManagerZarr,
    *,
    ion_block: int = 256,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a cleaned spectra matrix from a Zarr ion-image backend.

    Two-pass: pass 1 scans for nonzero pixels, pass 2 fills a
    pre-allocated ``(n_keep, n_ions)`` float32 array by ion-block.
    """
    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    coordinates = np.asarray(manager.coordinates)[:, :2].astype(np.int64, copy=False)
    flat_idx = manager.flat_idx

    if n_pixels == 0 or n_ions == 0:
        raise ValueError("Zarr source contains no data.")

    # Pass 1 – identify nonzero pixels
    is_nonzero = np.zeros(n_pixels, dtype=bool)
    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        dense_block = np.asarray(manager.dset[ion_start:ion_end])  # (block, H, W)
        pixel_block = dense_block.reshape(ion_end - ion_start, -1)[:, flat_idx]
        is_nonzero |= np.any(pixel_block != 0, axis=0)

    n_keep = int(is_nonzero.sum())
    if n_keep == 0:
        raise ValueError("No nonzero spectra in Zarr source.")

    n_zero = n_pixels - n_keep
    if n_zero:
        logger.info(
            f"Zarr cleaning: {n_pixels} pixels \u2192 {n_keep} kept "
            f"({n_zero} all-zero removed)"
        )

    keep_mask = is_nonzero
    out_coordinates = coordinates[keep_mask]

    # Pass 2 – fill output array column-blocks
    spectra_matrix = np.empty((n_keep, n_ions), dtype=np.float32)
    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        dense_block = np.asarray(manager.dset[ion_start:ion_end])
        pixel_block = dense_block.reshape(ion_end - ion_start, -1)[:, flat_idx]
        spectra_matrix[:, ion_start:ion_end] = pixel_block.T[keep_mask].astype(
            np.float32, copy=False
        )

    return spectra_matrix, out_coordinates


# ---------------------------------------------------------------------------
# Spatial image construction
# ---------------------------------------------------------------------------

def _infer_spatial_params(
    coordinates: np.ndarray,
    *,
    width: int | None = None,
    height: int | None = None,
) -> tuple[np.ndarray, np.ndarray, bool, int, int]:
    """Shared coordinate parsing for image builders.

    Returns ``(scatter_x, scatter_y, zero_based, width, height)``.
    """
    xy = coordinates[:, :2].astype(np.int64, copy=False)
    xs = xy[:, 0]
    ys = xy[:, 1]
    zero_based = bool(
        (xs.size > 0 and int(xs.min()) == 0)
        or (ys.size > 0 and int(ys.min()) == 0)
    )
    if width is None:
        width = int(xs.max()) + 1 if zero_based else int(xs.max())
    if height is None:
        height = int(ys.max()) + 1 if zero_based else int(ys.max())
    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid image dimensions: width={width}, height={height}")

    scatter_x = xs if zero_based else xs - 1
    scatter_y = ys if zero_based else ys - 1
    return scatter_x, scatter_y, zero_based, width, height


def build_umap_image(
    scaled_embedding: np.ndarray,
    coordinates: np.ndarray,
    *,
    width: int | None = None,
    height: int | None = None,
) -> np.ndarray:
    """Map a 3-component scaled embedding to an RGB spatial image."""
    if scaled_embedding.shape[1] != 3:
        raise ValueError(
            f"UMAP image requires 3-component embedding, "
            f"got {scaled_embedding.shape[1]}."
        )

    scatter_x, scatter_y, _, w, h = _infer_spatial_params(
        coordinates, width=width, height=height,
    )

    rgb = np.clip(scaled_embedding, 0.0, 1.0)
    rgb_uint8 = (rgb * 255).astype(np.uint8)

    image = np.zeros((h, w, 3), dtype=np.uint8)
    valid = (scatter_x >= 0) & (scatter_x < w) & (scatter_y >= 0) & (scatter_y < h)
    image[scatter_y[valid], scatter_x[valid]] = rgb_uint8[valid]

    return image


def build_kmeans_image(
    labels: np.ndarray,
    coordinates: np.ndarray,
    *,
    n_clusters: int,
    width: int | None = None,
    height: int | None = None,
    cmap_name: str = "tab20",
) -> np.ndarray:
    """Map cluster labels to a colour-coded spatial image."""
    scatter_x, scatter_y, _, w, h = _infer_spatial_params(
        coordinates, width=width, height=height,
    )

    palette = (matplotlib.colormaps[cmap_name](range(20))[:, :3] * 255).astype(np.uint8)
    if n_clusters > len(palette):
        raise ValueError(f"n_clusters must be <= {len(palette)} when using {cmap_name} colors.")
    colors = palette[labels]

    image = np.zeros((h, w, 3), dtype=np.uint8)
    valid = (scatter_x >= 0) & (scatter_x < w) & (scatter_y >= 0) & (scatter_y < h)
    image[scatter_y[valid], scatter_x[valid]] = colors[valid]

    return image


# ---------------------------------------------------------------------------
# High-level orchestration
# ---------------------------------------------------------------------------


def _validate_parameters(
    n_components: int,
    umap_sample: int | None,
    transform_batch_size: int | None,
    save_jpg: bool,
    save_matrix: bool,
) -> None:
    if n_components != 3:
        raise ValueError("n_components must be 3 for RGB image output.")
    if umap_sample is not None and umap_sample <= 0:
        raise ValueError("umap_sample must be a positive integer or None.")
    if transform_batch_size is not None and transform_batch_size <= 0:
        raise ValueError("transform_batch_size must be a positive integer or None.")
    if not save_jpg and not save_matrix:
        raise ValueError("At least one of save_jpg or save_matrix must be True.")

def _is_continuous(data_manager: Any) -> bool:
    """Heuristic check for continuous/shared m/z axis."""
    if data_manager.ms.meta.continuous is True:
        return True
    return False

def _get_spatial_dimensions(data_manager: Any) -> tuple[int | None, int | None]:
    """Extract (width, height) from data_manager metadata."""
    if isinstance(data_manager, IonImageDataManagerZarr):
        return int(data_manager.width), int(data_manager.height)
    meta = getattr(getattr(data_manager, "ms", None), "meta", None)
    w = getattr(meta, "max_count_of_pixels_x", None)
    h = getattr(meta, "max_count_of_pixels_y", None)
    return (int(w) if w else None), (int(h) if h else None)


def run_analysis(
    data_manager: Any,
    *,
    output_dir: str | None = "./feature_analysis_output",
    save_images: bool = True,
    save_matrix: bool = False,
    n_components: int = 3,
    umap_metric: str = "cosine",
    n_clusters: int = 10,
    random_state: int = 0,
    umap_sample: int | None = None,
    transform_batch_size: int | None = None,
    dpi: int = 300,
) -> bool:
    """Run the full UMAP/KMeans analysis pipeline on processed MSI data."""

    _validate_parameters(n_components, umap_sample, transform_batch_size, save_images, save_matrix)

    out = Path(output_dir) if output_dir is not None else Path("./feature_analysis_output")
    out.mkdir(parents=True, exist_ok=True)

    width, height = _get_spatial_dimensions(data_manager)

    # ---- 1. Build spectra matrix ----------------------------------------
    if isinstance(data_manager, IonImageDataManagerZarr):
        spectra_matrix, coordinates = zarr_build_spectra_matrix(data_manager)
    elif isinstance(data_manager, MSDataManagerImzML):
        if _is_continuous(data_manager) is False:
            logger.error("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")
            raise ValueError("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")
        spectra_matrix, coordinates = imzml_build_spectra_matrix(data_manager)
    else:
        logger.error(f"Unsupported data_manager type: {type(data_manager).__name__}")
        raise TypeError("Unsupported data_manager type.")

    # ---- 2. UMAP --------------------------------------------------------
    scaled_embedding = compute_umap_embeddings(
        spectra_matrix,
        n_components=n_components,
        metric=umap_metric,
        random_state=random_state,
        umap_sample=umap_sample,
        transform_batch_size=transform_batch_size,
    )

    # ---- 3. KMeans on the *same* scaled embedding -----------------------
    labels = kmeans_clustering(
        scaled_embedding,
        n_clusters=n_clusters,
        random_state=random_state,
    )

    # ---- 4. Spatial images ----------------------------------------------
    umap_img = build_umap_image(
        scaled_embedding,
        coordinates,
        width=width,
        height=height,
    )

    kmeans_img = build_kmeans_image(
        labels,
        coordinates,
        n_clusters=n_clusters,
        width=width,
        height=height,
    )

    # ---- 5. TIC ---------------------------------------------------------
    tic_img = read_tic_image(data_manager)

    # ---- 6. Save outputs ------------------------------------------------
    if save_images:
        plt.imsave(str(out / "umap_image.jpg"), umap_img, dpi=dpi)
        plt.imsave(str(out / "kmeans_image.jpg"), kmeans_img, dpi=dpi)

        vmax = np.nanpercentile(tic_img, 95)
        plt.imsave(
            str(out / "tic_image.jpg"),
            tic_img,
            cmap="inferno",
            vmax=vmax,
            dpi=dpi,
        )

    if save_matrix:
        result_group = zarr.open_group(str(out / "result.zarr"), mode="w")
        result_group.create_array("umap_image", data=umap_img)
        result_group.create_array("kmeans_image", data=kmeans_img)
        result_group.create_array("tic_image", data=tic_img)
        result_group.attrs["format"] = "massflow_feature_analysis"

    logger.info(f"UMAP/KMeans analysis complete. Output: {out}")
    return True
