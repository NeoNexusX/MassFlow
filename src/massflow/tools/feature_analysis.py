"""UMAP/KMeans feature analysis for MSI data."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import umap
import zarr
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler

from massflow.data_manager import MSDataManagerImzML
from massflow.data_manager import IonImageDataManagerZarr
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
# Data loading helpers 还需要完成另外三个data_manager的适配 TODO
# ---------------------------------------------------------------------------

def imzml_build_spectra_matrix(
    manager: Any,
    *,
    matrix_path: str | Path,
) -> tuple[np.memmap, np.ndarray]:
    """Build a cleaned spectra matrix from ``matrix_generator``."""

    # ------------------------------------------------------------------
    # Pass 1 – lightweight metadata scan (no intensity storage)
    # ------------------------------------------------------------------
    nonzero_flags: list[np.ndarray] = []   # bool per row
    xy_coords: list[np.ndarray] = []       # (batch_n, 2) int64
    n_features: int | None = None
    total_rows = 0

    for _, intensity_matrix, _lengths, coordinates in manager.matrix_generator(include_mz=False):
        intensity = np.asarray(intensity_matrix)
        coordinates = np.asarray(coordinates)

        if n_features is None:
            n_features = int(intensity.shape[1])
        assert n_features is not None

        # Record nonzero flags (no intensity data saved)
        batch_matrix = intensity[:, :n_features]
        nonzero_flags.append(np.any(batch_matrix != 0, axis=1))

        # Record xy only
        xy_coords.append(coordinates[:, :2].astype(np.int64, copy=False))

        total_rows += int(intensity.shape[0])

    if total_rows == 0:
        raise ValueError("No spectra were read from the source.")
    assert n_features is not None

    # Concatenate lightweight arrays
    all_nonzero = np.concatenate(nonzero_flags)        # (total_rows,)
    all_xy = np.vstack(xy_coords)                      # (total_rows, 2)

    # Determine which rows to keep:
    #   1. For each (x, y), find the *last* global index.
    #   2. A row is kept iff it is that last occurrence AND it is nonzero.
    last_occurrence: dict[tuple[int, int], int] = {}
    for i in range(total_rows):
        last_occurrence[(int(all_xy[i, 0]), int(all_xy[i, 1]))] = i

    is_last = np.zeros(total_rows, dtype=bool)
    for idx in last_occurrence.values():
        is_last[idx] = True

    keep_mask = is_last & all_nonzero   # bool (total_rows,)
    n_keep = int(keep_mask.sum())

    if n_keep == 0:
        logger.error("No nonzero spectra remain after removing duplicates and all-zero rows.")
        raise ValueError("No nonzero spectra remain after removing duplicates and all-zero rows.")

    # Free pass-1 intermediates
    del nonzero_flags, xy_coords, all_nonzero, is_last, last_occurrence

    # ------------------------------------------------------------------
    # Pass 2 – stream again, write only kept rows into memmap
    # ------------------------------------------------------------------
    spectra_matrix = np.memmap(
        str(matrix_path),
        dtype=np.float32,
        mode="w+",
        shape=(n_keep, n_features),
    )
    out_coordinates = np.empty((n_keep, 2), dtype=np.int64)

    write_cursor = 0
    batch_offset = 0

    for _, intensity_matrix, _lengths, coordinates in manager.matrix_generator(include_mz=False):
        intensity = np.asarray(intensity_matrix)
        coordinates = np.asarray(coordinates)
        batch_size_actual = int(intensity.shape[0])

        batch_keep = keep_mask[batch_offset : batch_offset + batch_size_actual]
        n_batch_keep = int(batch_keep.sum())

        if n_batch_keep > 0:
            kept_rows = intensity[batch_keep, :n_features].astype(np.float32, copy=False)
            spectra_matrix[write_cursor : write_cursor + n_batch_keep] = kept_rows
            out_coordinates[write_cursor : write_cursor + n_batch_keep] = (
                coordinates[batch_keep, :2].astype(np.int64, copy=False)
            )
            write_cursor += n_batch_keep

        batch_offset += batch_size_actual

    spectra_matrix.flush()

    return spectra_matrix, out_coordinates


def zarr_build_spectra_matrix(
    manager: IonImageDataManagerZarr,
    *,
    matrix_path: str | Path,
    ion_block: int = 256,
) -> tuple[np.memmap, np.ndarray]:
    """Build a cleaned spectra matrix from a Zarr ion-image backend."""
    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    coordinates = np.asarray(manager.coordinates)[:, :2].astype(np.int64, copy=False)
    flat_idx = manager.flat_idx

    if n_pixels == 0 or n_ions == 0:
        raise ValueError("Zarr source contains no data.")

    # ------------------------------------------------------------------
    # Pass 1 – identify nonzero pixels (no full matrix stored)
    # ------------------------------------------------------------------
    # Accumulate a per-pixel "has any nonzero intensity" flag by streaming
    # ion blocks, avoiding loading the full (n_pixels, n_ions) matrix.
    is_nonzero = np.zeros(n_pixels, dtype=bool)

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        # dset shape: (n_ions, H, W) — read a contiguous ion block
        dense_block = np.asarray(manager.dset[ion_start:ion_end])  # (block, H, W)
        block_size = ion_end - ion_start
        # Extract pixel vectors: (block, n_pixels) via flat_idx
        pixel_block = dense_block.reshape(block_size, -1)[:, flat_idx]
        # Any nonzero across these ions marks the pixel as nonzero
        is_nonzero |= np.any(pixel_block != 0, axis=0)

    n_keep = int(is_nonzero.sum())
    if n_keep == 0:
        raise ValueError("No nonzero spectra in Zarr source.")

    n_zero = n_pixels - n_keep
    if n_zero:
        logger.info(
            f"Zarr cleaning: {n_pixels} pixels → {n_keep} kept "
            f"({n_zero} all-zero removed)"
        )

    keep_mask = is_nonzero  # bool (n_pixels,)
    out_coordinates = coordinates[keep_mask]

    # ------------------------------------------------------------------
    # Pass 2 – fill memmap column-blocks
    # ------------------------------------------------------------------
    spectra_matrix = np.memmap(
        str(matrix_path),
        dtype=np.float32,
        mode="w+",
        shape=(n_keep, n_ions),
    )

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        dense_block = np.asarray(manager.dset[ion_start:ion_end])
        block_size = ion_end - ion_start
        pixel_block = dense_block.reshape(block_size, -1)[:, flat_idx]  # (block, n_pixels)
        # Transpose to (n_pixels, block), filter, write columns
        spectra_matrix[:, ion_start:ion_end] = pixel_block.T[keep_mask].astype(
            np.float32, copy=False
        )

    spectra_matrix.flush()

    logger.info(
        f"Built spectra matrix from Zarr: {n_keep} pixels × {n_ions} features"
        f" (memmap: {matrix_path})"
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

    cmap = matplotlib.colormaps[cmap_name]
    norm = Normalize(vmin=0, vmax=max(n_clusters - 1, 1))
    colors = (cmap(norm(labels))[:, :3] * 255).astype(np.uint8)

    image = np.zeros((h, w, 3), dtype=np.uint8)
    valid = (scatter_x >= 0) & (scatter_x < w) & (scatter_y >= 0) & (scatter_y < h)
    image[scatter_y[valid], scatter_x[valid]] = colors[valid]

    return image


# ---------------------------------------------------------------------------
# High-level orchestration
# ---------------------------------------------------------------------------


def _validate_parameters(
    umap_sample: int | None,
    transform_batch_size: int | None,
) -> None:
    if umap_sample is not None and umap_sample <= 0:
        raise ValueError("umap_sample must be a positive integer or None.")
    if transform_batch_size is not None and transform_batch_size <= 0:
        raise ValueError("transform_batch_size must be a positive integer or None.")

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


def run_umap_kmeans_analysis(
    data_manager: Any,
    output_dir: str | None = "./temp",
    *,
    n_components: int = 3,
    umap_metric: str = "cosine",
    n_clusters: int = 10,
    random_state: int = 0,
    umap_sample: int | None = None,
    transform_batch_size: int | None = None,
    dpi: int = 300,
) -> dict:
    """Run the full UMAP/KMeans analysis pipeline on processed MSI data."""

    _validate_parameters(umap_sample, transform_batch_size)

    out = Path(output_dir) if output_dir is not None else Path("./temp")
    out.mkdir(parents=True, exist_ok=True)

    width, height = _get_spatial_dimensions(data_manager)

    # Temporary memmap file for the spectra matrix – cleaned up in finally.
    fd, tmp_matrix_name = tempfile.mkstemp(
        prefix="_tmp_spectra_matrix_",
        suffix=".dat",
        dir=str(out),
    )
    os.close(fd)
    tmp_matrix_path = Path(tmp_matrix_name)
    spectra_matrix: np.memmap | None = None

    try:
        # ---- 1. Build spectra matrix ----------------------------------------
        if isinstance(data_manager, IonImageDataManagerZarr):
            spectra_matrix, coordinates = zarr_build_spectra_matrix(
                data_manager,
                matrix_path=tmp_matrix_path,
            )
        elif isinstance(data_manager, MSDataManagerImzML):
            if _is_continuous(data_manager) is False:
                logger.error("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")
                raise ValueError("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")
            spectra_matrix, coordinates = imzml_build_spectra_matrix(
                data_manager,
                matrix_path=tmp_matrix_path,
            )
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

        # ---- 5. Save outputs ------------------------------------------------
        plt.imsave(str(out / "umap_image.png"), umap_img, dpi=dpi)
        plt.imsave(str(out / "kmeans_image.png"), kmeans_img, dpi=dpi)

        img_h, img_w = umap_img.shape[:2]
        metadata = {
            "image_height": img_h,
            "image_width": img_w,
            "n_pixels": int(spectra_matrix.shape[0]),
            "n_features": int(spectra_matrix.shape[1]),
            "n_clusters": int(n_clusters),
            "n_components": int(n_components),
            "umap_metric": umap_metric,
            "random_state": int(random_state),
        }

        result_group = zarr.open_group(str(out / "result.zarr"), mode="w")
        result_group.attrs.update(metadata)
        result_group.create_array("umap_embedding", data=scaled_embedding)
        result_group.create_array("kmeans_labels", data=np.asarray(labels, dtype=np.int32))
        result_group.create_array("coordinates", data=coordinates)
        result_group.create_array("umap_image", data=umap_img)
        result_group.create_array("kmeans_image", data=kmeans_img)

        logger.info(f"UMAP/KMeans analysis complete. Output: {out}")
        return metadata

    finally:
        # Release memmap reference and delete the temp file.
        del spectra_matrix
        if tmp_matrix_path.exists():
            try:
                os.remove(tmp_matrix_path)
            except OSError:
                logger.warning(f"Failed to remove temp memmap file: {tmp_matrix_path}")
