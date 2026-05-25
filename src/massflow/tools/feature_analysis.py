"""UMAP/KMeans feature analysis for MSI data."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import umap
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler

from massflow.data_manager import MSDataManagerImzML
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
) -> np.ndarray:
    """Compute UMAP embeddings and MinMax-scale to [0, 1]."""
    logger.info(
        f"Computing UMAP: n_pixels={matrix.shape[0]}, "
        f"n_features={matrix.shape[1]}, n_components={n_components}, "
        f"metric={metric}"
    )
    matrix = np.nan_to_num(matrix.astype(np.float32, copy=False), copy=False)

    embeddings = umap.UMAP(
        metric=metric,
        n_components=n_components,
        random_state=random_state,
    ).fit_transform(matrix)

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
    batch_size: int = 256,
    max_threads: int = 2,
    matrix_path: str | Path,
) -> tuple[np.memmap, np.ndarray, np.ndarray]:
    """Build a cleaned spectra matrix from ``matrix_generator``."""
    gen_kwargs = dict(batch_size=batch_size, include_mz=True, max_threads=max_threads)

    # ------------------------------------------------------------------
    # Pass 1 – lightweight metadata scan (no intensity storage)
    # ------------------------------------------------------------------
    nonzero_flags: list[np.ndarray] = []   # bool per row
    xy_coords: list[np.ndarray] = []       # (batch_n, 2) int64
    mz_axis: np.ndarray | None = None
    n_features: int | None = None
    total_rows = 0

    for mz_data, intensity_matrix, lengths, coordinates in manager.matrix_generator(**gen_kwargs):
        mz_arr = np.asarray(mz_data) if mz_data is not None else None
        intensity = np.asarray(intensity_matrix)
        lengths_arr = np.asarray(lengths, dtype=np.int64)
        coordinates = np.asarray(coordinates)

        # Determine / validate n_features
        if mz_arr is not None and mz_arr.ndim == 1:
            batch_n_features = int(mz_arr.size)
        else:
            unique_lengths = np.unique(lengths_arr)
            if unique_lengths.size != 1:
                raise ValueError(
                    "Variable-length spectra detected. The data must be "
                    "aligned to a uniform m/z axis before UMAP/KMeans "
                    "analysis. Please run alignment/binning first."
                )
            batch_n_features = int(unique_lengths[0])

        if n_features is None:
            n_features = batch_n_features
        elif batch_n_features != n_features:
            raise ValueError(
                f"Inconsistent feature count across batches: "
                f"expected {n_features}, got {batch_n_features}."
            )

        # Capture mz_axis from the first batch
        if mz_axis is None:
            if mz_arr is not None and mz_arr.ndim == 1:
                mz_axis = mz_arr.copy()
            elif mz_arr is not None and mz_arr.ndim == 2:
                mz_axis = mz_arr[0, :n_features].copy()
            else:
                mz_axis = np.arange(n_features, dtype=np.float64)

        # Record nonzero flags (no intensity data saved)
        batch_matrix = intensity[:, :n_features]
        nonzero_flags.append(np.any(batch_matrix != 0, axis=1))

        # Record xy only
        xy_coords.append(coordinates[:, :2].astype(np.int64, copy=False))

        total_rows += int(intensity.shape[0])

    if total_rows == 0:
        raise ValueError("No spectra were read from the source.")
    assert n_features is not None and mz_axis is not None

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

    for mz_data, intensity_matrix, lengths, coordinates in manager.matrix_generator(**gen_kwargs):
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

    return spectra_matrix, out_coordinates, mz_axis


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


def _is_continuous(data_manager: Any) -> bool:
    """Heuristic check for continuous/shared m/z axis."""
    if isinstance(data_manager, MSDataManagerImzML):
        if data_manager.ms.meta.continuous is True:
            return True
    return False

def _get_spatial_dimensions(data_manager: Any) -> tuple[int | None, int | None]:
    """Extract (width, height) from data_manager metadata."""
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
    batch_size: int = 256,
    max_threads: int = 2,
    dpi: int = 300,
) -> dict:
    """Run the full UMAP/KMeans analysis pipeline on processed MSI data."""

    out = Path(output_dir) if output_dir is not None else Path("./temp")
    out.mkdir(parents=True, exist_ok=True)

    if _is_continuous(data_manager) is False:
        logger.error("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")
        raise ValueError("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")

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
        # ---- 1. Build spectra matrix (includes dedup + zero-row removal) ----
        spectra_matrix, coordinates, _ = imzml_build_spectra_matrix(
            data_manager,
            batch_size=batch_size,
            max_threads=max_threads,
            matrix_path=tmp_matrix_path,
        )

        # ---- 2. UMAP --------------------------------------------------------
        scaled_embedding = compute_umap_embeddings(
            spectra_matrix,
            n_components=n_components,
            metric=umap_metric,
            random_state=random_state,
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

        np.save(str(out / "umap_embedding.npy"), scaled_embedding)
        np.save(str(out / "kmeans_labels.npy"), labels)
        np.save(str(out / "coordinates.npy"), coordinates)

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
            "files": {
                "umap_image": "umap_image.png",
                "kmeans_image": "kmeans_image.png",
                "umap_embedding": "umap_embedding.npy",
                "kmeans_labels": "kmeans_labels.npy",
                "coordinates": "coordinates.npy",
            },
        }

        with open(out / "metadata.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

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
