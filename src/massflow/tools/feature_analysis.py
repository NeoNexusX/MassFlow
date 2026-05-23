"""UMAP/KMeans feature analysis for processed MSI data.

Provides both low-level algorithm functions and a high-level
``run_umap_kmeans_analysis`` entry point that orchestrates the full
pipeline: data loading → deduplication → UMAP → KMeans → spatial image
generation → output persistence.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import umap
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler

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
    """Compute UMAP embeddings and MinMax-scale to [0, 1].

    Args:
        matrix: Spectra matrix, shape ``(n_pixels, n_features)``.
        n_components: Target embedding dimensionality.
        metric: Distance metric for UMAP.
        random_state: Random seed for reproducibility.

    Returns:
        Scaled embedding of shape ``(n_pixels, n_components)``, float32.
    """
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
    """Perform KMeans clustering on embeddings.

    Args:
        embeddings: Input array, shape ``(n_samples, n_dims)``.
        n_clusters: Number of clusters.
        random_state: Random seed for reproducibility.

    Returns:
        Cluster labels of shape ``(n_samples,)``.
    """
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

def build_spectra_matrix(
    manager: Any,
    *,
    batch_size: int = 256,
    max_threads: int = 2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a dense spectra matrix from ``matrix_generator``.

    Reads batches of spectra from *manager* and assembles them into a
    single ``(n_pixels, n_features)`` matrix.

    Args:
        manager: An ``MSDataManagerImzML`` instance (must already have
            ``load_head_data`` called, or the caller should ensure it
            is ready to iterate).
        batch_size: Number of spectra per generator batch.
        max_threads: Thread count for the underlying I/O reader.

    Returns:
        A 3-tuple of:

        - **spectra_matrix** – shape ``(n_pixels, n_features)``, float32
        - **coordinates** – shape ``(n_pixels, >=2)``
        - **mz_axis** – shape ``(n_features,)``

    Raises:
        ValueError: If no data is read, or the data is not uniform-length
            (e.g. variable-length processed spectra that have not been
            aligned).
    """
    matrix_parts: list[np.ndarray] = []
    coord_parts: list[np.ndarray] = []
    mz_axis: np.ndarray | None = None
    expected_n_features: int | None = None

    for mz_data, intensity_matrix, lengths, coordinates in manager.matrix_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    ):
        mz_arr = np.asarray(mz_data) if mz_data is not None else None
        intensity = np.asarray(intensity_matrix)
        lengths = np.asarray(lengths, dtype=np.int64)
        coordinates = np.asarray(coordinates)

        # Determine n_features from this batch
        if mz_arr is not None and mz_arr.ndim == 1:
            # Continuous / shared mz: all spectra share the same mz axis
            n_features = int(mz_arr.size)
        else:
            # Processed with per-spectrum mz — require uniform length
            unique_lengths = np.unique(lengths)
            if unique_lengths.size != 1:
                raise ValueError(
                    "Variable-length spectra detected. The data must be "
                    "aligned to a uniform m/z axis before UMAP/KMeans "
                    "analysis. Please run alignment/binning first."
                )
            n_features = int(unique_lengths[0])

        # Consistency check across batches
        if expected_n_features is None:
            expected_n_features = n_features
        elif n_features != expected_n_features:
            raise ValueError(
                f"Inconsistent feature count across batches: "
                f"expected {expected_n_features}, got {n_features}."
            )

        # Capture mz_axis from first batch
        if mz_axis is None:
            if mz_arr is not None and mz_arr.ndim == 1:
                mz_axis = mz_arr.copy()
            elif mz_arr is not None and mz_arr.ndim == 2:
                mz_axis = mz_arr[0, :n_features].copy()
            else:
                mz_axis = np.arange(n_features, dtype=np.float64)

        # Trim intensity to n_features columns (matrix may have padding)
        batch_matrix = intensity[:, :n_features].astype(np.float32, copy=False)
        matrix_parts.append(batch_matrix)
        coord_parts.append(coordinates)

    if not matrix_parts:
        raise ValueError("No spectra were read from the source.")

    spectra_matrix = np.vstack(matrix_parts)
    all_coordinates = np.vstack(coord_parts)

    assert mz_axis is not None  # guaranteed by non-empty matrix_parts

    logger.info(
        f"Built spectra matrix: {spectra_matrix.shape[0]} pixels × "
        f"{spectra_matrix.shape[1]} features"
    )

    return spectra_matrix, all_coordinates, mz_axis


def deduplicate_by_coordinates(
    spectra_matrix: np.ndarray,
    coordinates: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Deduplicate pixels with identical (x, y) coordinates.

    When duplicate coordinates exist, the **last** spectrum is kept
    (consistent with imzML overwrite semantics).

    Args:
        spectra_matrix: shape ``(n_pixels, n_features)``
        coordinates: shape ``(n_pixels, >=2)``

    Returns:
        Deduplicated ``(spectra_matrix, coordinates)`` with the same
        column layout as the inputs.
    """
    n_before = spectra_matrix.shape[0]

    # Build mapping: (x, y) → last index
    coord_to_idx: dict[tuple[int, int], int] = {}
    for i in range(coordinates.shape[0]):
        key = (int(coordinates[i, 0]), int(coordinates[i, 1]))
        coord_to_idx[key] = i

    unique_indices = np.array(sorted(coord_to_idx.values()), dtype=np.intp)

    n_after = unique_indices.shape[0]
    if n_after < n_before:
        logger.info(
            f"Deduplicated coordinates: {n_before} → {n_after} pixels "
            f"({n_before - n_after} duplicates removed)"
        )

    return spectra_matrix[unique_indices], coordinates[unique_indices]


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
    """Map a 3-component scaled embedding to an RGB spatial image.

    Args:
        scaled_embedding: shape ``(n_pixels, 3)``, values in [0, 1].
        coordinates: shape ``(n_pixels, >=2)``, (x, y) columns.
        width: Image width; inferred from coordinates if *None*.
        height: Image height; inferred from coordinates if *None*.

    Returns:
        uint8 RGB image of shape ``(height, width, 3)``.
    """
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
    """Map cluster labels to a colour-coded spatial image.

    Args:
        labels: shape ``(n_pixels,)``
        coordinates: shape ``(n_pixels, >=2)``
        n_clusters: Total number of clusters (for colormap normalisation).
        width: Image width; inferred from coordinates if *None*.
        height: Image height; inferred from coordinates if *None*.
        cmap_name: Matplotlib colormap name.

    Returns:
        uint8 RGB image of shape ``(height, width, 3)``.
    """
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


def _get_spatial_dimensions(manager: Any) -> tuple[int | None, int | None]:
    """Extract (width, height) from manager metadata."""
    meta = getattr(getattr(manager, "ms", None), "meta", None)
    w = getattr(meta, "max_count_of_pixels_x", None)
    h = getattr(meta, "max_count_of_pixels_y", None)
    return (int(w) if w else None), (int(h) if h else None)


def run_umap_kmeans_analysis(
    manager: Any,
    output_dir: str | Path,
    *,
    n_components: int = 3,
    umap_metric: str = "cosine",
    n_clusters: int = 10,
    random_state: int = 0,
    batch_size: int = 256,
    max_threads: int = 2,
    dpi: int = 300,
) -> dict:
    """Run the full UMAP/KMeans analysis pipeline on processed MSI data.

    Args:
        manager: A data manager instance that exposes ``matrix_generator``
            (e.g. ``MSDataManagerImzML``).  Must already have data loaded
            (i.e. ``load_head_data`` called).
        output_dir: Directory for all output files.  Created if absent.
        n_components: UMAP target dimensionality (default 3 for RGB).
        umap_metric: Distance metric passed to UMAP.
        n_clusters: Number of KMeans clusters.
        random_state: Random seed for reproducibility.
        batch_size: Batch size for ``matrix_generator``.
        max_threads: Thread count for the underlying I/O reader.
        dpi: DPI for saved PNG images.

    Returns:
        A metadata ``dict`` summarising results and output file paths.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    width, height = _get_spatial_dimensions(manager)

    # ---- 1. Build spectra matrix --------------------------------------------
    spectra_matrix, coordinates, mz_axis = build_spectra_matrix(
        manager,
        batch_size=batch_size,
        max_threads=max_threads,
    )
    n_pixels_before = spectra_matrix.shape[0]

    # ---- 2. Deduplicate -----------------------------------------------------
    spectra_matrix, coordinates = deduplicate_by_coordinates(
        spectra_matrix, coordinates,
    )

    # ---- 3. UMAP ------------------------------------------------------------
    scaled_embedding = compute_umap_embeddings(
        spectra_matrix,
        n_components=n_components,
        metric=umap_metric,
        random_state=random_state,
    )

    # ---- 4. KMeans on the *same* scaled embedding ---------------------------
    labels = kmeans_clustering(
        scaled_embedding,
        n_clusters=n_clusters,
        random_state=random_state,
    )

    # ---- 5. Spatial images --------------------------------------------------
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

    # ---- 6. Save outputs ----------------------------------------------------
    plt.imsave(str(out / "umap_image.png"), umap_img, dpi=dpi)
    plt.imsave(str(out / "kmeans_image.png"), kmeans_img, dpi=dpi)

    np.save(str(out / "umap_embedding.npy"), scaled_embedding)
    np.save(str(out / "kmeans_labels.npy"), labels)
    np.save(str(out / "coordinates.npy"), coordinates)
    np.save(str(out / "mz_axis.npy"), mz_axis)

    img_h, img_w = umap_img.shape[:2]
    metadata = {
        "image_height": img_h,
        "image_width": img_w,
        "n_pixels": int(spectra_matrix.shape[0]),
        "n_pixels_before_dedup": int(n_pixels_before),
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
            "mz_axis": "mz_axis.npy",
        },
    }

    with open(out / "metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    logger.info(f"UMAP/KMeans analysis complete. Output: {out}")
    return metadata
