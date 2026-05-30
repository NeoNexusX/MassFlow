"""UMAP/KMeans feature analysis for MSI data."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterator

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


def compute_streaming_umap_embeddings(
    data_manager: Any,
    *,
    n_components: int = 3,
    metric: str = "cosine",
    random_state: int | None = None,
    sample_ratio: float = 0.2,
    transform_batch_size: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample-fit + batch-transform UMAP without building the full matrix."""
    reducer = umap.UMAP(
        metric=metric,
        n_components=n_components,
        random_state=random_state,
    )

    # ---- fit on sample ----
    if isinstance(data_manager, MSDataManagerImzML):
        sample_matrix = _read_imzml_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        )
    elif isinstance(data_manager, IonImageDataManagerZarr):
        sample_matrix = _read_zarr_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        )
    else:
        raise TypeError(f"Unsupported data_manager type: {type(data_manager).__name__}")

    sample_matrix = np.nan_to_num(sample_matrix, copy=False)
    logger.info(
        f"Streaming UMAP fit: sample_ratio={sample_ratio}, "
        f"sample_matrix shape={sample_matrix.shape}"
    )
    reducer.fit(sample_matrix)
    del sample_matrix

    # ---- transform in batches ----
    embeddings_parts: list[np.ndarray] = []
    coordinates_parts: list[np.ndarray] = []

    if isinstance(data_manager, MSDataManagerImzML):
        for batch_matrix, batch_coords in _iter_imzml_matrix_batches(
            data_manager, batch_size=transform_batch_size,
        ):
            batch_matrix = np.nan_to_num(batch_matrix, copy=False)
            embeddings_parts.append(
                reducer.transform(batch_matrix).astype(np.float32, copy=False)  # type: ignore
            )
            coordinates_parts.append(batch_coords)
    else:  # IonImageDataManagerZarr
        for batch_matrix, batch_coords in _iter_zarr_matrix_batches(
            data_manager,
        ):
            batch_matrix = np.nan_to_num(batch_matrix, copy=False)
            embeddings_parts.append(
                reducer.transform(batch_matrix).astype(np.float32, copy=False)  # type: ignore
            )
            coordinates_parts.append(batch_coords)

    if not embeddings_parts:
        raise ValueError("No non-zero spectra found during transform.")

    embeddings = np.vstack(embeddings_parts)
    coordinates = np.vstack(coordinates_parts)

    # imzML may have duplicate xy across batches; Zarr does not.
    if isinstance(data_manager, MSDataManagerImzML):
        embeddings, coordinates = _deduplicate_last_by_xy(embeddings, coordinates)

    scaled = MinMaxScaler().fit_transform(embeddings).astype(np.float32)  # type: ignore
    return scaled, coordinates


def kmeans_clustering(
    embeddings: np.ndarray,
    *,
    n_clusters: int = 10,
    random_state: int | None = None,
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
#
# Backend helpers are intentionally grouped by data-manager type.  Each backend
# should provide:
#   1. a sample reader for UMAP fit
#   2. a batch iterator for UMAP transform


# --- imzML backend ---------------------------------------------------------
#
# imzML reads spectra from the raw spectrum list for sampling, then uses the
# data manager's matrix_generator for sequential batch transforms.
def _read_imzml_sample_matrix(
    manager: MSDataManagerImzML,
    sample_ratio: float,
    *,
    random_state: int | None = None,
) -> np.ndarray:
    """Randomly sample spectra from *manager.ms* by ratio and return a matrix."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    n_spectra = len(manager.ms)
    n_features = len(manager.ms.shared_mz_list)  # type: ignore
    sample_count = max(1, int(np.ceil(n_spectra * sample_ratio)))

    rng = np.random.default_rng(random_state)
    sample_indices = np.sort(rng.choice(n_spectra, size=sample_count, replace=False))

    rows: list[np.ndarray] = []
    for idx in sample_indices:
        spectrum = manager.ms[int(idx)]
        intensity = np.asarray(spectrum.intensity, dtype=np.float32)
        spectrum.clear_data()
        if np.any(intensity != 0):
            rows.append(intensity[:n_features])

    if not rows:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")

    return np.vstack(rows)


def _iter_imzml_matrix_batches(
    manager: MSDataManagerImzML,
    *,
    batch_size: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield ``(batch_matrix, batch_coordinates)`` from *manager.matrix_generator*."""
    for _, intensity_matrix, _lengths, coordinates in manager.matrix_generator(
        batch_size=batch_size, include_mz=False,
    ):
        intensity = np.asarray(intensity_matrix, dtype=np.float32)
        coords = np.asarray(coordinates)[:, :2]

        nonzero_mask = np.any(intensity != 0, axis=1)
        if not np.any(nonzero_mask):
            continue

        yield intensity[nonzero_mask], coords[nonzero_mask]


def _deduplicate_last_by_xy(
    values: np.ndarray,
    coordinates: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Keep only the *last* occurrence of each unique (x, y) coordinate."""
    xy = coordinates[:, :2].astype(np.int64, copy=False)
    n = xy.shape[0]

    last: dict[tuple[int, int], int] = {}
    for i in range(n):
        last[(int(xy[i, 0]), int(xy[i, 1]))] = i

    keep = np.array(sorted(last.values()), dtype=np.intp)
    return values[keep], coordinates[keep]


# --- Zarr ion-image backend ------------------------------------------------
#
# Zarr stores ion images as (ion, y, x).  These helpers read selected pixel
# columns from ion blocks and present them as spectra rows.
def _read_zarr_pixel_columns(
    manager: IonImageDataManagerZarr,
    pixel_flat_idx: np.ndarray,
    out: np.ndarray,
    ion_block: int,
) -> None:
    """Fill *out* ``(n_pixels_subset, n_ions)`` by reading ion-blocks from *manager.dset*."""
    n_ions = int(manager.n_ions)
    h, w = int(manager.height), int(manager.width)
    hw = h * w

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        dense_block = np.asarray(manager.dset[ion_start:ion_end])  # (block, H, W)
        flat = dense_block.reshape(ion_end - ion_start, hw)        # (block, H*W)
        out[:, ion_start:ion_end] = flat[:, pixel_flat_idx].T


def _read_zarr_sample_matrix(
    manager: IonImageDataManagerZarr,
    sample_ratio: float,
    *,
    random_state: int | None = None,
    ion_block: int = 256,
) -> np.ndarray:
    """Randomly sample pixels from a Zarr backend by ratio."""
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


def _iter_zarr_matrix_batches(
    manager: IonImageDataManagerZarr,
    *,
    pixel_batch_size: int = 20000,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield ``(batch_matrix, batch_coordinates)`` by pixel-batch from a Zarr backend."""
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

        yield batch_matrix[nonzero_mask], batch_coords[nonzero_mask]


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
    sample_ratio: float,
    transform_batch_size: int | None,
    save_images: bool,
    save_matrix: bool,
) -> None:
    if n_components != 3:
        raise ValueError("n_components must be 3 for RGB image output.")
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")
    if transform_batch_size is not None and transform_batch_size <= 0:
        raise ValueError("transform_batch_size must be a positive integer or None.")
    if not save_images and not save_matrix:
        raise ValueError("At least one of save_images or save_matrix must be True.")

def _is_continuous(data_manager: Any) -> bool:
    """Heuristic check for continuous/shared m/z axis."""
    meta = getattr(getattr(data_manager, "ms", None), "meta", None)
    continuous = getattr(meta, "continuous", None)
    if continuous is True:
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
    random_state: int | None = None,
    sample_ratio: float = 0.3,
    transform_batch_size: int | None = 2048,
    dpi: int = 300,
) -> bool:
    """Run the full UMAP/KMeans analysis pipeline on processed MSI data."""

    # --- 1. Validate parameters and prepare output directory ------------
    _validate_parameters(n_components, sample_ratio, transform_batch_size, save_images, save_matrix)

    out = Path(output_dir) if output_dir is not None else Path("./feature_analysis_output")
    out.mkdir(parents=True, exist_ok=True)

    if isinstance(data_manager, MSDataManagerImzML) and _is_continuous(data_manager) is False:
        logger.error("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")
        raise ValueError("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")

    width, height = _get_spatial_dimensions(data_manager)
    batch_size = 512 if transform_batch_size is None else int(transform_batch_size)

    # --- 2. UMAP embedding with streaming fit/transform ----------------
    scaled_embedding, coordinates = compute_streaming_umap_embeddings(
        data_manager,
        n_components=n_components,
        metric=umap_metric,
        random_state=random_state,
        sample_ratio=sample_ratio,
        transform_batch_size=batch_size,
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
        plt.imsave(str(out / "umap_image.png"), umap_img, dpi=dpi)
        plt.imsave(str(out / "kmeans_image.png"), kmeans_img, dpi=dpi)

        vmax = np.nanpercentile(tic_img, 95)
        plt.imsave(
            str(out / "tic_image.png"),
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
