"""UMAP/KMeans feature analysis for MSI data."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterator, Literal

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import umap
import zarr
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
    MSDataManagerZarr,
    ZARR_MODE_CONTINUOUS,
)
from massflow.segmentation.segmentation_pipeline import TAB20_COLOR_NAMES, UmapKmeansResult
from massflow.tools.funs import zero_based_xy
from massflow.tools.logger import get_logger

logger = get_logger("massflow.segmentation.feature_analysis")


def _nan_to_num_inplace_chunks(matrix: np.ndarray, row_chunk_size: int = 128) -> np.ndarray:
    """Replace NaN/Inf values in-place without one large temporary mask."""
    if matrix.ndim <= 1:
        np.nan_to_num(matrix, copy=False)
        return matrix

    for row_start in range(0, matrix.shape[0], row_chunk_size):
        np.nan_to_num(matrix[row_start:row_start + row_chunk_size], copy=False)
    return matrix


# ---------------------------------------------------------------------------
# Low-level algorithm functions
# ---------------------------------------------------------------------------

# --- UMAP ------------------------------------------------------------------
def enable_umap_parallel_transform(reducer: umap.UMAP) -> None:
    """Monkey-patch UMAP reducer to enable parallel batch transforms."""
    index = getattr(reducer, "_knn_search_index", None) # pylint: disable=protected-access
    if index is None:
        raise ValueError("UMAP transform requires at least two fitted samples.")
    index.parallel_batch_queries = True # type: ignore

    for attr in ("_search_function", "_deheap_function", "_rerank_function"):
        if hasattr(index, attr):
            delattr(index, attr)

    index.prepare() # type: ignore


def _create_umap_reducer(
    *,
    n_components: int,
    metric: str,
    random_state: int | None,
) -> umap.UMAP:
    """Create the UMAP reducer used by feature analysis."""
    return umap.UMAP(
        metric=metric,
        output_metric="euclidean",
        n_components=n_components,
        random_state=random_state,
        n_jobs=-1,
        force_approximation_algorithm=True,
    )


def compute_streaming_umap_embeddings(
    data_manager: Any,
    *,
    n_components: int = 3,
    metric: str = "cosine",
    random_state: int | None = None,
    sample_ratio: float = 0.2,
    transform_batch_size: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute UMAP embeddings with sample-fit and batch-transform."""
    reducer = _create_umap_reducer(
        n_components=n_components,
        metric=metric,
        random_state=random_state,
    )
    sample_matrix, full_coordinates = _read_umap_sample_matrix(
        data_manager,
        sample_ratio,
        random_state=random_state,
    )

    sample_matrix = _nan_to_num_inplace_chunks(sample_matrix)
    logger.info(
        f"Streaming UMAP fit: sample_ratio={sample_ratio}, "
        f"sample_matrix shape={sample_matrix.shape}"
    )
    if sample_matrix.shape[0] < 2:
        raise ValueError("UMAP transform requires at least two non-zero fitted samples.")
    reducer.fit(sample_matrix)
    if sample_ratio == 1.0:
        embeddings = reducer.transform(sample_matrix).astype(np.float32, copy=False) # type: ignore
        coordinates = full_coordinates
        del sample_matrix
        return _finalize_umap_embeddings(data_manager, embeddings, coordinates) # type: ignore
    del sample_matrix

    if random_state is None:
        enable_umap_parallel_transform(reducer)

    embeddings, coordinates = _transform_umap_batches(
        data_manager,
        reducer,
        batch_size=transform_batch_size,
    )
    return _finalize_umap_embeddings(data_manager, embeddings, coordinates)


def _read_umap_sample_matrix(
    data_manager: Any,
    sample_ratio: float,
    *,
    random_state: int | None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Read the matrix used to fit UMAP and optional full coordinates."""
    if sample_ratio == 1.0:
        sample_matrix, coordinates = _read_full_matrix(data_manager, batch_size=4096)
        return sample_matrix, coordinates
    if isinstance(data_manager, MSDataManagerImzML):
        return _read_imzml_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        ), None
    if isinstance(data_manager, MSDataManagerZarr):
        return _read_spectra_zarr_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        ), None
    if isinstance(data_manager, _ION_IMAGE_MANAGERS):
        return _read_ion_image_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        ), None
    raise TypeError(f"Unsupported data_manager type: {type(data_manager).__name__}")


def _iter_umap_matrix_batches(
    data_manager: Any,
    *,
    batch_size: int,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Iterate nonzero spectra batches for UMAP transform."""
    if isinstance(data_manager, (MSDataManagerImzML, MSDataManagerZarr)):
        return _iter_spectra_matrix_batches(data_manager, batch_size=batch_size)
    if isinstance(data_manager, _ION_IMAGE_MANAGERS):
        return _iter_ion_image_matrix_batches(data_manager, batch_size=batch_size)
    raise TypeError(f"Unsupported data_manager type: {type(data_manager).__name__}")


def _estimate_umap_transform_batches(data_manager: Any, batch_size: int) -> int:
    """Estimate source batch count for UMAP transform progress logging."""
    if isinstance(data_manager, MSDataManagerImzML):
        total_items = len(data_manager.ms)
    elif isinstance(data_manager, MSDataManagerZarr):
        total_items = len(data_manager.ms)
        if total_items == 0:
            reader = data_manager._require_reader()  # pylint: disable=protected-access
            total_items = reader.n_pixels
    elif isinstance(data_manager, _ION_IMAGE_MANAGERS):
        total_items = int(data_manager.n_pixels)
    else:
        raise TypeError(f"Unsupported data_manager type: {type(data_manager).__name__}")
    return max(1, (total_items + batch_size - 1) // batch_size)


def _transform_umap_batches(
    data_manager: Any,
    reducer: umap.UMAP,
    *,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform nonzero spectra batches with a fitted UMAP reducer."""
    embeddings_parts: list[np.ndarray] = []
    coordinates_parts: list[np.ndarray] = []
    total_batches = _estimate_umap_transform_batches(data_manager, batch_size)
    transformed_batches = 0

    for batch_idx, (batch_matrix, batch_coords) in enumerate(
        _iter_umap_matrix_batches(data_manager, batch_size=batch_size),
        start=1,
    ):
        transformed_batches = batch_idx
        batch_matrix = _nan_to_num_inplace_chunks(batch_matrix)
        embeddings_parts.append(reducer.transform(batch_matrix).astype(np.float32, copy=False)) # type: ignore
        coordinates_parts.append(batch_coords)
        if total_batches <= 10 or batch_idx % max(1, total_batches // 10) == 0 or batch_idx == total_batches:
            progress = (batch_idx / total_batches) * 100
            logger.info(
                f"Transforming spectra with UMAP: transformed {batch_idx}/{total_batches} "
                f"batches ({progress:.1f}%)."
            )

    if not embeddings_parts:
        raise ValueError("No non-zero spectra found during transform.")

    logger.info(f"UMAP transform complete: transformed {transformed_batches} nonzero batches.")
    return np.vstack(embeddings_parts), np.vstack(coordinates_parts)


def _finalize_umap_embeddings(
    data_manager: Any,
    embeddings: np.ndarray,
    coordinates: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Deduplicate spectra coordinates and scale UMAP embeddings to RGB range."""
    if coordinates is None:
        raise ValueError("UMAP coordinates are missing.")
    if isinstance(data_manager, (MSDataManagerImzML, MSDataManagerZarr)):
        embeddings, coordinates = _deduplicate_last_by_xy(embeddings, coordinates)

    scaled = MinMaxScaler().fit_transform(embeddings).astype(np.float32)  # type: ignore
    return scaled, coordinates


# --- KMeans ----------------------------------------------------------------
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
# Each backend family provides:
#   1. a sample reader for UMAP fit
#   2. a batch iterator for UMAP transform
# Spectrum managers (imzML / pixel-major Zarr) stream rows through
# ``matrix_generator``; ion-image managers (Binary / H5 / Zarr) read pixel
# columns out of their (ion, pixel) blocks.

_ION_IMAGE_MANAGERS = (IonImageDataManagerZarr, IonImageDataManagerH5, IonImageDataManagerBinary)


def _nonzero_rows(
    matrix: np.ndarray,
    coordinates: np.ndarray,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Drop all-zero spectra; return None when nothing is left."""
    nonzero_mask = np.any(matrix != 0, axis=1)
    if not np.any(nonzero_mask):
        return None
    if np.all(nonzero_mask):
        return matrix, coordinates
    return matrix[nonzero_mask], coordinates[nonzero_mask]


def _read_full_matrix(
    data_manager: Any,
    *,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a full nonzero spectra matrix and coordinate matrix from the active backend."""
    matrix_parts: list[np.ndarray] = []
    coordinates_parts: list[np.ndarray] = []
    for batch_matrix, batch_coords in _iter_umap_matrix_batches(data_manager, batch_size=batch_size):
        matrix_parts.append(batch_matrix)
        coordinates_parts.append(batch_coords)

    if not matrix_parts:
        raise ValueError("No non-zero spectra found.")

    return np.vstack(matrix_parts), np.vstack(coordinates_parts)


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


# --- Spectrum backends (imzML / pixel-major Zarr) ---------------------------


def _iter_spectra_matrix_batches(
    manager: MSDataManagerImzML | MSDataManagerZarr,
    *,
    batch_size: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield nonzero ``(batch_matrix, batch_coordinates)`` from *manager.matrix_generator*."""
    for _, intensity_matrix, _lengths, coordinates in manager.matrix_generator(
        batch_size=batch_size, include_mz=False,
    ):
        batch = _nonzero_rows(
            np.asarray(intensity_matrix, dtype=np.float32),
            np.asarray(coordinates)[:, :2],
        )
        if batch is not None:
            yield batch


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


def _read_spectra_zarr_sample_matrix(
    manager: MSDataManagerZarr,
    sample_ratio: float,
    *,
    random_state: int | None = None,
) -> np.ndarray:
    """Randomly sample spectra from a pixel-major Zarr store in indexed batches."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    source_indices = manager._target_indices()  # pylint: disable=protected-access
    n_spectra = int(source_indices.size)
    if n_spectra == 0:
        raise ValueError("No spectra found in the selected spatial region.")

    sample_count = max(1, int(np.ceil(n_spectra * sample_ratio)))
    rng = np.random.default_rng(random_state)
    sample_positions = np.sort(
        rng.choice(n_spectra, size=sample_count, replace=False)
    )
    sample_indices = source_indices[sample_positions]

    reader = manager._require_reader()  # pylint: disable=protected-access
    n_features = reader.n_ions
    sample_matrix = np.empty((sample_count, n_features), dtype=np.float32)
    for batch_start in range(0, sample_count, 4096):
        batch_end = min(sample_count, batch_start + 4096)
        _, intensity_flat, _lengths, _coordinates = reader.read_flat_batch(
            sample_indices[batch_start:batch_end],
            include_mz=False,
        )
        sample_matrix[batch_start:batch_end] = np.asarray(
            intensity_flat,
            dtype=np.float32,
        ).reshape(batch_end - batch_start, n_features)

    nonzero_mask = np.any(sample_matrix != 0, axis=1)
    if not np.all(nonzero_mask):
        sample_matrix = sample_matrix[nonzero_mask]
    if sample_matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")
    return sample_matrix


# --- Ion-image backends (Binary / H5 / Zarr) --------------------------------


def _read_ion_image_pixel_columns(
    manager: IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr,
    pixel_indexer: slice | np.ndarray,
    out: np.ndarray,
    ion_block: int,
) -> None:
    """Fill *out* ``(n_selected_pixels, n_ions)`` with the selected pixels' spectra.

    ``pixel_indexer`` indexes the compact pixel dimension (0..n_pixels-1).
    """
    n_ions = int(manager.n_ions)
    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        if isinstance(manager, IonImageDataManagerBinary):
            # Slice the (ion, pixel) memmap directly instead of materialising whole ion rows.
            block = manager.intensity_matrix[ion_start:ion_end, pixel_indexer]
        else:
            block = manager.get_ion_block(ion_start, ion_end)[:, pixel_indexer]
        out[:, ion_start:ion_end] = np.asarray(block).T


def _read_ion_image_sample_matrix(
    manager: IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr,
    sample_ratio: float,
    *,
    random_state: int | None = None,
    ion_block: int = 256,
) -> np.ndarray:
    """Randomly sample nonzero pixel spectra from an ion-image backend by ratio."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    n_pixels = int(manager.n_pixels)
    sample_count = max(1, int(np.ceil(n_pixels * sample_ratio)))

    rng = np.random.default_rng(random_state)
    pixel_indices = np.sort(rng.choice(n_pixels, size=sample_count, replace=False))

    sample_matrix = np.empty((sample_count, int(manager.n_ions)), dtype=np.float32)
    _read_ion_image_pixel_columns(manager, pixel_indices, sample_matrix, ion_block)
    sample_matrix = sample_matrix[np.any(sample_matrix != 0, axis=1)]

    if sample_matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")

    return sample_matrix


def _iter_ion_image_matrix_batches(
    manager: IonImageDataManagerBinary | IonImageDataManagerH5 | IonImageDataManagerZarr,
    *,
    batch_size: int = 256,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield nonzero ``(batch_matrix, batch_coordinates)`` by pixel batch from an ion-image backend."""
    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    coordinates = np.asarray(manager.coordinates)

    for pixel_start in range(0, n_pixels, batch_size):
        pixel_end = min(n_pixels, pixel_start + batch_size)
        batch_matrix = np.empty((pixel_end - pixel_start, n_ions), dtype=np.float32)
        _read_ion_image_pixel_columns(manager, slice(pixel_start, pixel_end), batch_matrix, ion_block)
        batch = _nonzero_rows(batch_matrix, coordinates[pixel_start:pixel_end, :2])
        if batch is not None:
            yield batch


# ---------------------------------------------------------------------------
# Spatial image construction
# ---------------------------------------------------------------------------

def _infer_spatial_params(
    coordinates: np.ndarray,
    *,
    width: int | None = None,
    height: int | None = None,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Shared coordinate parsing for image builders.

    Returns ``(scatter_x, scatter_y, width, height)`` with 0-based scatter indices.
    """
    scatter_x, scatter_y = zero_based_xy(coordinates)
    if width is None:
        width = int(scatter_x.max()) + 1
    if height is None:
        height = int(scatter_y.max()) + 1
    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid image dimensions: width={width}, height={height}")
    return scatter_x, scatter_y, width, height


def _scatter_to_image(
    values: np.ndarray,
    coordinates: np.ndarray,
    *,
    width: int | None,
    height: int | None,
    fill_value: int = 0,
) -> np.ndarray:
    """Place per-pixel ``values`` (n, ...) into an ``(height, width, ...)`` image."""
    scatter_x, scatter_y, w, h = _infer_spatial_params(
        coordinates, width=width, height=height,
    )
    image = np.full((h, w, *values.shape[1:]), fill_value, dtype=values.dtype)
    valid = (scatter_x >= 0) & (scatter_x < w) & (scatter_y >= 0) & (scatter_y < h)
    image[scatter_y[valid], scatter_x[valid]] = values[valid]
    return image


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

    rgb_uint8 = (np.clip(scaled_embedding, 0.0, 1.0) * 255).astype(np.uint8)
    return _scatter_to_image(rgb_uint8, coordinates, width=width, height=height)


def build_kmeans_image(
    labels: np.ndarray,
    coordinates: np.ndarray,
    *,
    n_clusters: int,
    width: int | None = None,
    height: int | None = None,
) -> np.ndarray:
    """Map cluster labels to a colour-coded spatial image."""
    colors = _tab20_rgb_palette(n_clusters)[labels]
    return _scatter_to_image(colors, coordinates, width=width, height=height)


def build_kmeans_label_image(
    labels: np.ndarray,
    coordinates: np.ndarray,
    *,
    width: int | None = None,
    height: int | None = None,
) -> np.ndarray:
    """Map cluster labels to a spatial integer-label image."""
    return _scatter_to_image(
        np.asarray(labels, dtype=np.int32),
        coordinates,
        width=width,
        height=height,
        fill_value=-1,
    )


def _tab20_label_colors(n_clusters: int) -> dict[int, str]:
    """Build the human-readable label-to-colour mapping for KMeans results."""
    return {label: TAB20_COLOR_NAMES[label] for label in range(n_clusters)}


def _tab20_rgb_palette(n_clusters: int) -> np.ndarray:
    """Return the tab20 RGB palette used by KMeans labels."""
    tab20_colors = matplotlib.colormaps["tab20"].colors  # type: ignore
    palette = (np.asarray(tab20_colors)[:, :3] * 255).astype(np.uint8)
    if n_clusters > len(palette):
        raise ValueError(f"n_clusters must be <= {len(palette)}.")
    return palette


# ---------------------------------------------------------------------------
# High-level orchestration
# ---------------------------------------------------------------------------


def _validate_umap_parameters(
    n_components: int,
    sample_ratio: float,
    transform_batch_size: int | None,
    random_state: int | None,
) -> None:
    if n_components != 3:
        raise ValueError("n_components must be 3 for RGB image output.")
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")
    if transform_batch_size is not None and transform_batch_size <= 0:
        raise ValueError("transform_batch_size must be a positive integer or None.")
    if random_state is not None:
        logger.warning(
            "Random state is set; UMAP forces single-threaded execution (n_jobs=1), "
            "so runs may be slower. Set random_state=None to allow parallelism."
        )


def _validate_parameters(
    n_components: int,
    sample_ratio: float,
    transform_batch_size: int | None,
    random_state: int | None,
    n_clusters: int,
) -> None:
    _validate_umap_parameters(
        n_components,
        sample_ratio,
        transform_batch_size,
        random_state,
    )
    if n_clusters <= 0 or n_clusters > 20:
        raise ValueError("n_clusters must be in the range [1, 20] for tab20 color mapping.")

def _is_continuous(data_manager: Any) -> bool:
    """Heuristic check for continuous/shared m/z axis."""
    if isinstance(data_manager, MSDataManagerZarr):
        return data_manager.storage_mode == ZARR_MODE_CONTINUOUS
    meta = getattr(getattr(data_manager, "ms", None), "meta", None)
    return getattr(meta, "continuous", None) is True

def _get_spatial_dimensions(data_manager: Any) -> tuple[int | None, int | None]:
    """Extract (width, height) from data_manager metadata."""
    if isinstance(data_manager, _ION_IMAGE_MANAGERS):
        return int(data_manager.width), int(data_manager.height)
    meta = getattr(getattr(data_manager, "ms", None), "meta", None)
    w = getattr(meta, "max_count_of_pixels_x", None)
    h = getattr(meta, "max_count_of_pixels_y", None)
    return (int(w) if w else None), (int(h) if h else None)


def _run_umap_stage(
    data_manager: Any,
    *,
    n_components: int,
    sample_ratio: float,
    transform_batch_size: int | None,
    umap_metric: str,
    random_state: int | None,
) -> tuple[np.ndarray, np.ndarray, int | None, int | None]:
    if (
        isinstance(data_manager, (MSDataManagerImzML, MSDataManagerZarr))
        and _is_continuous(data_manager) is False
    ):
        logger.error("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")
        raise ValueError("Data must have a continuous/shared m/z axis for UMAP/KMeans analysis.")

    width, height = _get_spatial_dimensions(data_manager)
    batch_size = 1024 if transform_batch_size is None else int(transform_batch_size)
    scaled_embedding, coordinates = compute_streaming_umap_embeddings(
        data_manager,
        n_components=n_components,
        metric=umap_metric,
        random_state=random_state,
        sample_ratio=sample_ratio,
        transform_batch_size=batch_size,
    )
    return scaled_embedding, coordinates, width, height


def _save_rgb_image(
    image: np.ndarray,
    output_path: str | Path,
    *,
    dpi: int,
    legend_handles: list[mpatches.Patch] | None = None,
    show: bool = False,
) -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    pil_kwargs: dict[str, Any] | None = None
    if output.suffix.lower() in {".jpg", ".jpeg"}:
        pil_kwargs = {"quality": 95, "subsampling": 0, "optimize": True}

    if legend_handles is None:
        plt.imsave(
            str(output),
            image,
            dpi=dpi,
            pil_kwargs=pil_kwargs,
        )
        return

    fig, ax = plt.subplots()
    ax.imshow(image)
    ax.axis("off")
    ncol = max(1, math.ceil(len(legend_handles) / 10))
    ax.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        ncol=ncol,
    )
    fig.savefig(
        str(output),
        bbox_inches="tight",
        dpi=dpi,
        pil_kwargs=pil_kwargs,
    )
    if show:
        plt.show()
    plt.close(fig)


def _save_umap_result_zarr(
    data_manager: Any,
    output_path: str | Path,
    scaled_embedding: np.ndarray,
    coordinates: np.ndarray,
) -> Path:
    """Save UMAP arrays inside the source MSI Zarr or a standalone Zarr group."""
    scatter_x, scatter_y, *_ = _infer_spatial_params(coordinates)
    zero_based_coordinates = np.column_stack((scatter_x, scatter_y)).astype(
        np.uint32,
        copy=False,
    )

    if isinstance(data_manager, MSDataManagerZarr):
        zarr_path = Path(data_manager.zarr_filepath)
        root = zarr.open_group(str(zarr_path), mode="r+")
        analysis_group = root.require_group("analysis")
        result_group = analysis_group.require_group("umap")
        overwrite = True
    else:
        zarr_path = Path(output_path).with_suffix(".zarr")
        result_group = zarr.open_group(str(zarr_path), mode="w", zarr_format=3)
        overwrite = False

    result_group.create_array(
        "scaled_embedding",
        data=scaled_embedding,
        overwrite=overwrite,
    )
    result_group.create_array(
        "coordinates",
        data=zero_based_coordinates,
        overwrite=overwrite,
    )
    return zarr_path


def plot_umap_image(
    data_manager: Any,
    *,
    output_path: str | Path = "./umap_image.jpg",
    dpi: int = 600,
    save_matrix: Literal["zarr"] | None = None,
    sample_ratio: float = 0.1,
    transform_batch_size: int | None = 1024,
    n_components: int = 3,
    umap_metric: str = "cosine",
    random_state: int | None = None,
) -> np.ndarray:
    """Run UMAP analysis, save the spatial UMAP image, and return it."""
    if save_matrix not in {None, "zarr"}:
        raise ValueError("save_matrix must be None or 'zarr'.")
    _validate_umap_parameters(
        n_components,
        sample_ratio,
        transform_batch_size,
        random_state,
    )
    scaled_embedding, coordinates, width, height = _run_umap_stage(
        data_manager,
        n_components=n_components,
        sample_ratio=sample_ratio,
        transform_batch_size=transform_batch_size,
        umap_metric=umap_metric,
        random_state=random_state,
    )
    umap_img = build_umap_image(
        scaled_embedding,
        coordinates,
        width=width,
        height=height,
    )
    _save_rgb_image(umap_img, output_path, dpi=dpi)
    logger.info(f"UMAP image saved to: {output_path}")
    if save_matrix == "zarr":
        zarr_path = _save_umap_result_zarr(
            data_manager,
            output_path,
            scaled_embedding,
            coordinates,
        )
        result_location = (
            zarr_path / "analysis" / "umap"
            if isinstance(data_manager, MSDataManagerZarr)
            else zarr_path
        )
        logger.info(f"UMAP result matrices saved to: {result_location}")
    return umap_img


def plot_umap_kmeans_image(
    data_manager: Any,
    *,
    output_dir: str | Path = "./umap_kmeans_output",
    sample_ratio: float = 0.1,
    transform_batch_size: int | None = 1024,
    dpi: int = 600,
    save_images: bool = True,
    save_matrix: Literal["npz", "zarr"] | None = None,
    n_components: int = 3,
    umap_metric: str = "cosine",
    n_clusters: int = 10,
    random_state: int | None = None,
    show: bool = False,
) -> UmapKmeansResult:
    """Run the full UMAP/KMeans analysis pipeline on processed MSI data."""

    # --- 1. Validate parameters and prepare output directory ------------
    _validate_parameters(n_components, sample_ratio, transform_batch_size, random_state, n_clusters)
    if save_matrix not in {None, "npz", "zarr"}:
        raise ValueError("save_matrix must be None, 'npz', or 'zarr'.")

    out = Path(output_dir) if output_dir is not None else Path("./umap_kmeans_output")
    out.mkdir(parents=True, exist_ok=True)

    # --- 2. UMAP embedding ---------------------------------------------
    scaled_embedding, coordinates, width, height = _run_umap_stage(
        data_manager,
        n_components=n_components,
        sample_ratio=sample_ratio,
        transform_batch_size=transform_batch_size,
        umap_metric=umap_metric,
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
    kmeans_label_img = build_kmeans_label_image(
        labels,
        coordinates,
        width=width,
        height=height,
    )
    label_colors = _tab20_label_colors(n_clusters)
    result = UmapKmeansResult(
        umap_image=umap_img,
        kmeans_image=kmeans_img,
        kmeans_label_image=kmeans_label_img,
        label_colors=label_colors,
    )

    # ---- 5. Save outputs ------------------------------------------------
    if save_images:
        _save_rgb_image(umap_img, out / "umap_image.png", dpi=dpi)
        palette = _tab20_rgb_palette(n_clusters)
        legend_handles = [
            mpatches.Patch(color=tuple(palette[i] / 255), label=str(i))
            for i in range(n_clusters)
        ]
        _save_rgb_image(kmeans_img, out / "kmeans_image.png", dpi=dpi, legend_handles=legend_handles, show=show)

    if save_matrix == "zarr":
        result_group = zarr.open_group(str(out / "umap_kmeans_result.zarr"), mode="w")
        result_group.create_array("umap_image", data=result.umap_image)
        result_group.create_array("kmeans_image", data=result.kmeans_image)
        result_group.create_array("kmeans_label_image", data=result.kmeans_label_image)
        result_group.attrs["format"] = "massflow_feature_analysis"
        result_group.attrs["label_colors"] = {
            str(label): color_name for label, color_name in result.label_colors.items()
        }
        (out / "kmeans_label_colors.json").write_text(
            json.dumps(result.label_colors, indent=2),
            encoding="utf-8",
        )
    elif save_matrix == "npz":
        color_items = sorted(result.label_colors.items())
        np.savez_compressed(
            out / "umap_kmeans_result.npz",
            umap_image=result.umap_image,
            kmeans_image=result.kmeans_image,
            kmeans_label_image=result.kmeans_label_image,
            label_color_labels=np.asarray([label for label, _ in color_items], dtype=np.int32),
            label_color_names=np.asarray([color_name for _, color_name in color_items]),
        )
        (out / "kmeans_label_colors.json").write_text(
            json.dumps(result.label_colors, indent=2),
            encoding="utf-8",
        )

    logger.info(f"UMAP/KMeans analysis complete. Output: {out}")
    result.print_label_colors()
    return result
