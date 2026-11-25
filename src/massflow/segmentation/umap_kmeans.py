"""UMAP/KMeans feature analysis for MSI data."""

from __future__ import annotations

import json
import math
from pathlib import Path
from time import perf_counter
from typing import Any, Iterator, Literal

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import umap
import zarr
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler

from massflow.concurrency import (
    get_allocation,
    phase,
    resolve_compute_workers,
)

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
    MSDataManagerZarr,
    ZARR_MODE_CONTINUOUS,
)
from massflow.segmentation.segmentation_pipeline import (
    KMeansInit,
    TAB20_COLOR_NAMES,
    UmapKmeansResult,
)
from massflow.tools.logger import get_logger
from massflow.tools.read_matrix import filter_target_indices

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
    n_neighbors: int,
    min_dist: float,
    spread: float,
    n_epochs: int | None,
    learning_rate: float,
    umap_init: str | np.ndarray,
    repulsion_strength: float,
    negative_sample_rate: int,
    random_state: int | None,
) -> umap.UMAP:
    """Create the UMAP reducer used by feature analysis."""
    # -1 keeps UMAP's own default (use the calling thread's numba mask); a
    # positive allocation is applied by UMAP itself during fit.
    umap_jobs = get_allocation("compute", default=-1)
    if umap_jobs > 0:
        umap_jobs = resolve_compute_workers(umap_jobs)
    return umap.UMAP(
        metric=metric,
        output_metric="euclidean",
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        n_epochs=n_epochs,
        learning_rate=learning_rate,
        # UMAP validates ndarray initial coordinates as well as method names.
        init=umap_init,  # pyright: ignore[reportArgumentType]
        repulsion_strength=repulsion_strength,
        negative_sample_rate=negative_sample_rate,
        random_state=random_state,
        n_jobs=umap_jobs,
        force_approximation_algorithm=True,
    )


@phase("umap")
def compute_streaming_umap_embeddings(
    data_manager: Any,
    *,
    n_components: int = 3,
    metric: str = "cosine",
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    spread: float = 1.0,
    n_epochs: int | None = None,
    learning_rate: float = 1.0,
    umap_init: str | np.ndarray = "spectral",
    repulsion_strength: float = 1.0,
    negative_sample_rate: int = 5,
    random_state: int | None = None,
    sample_ratio: float = 1.0,
    transform_batch_size: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute UMAP embeddings with sample-fit and batch-transform.

    ``n_epochs=None`` lets UMAP choose the training length. An array passed as
    ``umap_init`` supplies one initial embedding vector per fitted row.
    """
    reducer = _create_umap_reducer(
        n_components=n_components,
        metric=metric,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        n_epochs=n_epochs,
        learning_rate=learning_rate,
        umap_init=umap_init,
        repulsion_strength=repulsion_strength,
        negative_sample_rate=negative_sample_rate,
        random_state=random_state,
    )
    spectra_input = isinstance(data_manager, (MSDataManagerImzML, MSDataManagerZarr))
    read_started = perf_counter()
    sample_matrix, full_coordinates, sample_indices = _read_umap_sample_matrix(
        data_manager,
        sample_ratio,
        random_state=random_state,
    )

    sample_matrix = _nan_to_num_inplace_chunks(sample_matrix)
    if spectra_input:
        logger.info(f"UMAP input preparation took {perf_counter() - read_started:.2f}s.")
    logger.info(
        f"Streaming UMAP fit: sample_ratio={sample_ratio}, "
        f"sample_matrix shape={sample_matrix.shape}"
    )
    if sample_matrix.shape[0] < 2:
        raise ValueError("UMAP transform requires at least two non-zero fitted samples.")
    fit_started = perf_counter()
    reducer.fit(sample_matrix)
    if spectra_input:
        logger.info(f"UMAP fit took {perf_counter() - fit_started:.2f}s.")
    if sample_ratio == 1.0:
        if spectra_input:
            embeddings = reducer.embedding_.astype(np.float32, copy=False)
        else:
            embeddings = reducer.transform(sample_matrix).astype(np.float32, copy=False) # type: ignore
        coordinates = full_coordinates
        del sample_matrix
        return _finalize_umap_embeddings(data_manager, embeddings, coordinates) # type: ignore
    del sample_matrix

    transform_started = perf_counter()
    if random_state is None and sample_indices is None:
        enable_umap_parallel_transform(reducer)

    embeddings, coordinates = _transform_umap_batches(
        data_manager,
        reducer,
        batch_size=transform_batch_size,
        sample_indices=sample_indices,
    )
    if spectra_input:
        logger.info(
            f"UMAP transform (including index preparation) took "
            f"{perf_counter() - transform_started:.2f}s."
        )
    return _finalize_umap_embeddings(data_manager, embeddings, coordinates)


def _read_umap_sample_matrix(
    data_manager: Any,
    sample_ratio: float,
    *,
    random_state: int | None,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    """Read training spectra, optional full coordinates and fitted source IDs."""
    if sample_ratio == 1.0:
        sample_matrix, coordinates = _read_full_matrix(data_manager, batch_size=4096)
        return sample_matrix, coordinates, None
    if isinstance(data_manager, MSDataManagerImzML):
        matrix, indices = _read_imzml_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        )
        return matrix, None, indices
    if isinstance(data_manager, MSDataManagerZarr):
        matrix, indices = _read_spectra_zarr_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        )
        return matrix, None, indices
    if isinstance(data_manager, IonImageDataManagerZarr):
        return _read_zarr_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        ), None, None
    if isinstance(data_manager, IonImageDataManagerH5):
        return _read_dense_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        ), None, None
    if isinstance(data_manager, IonImageDataManagerBinary):
        return _read_binary_sample_matrix(
            data_manager, sample_ratio, random_state=random_state,
        ), None, None
    raise TypeError(f"Unsupported data_manager type: {type(data_manager).__name__}")


def _iter_umap_matrix_batches(
    data_manager: Any,
    *,
    batch_size: int,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Iterate nonzero spectra batches for UMAP transform."""
    if isinstance(data_manager, MSDataManagerImzML):
        return _iter_imzml_matrix_batches(data_manager, batch_size=batch_size)
    if isinstance(data_manager, MSDataManagerZarr):
        return _iter_spectra_zarr_matrix_batches(
            data_manager, batch_size=batch_size,
        )
    if isinstance(data_manager, IonImageDataManagerZarr):
        return _iter_zarr_matrix_batches(data_manager, batch_size=batch_size)
    if isinstance(data_manager, IonImageDataManagerH5):
        return _iter_dense_matrix_batches(data_manager, batch_size=batch_size)
    if isinstance(data_manager, IonImageDataManagerBinary):
        return _iter_binary_matrix_batches(data_manager, batch_size=batch_size)
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
    elif isinstance(data_manager, (IonImageDataManagerZarr, IonImageDataManagerH5, IonImageDataManagerBinary)):
        total_items = int(data_manager.n_pixels)
    else:
        raise TypeError(f"Unsupported data_manager type: {type(data_manager).__name__}")
    return max(1, (total_items + batch_size - 1) // batch_size)


def _transform_umap_batches(
    data_manager: Any,
    reducer: umap.UMAP,
    *,
    batch_size: int,
    sample_indices: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform nonzero spectra batches with a fitted UMAP reducer."""
    if sample_indices is not None:
        return _transform_unseen_spectra(
            data_manager, reducer, sample_indices, batch_size=batch_size,
        )
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


def _transform_unseen_spectra(
    data_manager: MSDataManagerImzML | MSDataManagerZarr,
    reducer: umap.UMAP,
    sample_indices: np.ndarray,
    *,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Reuse fitted source rows and pack only unseen spectra into query batches."""
    if isinstance(data_manager, MSDataManagerImzML):
        source_indices = filter_target_indices(
            data_manager.parser.coordinates, data_manager.target_locs,
        )
    else:
        source_indices = data_manager._target_indices()  # pylint: disable=protected-access
    sample_order = np.argsort(sample_indices)
    sorted_samples = sample_indices[sample_order]
    fit_embedding = np.asarray(reducer.embedding_, dtype=np.float32)
    embeddings = np.empty((source_indices.size, fit_embedding.shape[1]), dtype=np.float32)
    batches = data_manager.matrix_generator(batch_size=batch_size, include_mz=False)
    batch = next(batches, None)
    if batch is None:
        raise ValueError("No non-zero spectra found during transform.")
    coordinates = np.empty((source_indices.size, 2), dtype=np.asarray(batch[3]).dtype)
    query = np.empty((batch_size, np.asarray(batch[1]).shape[1]), dtype=np.float32)
    query_positions = np.empty(batch_size, dtype=np.intp)
    source_cursor = output_cursor = filled = calls = transformed_rows = reused_rows = 0
    batch_number = 0
    total_batches = max(1, math.ceil(source_indices.size / batch_size))

    def flush(count: int) -> None:
        nonlocal calls, transformed_rows
        if calls == 0 and reducer.random_state is None:
            enable_umap_parallel_transform(reducer)
        _nan_to_num_inplace_chunks(query[:count])
        embeddings[query_positions[:count]] = reducer.transform(query[:count])
        calls += 1
        transformed_rows += count

    while batch is not None:
        _, intensity_matrix, _, batch_coordinates = batch
        batch_number += 1
        intensity = np.asarray(intensity_matrix, dtype=np.float32)
        end = source_cursor + intensity.shape[0]
        indices = source_indices[source_cursor:end]
        source_cursor = end
        nonzero = np.any(intensity != 0, axis=1)
        indices = indices[nonzero]
        if not np.all(nonzero):
            intensity = intensity[nonzero]
        output_end = output_cursor + indices.size
        coordinates[output_cursor:output_end] = np.asarray(batch_coordinates)[nonzero, :2]
        positions = np.searchsorted(sorted_samples, indices)
        fitted = positions < sorted_samples.size
        fitted[fitted] = sorted_samples[positions[fitted]] == indices[fitted]
        embeddings[output_cursor + np.flatnonzero(fitted)] = fit_embedding[
            sample_order[positions[fitted]]
        ]
        reused_rows += int(np.count_nonzero(fitted))
        unseen = np.flatnonzero(~fitted)
        offset = 0
        while offset < unseen.size:
            take = min(batch_size - filled, unseen.size - offset)
            selected = unseen[offset:offset + take]
            query[filled:filled + take] = intensity[selected]
            query_positions[filled:filled + take] = output_cursor + selected
            filled += take
            offset += take
            if filled == batch_size:
                flush(filled)
                filled = 0
        output_cursor = output_end
        if batch_number % max(1, total_batches // 10) == 0 or batch_number == total_batches:
            logger.info(f"UMAP source batches: {batch_number}/{total_batches}.")
        batch = next(batches, None)
    if filled:
        flush(filled)
    logger.info(
        f"UMAP reused {reused_rows} fit rows; transformed {transformed_rows} "
        f"unseen rows in {calls} batches."
    )
    return embeddings[:output_cursor], coordinates[:output_cursor]


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
    kmeans_init: KMeansInit = "k-means++",
    n_init: int | Literal["auto"] = "auto",
    max_iter: int = 300,
    tol: float = 1e-4,
    algorithm: Literal["lloyd", "elkan"] = "lloyd",
    random_state: int | None = None,
) -> np.ndarray:
    """Perform KMeans clustering on embeddings.

    With ``n_init="auto"``, scikit-learn chooses the number of runs from
    ``kmeans_init``. The default ``"k-means++"`` uses one run.
    """
    logger.info(f"Performing KMeans clustering with {n_clusters} clusters.")

    if isinstance(kmeans_init, str) or callable(kmeans_init):
        init = kmeans_init
    else:
        init = np.asarray(kmeans_init)
    with phase("kmeans"):
        clusters = KMeans(
            n_clusters=n_clusters,
            init=init, # type: ignore
            n_init=n_init,
            max_iter=max_iter,
            tol=tol,
            algorithm=algorithm,
            random_state=random_state,
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


def _read_full_matrix(
    data_manager: Any,
    *,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a full nonzero spectra matrix and coordinate matrix from the active backend."""
    if isinstance(data_manager, (MSDataManagerImzML, MSDataManagerZarr)):
        return _read_full_spectra_matrix(data_manager, batch_size=batch_size)
    if isinstance(data_manager, IonImageDataManagerZarr):
        batch_iter = _iter_zarr_matrix_batches(
            data_manager, batch_size=batch_size,
        )
    elif isinstance(data_manager, IonImageDataManagerH5):
        batch_iter = _iter_dense_matrix_batches(
            data_manager, batch_size=batch_size,
        )
    elif isinstance(data_manager, IonImageDataManagerBinary):
        batch_iter = _iter_binary_matrix_batches(
            data_manager, batch_size=batch_size,
        )
    else:
        raise TypeError(f"Unsupported data_manager type: {type(data_manager).__name__}")

    matrix_parts: list[np.ndarray] = []
    coordinates_parts: list[np.ndarray] = []
    for batch_matrix, batch_coords in batch_iter:
        matrix_parts.append(batch_matrix)
        coordinates_parts.append(batch_coords)

    if not matrix_parts:
        raise ValueError("No non-zero spectra found.")

    return np.vstack(matrix_parts), np.vstack(coordinates_parts)


def _read_full_spectra_matrix(
    data_manager: MSDataManagerImzML | MSDataManagerZarr,
    *,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Fill one allocation instead of retaining and stacking all spectra batches."""
    if isinstance(data_manager, MSDataManagerImzML):
        # Use the same spatial selection as matrix_generator, independent of ms.
        source_indices = filter_target_indices(
            data_manager.parser.coordinates, data_manager.target_locs,
        )
    else:
        source_indices = data_manager._target_indices()  # pylint: disable=protected-access
    n_spectra = source_indices.size
    del source_indices

    batches = _iter_umap_matrix_batches(data_manager, batch_size=batch_size)
    batch = next(batches, None)
    if batch is None:
        raise ValueError("No non-zero spectra found.")
    matrix = np.empty((n_spectra, batch[0].shape[1]), dtype=np.float32)
    coordinates = np.empty((n_spectra, batch[1].shape[1]), dtype=batch[1].dtype)
    row_start = 0
    while batch is not None:
        batch_matrix, batch_coords = batch
        row_end = row_start + batch_matrix.shape[0]
        matrix[row_start:row_end] = batch_matrix
        coordinates[row_start:row_end] = batch_coords
        row_start = row_end
        batch = next(batches, None)

    # Both arrays own their buffers; no views escape before shrinking zero rows.
    if row_start < n_spectra:
        matrix.resize((row_start, matrix.shape[1]))
        coordinates.resize((row_start, coordinates.shape[1]))
    return matrix, coordinates


# --- imzML backend ---------------------------------------------------------
#
# imzML samples source pixel IDs with a shared .ibd handle, then uses the
# data manager's matrix_generator for sequential batch transforms.
def _read_imzml_sample_matrix(
    manager: MSDataManagerImzML,
    sample_ratio: float,
    *,
    random_state: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Randomly sample spectra and return their matrix and source pixel IDs."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    available_indices = filter_target_indices(manager.parser.coordinates, manager.target_locs)
    n_spectra = available_indices.size
    n_features = len(manager.ms.shared_mz_list)  # type: ignore
    sample_count = max(1, int(np.ceil(n_spectra * sample_ratio)))

    rng = np.random.default_rng(random_state)
    sample_positions = np.sort(rng.choice(n_spectra, size=sample_count, replace=False))
    sample_indices = available_indices[sample_positions]

    sample_matrix = np.empty((sample_count, n_features), dtype=np.float32)
    source_indices = np.empty(sample_count, dtype=np.int64)
    row_count = 0
    reader = manager.reader
    with open(manager.ibd_filepath, "rb") as file_handle:
        for idx in sample_indices:
            _, raw_intensity = reader.read_spectrum_from_file(file_handle, int(idx))
            intensity = np.asarray(raw_intensity, dtype=np.float32)
            if np.any(intensity != 0):
                sample_matrix[row_count] = intensity[:n_features]
                source_indices[row_count] = idx
                row_count += 1

    if row_count == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")
    if row_count < sample_count:
        sample_matrix.resize((row_count, n_features))
        source_indices.resize(row_count)
    return sample_matrix, source_indices


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

        # Direct .ibd reads may be read-only views over bytes; UMAP sanitizes in place.
        if np.all(nonzero_mask) and intensity.flags.writeable:
            yield intensity, coords
        else:
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


# --- Pixel-major spectra Zarr backend --------------------------------------


def _read_spectra_zarr_sample_matrix(
    manager: MSDataManagerZarr,
    sample_ratio: float,
    *,
    random_state: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Read sampled spectra and their source IDs from a pixel-major Zarr store."""
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
        sample_indices = sample_indices[nonzero_mask]
    if sample_matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")
    return sample_matrix, sample_indices


def _iter_spectra_zarr_matrix_batches(
    manager: MSDataManagerZarr,
    *,
    batch_size: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield nonzero spectrum matrices from a pixel-major Zarr store."""
    for _, intensity_matrix, _lengths, coordinates in manager.matrix_generator(
        batch_size=batch_size,
        include_mz=False,
    ):
        intensity = np.asarray(intensity_matrix, dtype=np.float32)
        coords = np.asarray(coordinates)[:, :2]

        nonzero_mask = np.any(intensity != 0, axis=1)
        if not np.any(nonzero_mask):
            continue

        if np.all(nonzero_mask):
            yield intensity, coords
        else:
            yield intensity[nonzero_mask], coords[nonzero_mask]


# --- Compact Zarr ion-image backend ----------------------------------------


def _read_zarr_pixel_columns(
    manager: IonImageDataManagerZarr,
    pixel_indices: np.ndarray,
    out: np.ndarray,
    ion_block: int,
) -> None:
    """Fill pixel spectra by reading compact ion-major blocks."""
    n_ions = int(manager.n_ions)
    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        out[:, ion_start:ion_end] = manager.get_ion_block(
            ion_start, ion_end
        )[:, pixel_indices].T


def _read_zarr_sample_matrix(
    manager: IonImageDataManagerZarr,
    sample_ratio: float,
    *,
    random_state: int | None = None,
    ion_block: int = 256,
) -> np.ndarray:
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")
    sample_count = max(1, int(np.ceil(manager.n_pixels * sample_ratio)))
    rng = np.random.default_rng(random_state)
    pixel_indices = np.sort(
        rng.choice(manager.n_pixels, size=sample_count, replace=False)
    )
    matrix = np.empty((sample_count, manager.n_ions), dtype=np.float32)
    _read_zarr_pixel_columns(manager, pixel_indices, matrix, ion_block)
    matrix = matrix[np.any(matrix != 0, axis=1)]
    if matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")
    return matrix


def _iter_zarr_matrix_batches(
    manager: IonImageDataManagerZarr,
    *,
    batch_size: int = 256,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    coordinates = np.asarray(manager.coordinates)
    for pixel_start in range(0, manager.n_pixels, batch_size):
        pixel_end = min(manager.n_pixels, pixel_start + batch_size)
        pixel_indices = np.arange(pixel_start, pixel_end)
        matrix = np.empty((pixel_indices.size, manager.n_ions), dtype=np.float32)
        _read_zarr_pixel_columns(manager, pixel_indices, matrix, ion_block)
        nonzero = np.any(matrix != 0, axis=1)
        if np.any(nonzero):
            yield matrix[nonzero], coordinates[pixel_start:pixel_end, :2][nonzero]


# --- Dense ion-image backend (H5) ------------------------------------------
def _read_dense_pixel_columns(
    manager: IonImageDataManagerH5,
    pixel_indices: np.ndarray,
    out: np.ndarray,
    ion_block: int,
) -> None:
    """Fill ``out`` by reading compact ion blocks through the manager API."""
    n_ions = int(manager.n_ions)

    for ion_start in range(0, n_ions, ion_block):
        ion_end = min(n_ions, ion_start + ion_block)
        block = manager.get_ion_block(ion_start, ion_end)
        out[:, ion_start:ion_end] = block[:, pixel_indices].T


def _read_dense_sample_matrix(
    manager: IonImageDataManagerH5,
    sample_ratio: float,
    *,
    random_state: int | None = None,
    ion_block: int = 256,
) -> np.ndarray:
    """Randomly sample pixels from a dense (Zarr/H5) backend by ratio."""
    if sample_ratio <= 0 or sample_ratio > 1:
        raise ValueError("sample_ratio must be in (0, 1].")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    sample_count = max(1, int(np.ceil(n_pixels * sample_ratio)))

    rng = np.random.default_rng(random_state)
    pixel_indices = np.sort(rng.choice(n_pixels, size=sample_count, replace=False))

    sample_matrix = np.empty((sample_count, n_ions), dtype=np.float32)
    _read_dense_pixel_columns(manager, pixel_indices, sample_matrix, ion_block)

    nonzero_mask = np.any(sample_matrix != 0, axis=1)
    sample_matrix = sample_matrix[nonzero_mask]

    if sample_matrix.shape[0] == 0:
        raise ValueError("All sampled spectra are zero; cannot build sample matrix.")

    return sample_matrix


def _iter_dense_matrix_batches(
    manager: IonImageDataManagerH5,
    *,
    batch_size: int = 256,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield ``(batch_matrix, batch_coordinates)`` by pixel-batch from a dense backend."""
    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    coordinates = np.asarray(manager.coordinates)

    for pixel_start in range(0, n_pixels, batch_size):
        pixel_end = min(n_pixels, pixel_start + batch_size)
        pixel_indices = np.arange(pixel_start, pixel_end)
        batch_coords = coordinates[pixel_start:pixel_end, :2]

        batch_matrix = np.empty((pixel_end - pixel_start, n_ions), dtype=np.float32)
        _read_dense_pixel_columns(manager, pixel_indices, batch_matrix, ion_block)

        nonzero_mask = np.any(batch_matrix != 0, axis=1)
        if not np.any(nonzero_mask):
            continue

        if np.all(nonzero_mask):
            yield batch_matrix, batch_coords
        else:
            yield batch_matrix[nonzero_mask], batch_coords[nonzero_mask]


# --- Binary ion-image backend ----------------------------------------------
#
# Binary stores intensity as a memmap matrix of shape (n_ions, n_pixels).
# Unlike the dense backends (Zarr/H5), pixel_indexer here is a compact
# index into the pixel dimension (0..n_pixels-1), not a flat index into
# the H*W image space.
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


def _read_binary_sample_matrix(
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


def _iter_binary_matrix_batches(
    manager: IonImageDataManagerBinary,
    *,
    batch_size: int = 256,
    ion_block: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield ``(batch_matrix, batch_coordinates)`` by pixel-batch from binary data."""
    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    n_pixels = int(manager.n_pixels)
    n_ions = int(manager.n_ions)
    coordinates = np.asarray(manager.coordinates)

    for pixel_start in range(0, n_pixels, batch_size):
        pixel_end = min(n_pixels, pixel_start + batch_size)
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
    """Map a scaled embedding to an RGBA spatial image from its first three components.

    Components beyond the third are not drawn.  A two-component embedding fills
    the red and green channels, and a one-component embedding is drawn in grayscale.
    """
    n_components = scaled_embedding.shape[1]
    if n_components < 1:
        raise ValueError(f"UMAP image requires at least one component, got {n_components}.")
    if n_components > 3:
        logger.info(f"UMAP image shows the first 3 of {n_components} components.")

    scatter_x, scatter_y, _, w, h = _infer_spatial_params(
        coordinates, width=width, height=height,
    )

    n_drawn = min(n_components, 3)
    rgb = np.zeros((scaled_embedding.shape[0], 3), dtype=scaled_embedding.dtype)
    if n_components == 1:
        rgb[...] = scaled_embedding
    else:
        rgb[:, :n_drawn] = scaled_embedding[:, :n_drawn]
    np.clip(rgb, 0.0, 1.0, out=rgb)
    rgba_uint8 = np.full((scaled_embedding.shape[0], 4), 255, dtype=np.uint8)
    rgba_uint8[:, :3] = (rgb * 255).astype(np.uint8)

    image = np.zeros((h, w, 4), dtype=np.uint8)
    valid = (scatter_x >= 0) & (scatter_x < w) & (scatter_y >= 0) & (scatter_y < h)
    image[scatter_y[valid], scatter_x[valid]] = rgba_uint8[valid]

    return image


def build_kmeans_image(
    labels: np.ndarray,
    coordinates: np.ndarray,
    *,
    n_clusters: int,
    width: int | None = None,
    height: int | None = None,
) -> np.ndarray:
    """Map cluster labels to a colour-coded RGBA spatial image."""
    scatter_x, scatter_y, _, w, h = _infer_spatial_params(
        coordinates, width=width, height=height,
    )

    palette = _tab20_rgb_palette(n_clusters)
    colors = palette[labels]
    rgba = np.full((colors.shape[0], 4), 255, dtype=np.uint8)
    rgba[:, :3] = colors

    image = np.zeros((h, w, 4), dtype=np.uint8)
    valid = (scatter_x >= 0) & (scatter_x < w) & (scatter_y >= 0) & (scatter_y < h)
    image[scatter_y[valid], scatter_x[valid]] = rgba[valid]

    return image


def build_kmeans_label_image(
    labels: np.ndarray,
    coordinates: np.ndarray,
    *,
    width: int | None = None,
    height: int | None = None,
) -> np.ndarray:
    """Map cluster labels to a spatial integer-label image."""
    scatter_x, scatter_y, _, w, h = _infer_spatial_params(
        coordinates, width=width, height=height,
    )

    image = np.full((h, w), -1, dtype=np.int32)
    valid = (scatter_x >= 0) & (scatter_x < w) & (scatter_y >= 0) & (scatter_y < h)
    image[scatter_y[valid], scatter_x[valid]] = labels[valid]

    return image


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
    if n_components % 1 != 0 or n_components < 1:
        raise ValueError(f"n_components must be a positive integer, got {n_components}.")
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
    continuous = getattr(meta, "continuous", None)
    if continuous is True:
        return True
    return False

def _get_spatial_dimensions(data_manager: Any) -> tuple[int | None, int | None]:
    """Extract (width, height) from data_manager metadata."""
    if isinstance(data_manager, (IonImageDataManagerZarr, IonImageDataManagerH5, IonImageDataManagerBinary)):
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
    n_neighbors: int,
    min_dist: float,
    spread: float,
    n_epochs: int | None,
    learning_rate: float,
    umap_init: str | np.ndarray,
    repulsion_strength: float,
    negative_sample_rate: int,
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
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        n_epochs=n_epochs,
        learning_rate=learning_rate,
        umap_init=umap_init,
        repulsion_strength=repulsion_strength,
        negative_sample_rate=negative_sample_rate,
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
        transparent=True,
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
    sample_ratio: float = 1.0,
    transform_batch_size: int | None = 4096,
    n_components: int = 3,
    umap_metric: str = "cosine",
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    spread: float = 1.0,
    n_epochs: int | None = None,
    learning_rate: float = 1.0,
    umap_init: str | np.ndarray = "spectral",
    repulsion_strength: float = 1.0,
    negative_sample_rate: int = 5,
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
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        n_epochs=n_epochs,
        learning_rate=learning_rate,
        umap_init=umap_init,
        repulsion_strength=repulsion_strength,
        negative_sample_rate=negative_sample_rate,
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
    sample_ratio: float = 1.0,
    transform_batch_size: int | None = 1024,
    dpi: int = 600,
    save_images: bool = True,
    save_matrix: Literal["npz", "zarr"] | None = None,
    n_components: int = 3,
    umap_metric: str = "cosine",
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    spread: float = 1.0,
    n_epochs: int | None = None,
    learning_rate: float = 1.0,
    umap_init: str | np.ndarray = "spectral",
    repulsion_strength: float = 1.0,
    negative_sample_rate: int = 5,
    n_clusters: int = 10,
    kmeans_init: KMeansInit = "k-means++",
    n_init: int | Literal["auto"] = "auto",
    max_iter: int = 300,
    tol: float = 1e-4,
    algorithm: Literal["lloyd", "elkan"] = "lloyd",
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
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        n_epochs=n_epochs,
        learning_rate=learning_rate,
        umap_init=umap_init,
        repulsion_strength=repulsion_strength,
        negative_sample_rate=negative_sample_rate,
        random_state=random_state,
    )

    # ---- 3. KMeans on the *same* scaled embedding -----------------------
    labels = kmeans_clustering(
        scaled_embedding,
        n_clusters=n_clusters,
        kmeans_init=kmeans_init,
        n_init=n_init,
        max_iter=max_iter,
        tol=tol,
        algorithm=algorithm,
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
