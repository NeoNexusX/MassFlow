from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np
import umap
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler

from massflow.data_manager import MSDataManagerImzML

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FILEPATH = "tests/temp/processed_data.imzML"
OUTPUT_DIR = "tests/temp/umap_kmeans"


def imzml_to_matrix(
    data_manager: MSDataManagerImzML,
    *,
    batch_size: int = 512,
    max_threads: int = 2,
    max_features: int = 1000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Return X[pixel, mz_feature], coords, mz_axis, feature_indices, width, height."""
    mz_axis, keep_pixels, feature_indices = _scan_data(
        data_manager,
        batch_size=batch_size,
        max_threads=max_threads,
        max_features=max_features,
    )

    n_pixels = int(np.sum(keep_pixels))
    data = np.empty((n_pixels, feature_indices.size), dtype=np.float32)
    coords_out = np.empty((n_pixels, 2), dtype=np.int64)

    in_start = 0
    out_start = 0
    for batch_mz, matrix, lengths, coords in data_manager.matrix_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    ):
        matrix = _trim_matrix(matrix, lengths)
        _check_mz_axis(mz_axis, batch_mz)

        keep = keep_pixels[in_start : in_start + matrix.shape[0]]
        rows = matrix[keep][:, feature_indices]
        out_end = out_start + rows.shape[0]
        data[out_start:out_end] = rows.astype(np.float32, copy=False)
        coords_out[out_start:out_end] = np.asarray(coords)[keep, :2]
        in_start += matrix.shape[0]
        out_start = out_end

    width, height = _image_size(data_manager, coords_out)
    return data, coords_out, mz_axis[feature_indices], feature_indices, width, height


def compute_umap_embeddings(
    matrix: np.ndarray,
    *,
    n_components: int = 3,
    metric: str = "cosine",
    random_state: int = 0,
) -> np.ndarray:
    matrix = np.nan_to_num(matrix.astype(np.float32, copy=False), copy=False)
    embeddings = umap.UMAP(
        metric=metric,
        n_components=n_components,
        random_state=random_state,
    ).fit_transform(matrix)
    return MinMaxScaler().fit_transform(embeddings) # type: ignore[return-value]


def kmeans_clustering(embeddings: np.ndarray, *, n_clusters: int = 10, random_state: int = 0) -> np.ndarray:
    return KMeans(n_clusters=n_clusters, random_state=random_state, n_init="auto").fit_predict(embeddings)


def plot_results(
    embeddings: np.ndarray,
    clusters: np.ndarray,
    coords: np.ndarray,
    width: int,
    height: int,
    *,
    output_dir: str | Path = OUTPUT_DIR,
    suffix: str = "",
) -> dict[str, str]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    umap_img = np.zeros((height, width, 3), dtype=np.float32)
    cluster_img = np.full((height, width), -1, dtype=np.int32)

    xy = _zero_based_xy(coords)
    umap_img[xy[:, 1], xy[:, 0]] = embeddings[:, :3]
    cluster_img[xy[:, 1], xy[:, 0]] = clusters
    umap_img = np.clip(umap_img, 0, 1)

    umap_path = output_dir / "umap_rgb.png"
    cluster_path = output_dir / f"kmeans_clusters{suffix}.png"
    combined_path = output_dir / f"umap_kmeans_combined{suffix}.png"

    if not umap_path.exists():
        plt.imsave(umap_path, umap_img)
    plt.imsave(cluster_path, np.ma.masked_where(cluster_img < 0, cluster_img), cmap="tab20")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=200)
    axes[0].imshow(umap_img)
    axes[0].set_title("UMAP")
    axes[1].imshow(np.ma.masked_where(cluster_img < 0, cluster_img), cmap="tab20")
    axes[1].set_title("KMeans")
    for ax in axes:
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(combined_path, bbox_inches="tight")
    plt.close(fig)

    return {
        "umap_rgb": str(umap_path),
        "kmeans_clusters": str(cluster_path),
        "combined": str(combined_path),
    }


def _scan_data(
    data_manager: MSDataManagerImzML,
    *,
    batch_size: int,
    max_threads: int,
    max_features: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mz_axis = None
    coords_list = []
    nonzero_list = []
    feature_sum = None
    feature_sumsq = None
    n_rows = 0

    for batch_mz, matrix, lengths, coords in data_manager.matrix_generator(
        batch_size=batch_size,
        include_mz=True,
        max_threads=max_threads,
    ):
        matrix = _trim_matrix(matrix, lengths)
        batch_mz = np.asarray(batch_mz)
        if mz_axis is None:
            mz_axis = batch_mz
            feature_sum = np.zeros(matrix.shape[1], dtype=np.float64)
            feature_sumsq = np.zeros(matrix.shape[1], dtype=np.float64)
        else:
            _check_mz_axis(mz_axis, batch_mz)

        feature_sum += matrix.sum(axis=0, dtype=np.float64)
        feature_sumsq += np.einsum("ij,ij->j", matrix, matrix, dtype=np.float64)
        n_rows += matrix.shape[0]
        nonzero_list.append(np.any(matrix != 0, axis=1))
        coords_list.append(np.asarray(coords)[:, :3])

    if mz_axis is None or feature_sum is None or feature_sumsq is None:
        raise ValueError("No spectra were read from imzML.")

    coords_all = np.vstack(coords_list)
    keep_pixels = np.concatenate(nonzero_list)
    keep_pixels &= _keep_last_duplicate_coord(coords_all)

    feature_indices = np.arange(mz_axis.size)
    if 0 < max_features < mz_axis.size:
        mean = feature_sum / n_rows
        variance = np.maximum(feature_sumsq / n_rows - mean * mean, 0)
        feature_indices = np.sort(np.argpartition(variance, -max_features)[-max_features:])

    return mz_axis, keep_pixels, feature_indices.astype(np.intp)


def _trim_matrix(matrix: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    lengths = np.asarray(lengths)
    if np.min(lengths) != np.max(lengths):
        raise ValueError("UMAP requires aligned spectra with equal length.")
    return np.asarray(matrix)[:, : int(lengths[0])]


def _check_mz_axis(expected: np.ndarray, actual: np.ndarray) -> None:
    actual = np.asarray(actual)
    if expected.shape != actual.shape or not np.allclose(expected, actual, rtol=1e-7, atol=1e-8):
        raise ValueError("All batches must share the same m/z axis.")


def _keep_last_duplicate_coord(coords: np.ndarray) -> np.ndarray:
    _, reversed_indices = np.unique(coords[::-1], axis=0, return_index=True)
    keep_indices = coords.shape[0] - 1 - reversed_indices
    keep = np.zeros(coords.shape[0], dtype=bool)
    keep[keep_indices] = True
    return keep


def _image_size(data_manager: MSDataManagerImzML, coords: np.ndarray) -> tuple[int, int]:
    meta = data_manager.ms.meta
    width = getattr(meta, "max_count_of_pixels_x", None)
    height = getattr(meta, "max_count_of_pixels_y", None)
    if width and height:
        return int(width), int(height)

    xy = _zero_based_xy(coords)
    return int(xy[:, 0].max()) + 1, int(xy[:, 1].max()) + 1


def _zero_based_xy(coords: np.ndarray) -> np.ndarray:
    xy = np.asarray(coords)[:, :2].astype(np.int64, copy=False)
    return xy if np.any(xy == 0) else xy - 1


def main(
    *,
    input_path: str | Path = FILEPATH,
    output_dir: str | Path = OUTPUT_DIR,
    batch_size: int = 512,
    max_threads: int = 2,
    max_features: int = 0,
    n_clusters_1: int = 5,
    n_clusters_2: int = 10,
    random_state: int = 0,
) -> None:
    with MSDataManagerImzML(filepath=input_path) as dm:
        dm.load_head_data()
        matrix, coords, mz_axis, feature_indices, width, height = imzml_to_matrix(
            dm,
            batch_size=batch_size,
            max_threads=max_threads,
            max_features=max_features,
        )

    embeddings = compute_umap_embeddings(matrix, random_state=random_state)
    clusters_1 = kmeans_clustering(embeddings, n_clusters=n_clusters_1, random_state=random_state)
    clusters_2 = kmeans_clustering(embeddings, n_clusters=n_clusters_2, random_state=random_state)

    outputs_1 = plot_results(embeddings, clusters_1, coords, width, height, output_dir=output_dir, suffix=f"_k{n_clusters_1}")
    outputs_2 = plot_results(embeddings, clusters_2, coords, width, height, output_dir=output_dir, suffix=f"_k{n_clusters_2}")

    outputs = {
        "umap_rgb": outputs_1["umap_rgb"],
        f"kmeans_clusters_k{n_clusters_1}": outputs_1["kmeans_clusters"],
        f"combined_k{n_clusters_1}": outputs_1["combined"],
        f"kmeans_clusters_k{n_clusters_2}": outputs_2["kmeans_clusters"],
        f"combined_k{n_clusters_2}": outputs_2["combined"],
    }

    output_dir_path = Path(output_dir)
    np.save(output_dir_path / "umap_embeddings.npy", embeddings)
    np.save(output_dir_path / f"kmeans_labels_k{n_clusters_1}.npy", clusters_1)
    np.save(output_dir_path / f"kmeans_labels_k{n_clusters_2}.npy", clusters_2)
    np.save(output_dir_path / "coordinates.npy", coords)
    np.save(output_dir_path / "mz_axis.npy", mz_axis)
    np.save(output_dir_path / "feature_indices.npy", feature_indices)
    (output_dir_path / "metadata.json").write_text(
        json.dumps(
            {
                "matrix_shape": list(matrix.shape),
                "image_shape": [height, width],
                "n_clusters": [n_clusters_1, n_clusters_2],
                "outputs": outputs,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
