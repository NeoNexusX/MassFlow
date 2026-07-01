from __future__ import annotations

from typing import Any
import time
from collections.abc import Iterator
from pathlib import Path

import pytest
import zarr
import numpy as np
from PIL import Image

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor
from massflow.segmentation.umap_kmeans import plot_umap_kmeans_image
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tests.test_feature_analysis")


FILEPATH = "data/example.imzML"
ROUND = 5


def _require_imzml_pair(path: Path) -> None:
    if not path.exists() or not path.with_suffix(".ibd").exists():
        pytest.skip(f"Missing real imzML test data: {path} and matching .ibd")


@pytest.fixture(autouse=True, scope="module")
def cleanup_output_dirs():
    yield
    import shutil
    dirs = [
        "./all_benchmark",
        "./ratio_0_2_benchmark",
        "./ratio_0_4_benchmark",
        "./ratio_1_0_benchmark",
        "./validation_imzml",
        "./validation_zarr",
        "./temp_zarr_dm"
    ]
    for d in dirs:
        shutil.rmtree(d, ignore_errors=True)


def _run_feature_analysis(
    data_manager: Any,
    sample_ratio: float,
    transform_batch_size: int | None,
    output_dir: str,
    random_state: int | None = None,
) -> None:
    "Run UMAP/KMeans plotting analysis."
    plot_umap_kmeans_image(
        data_manager,
        output_dir=output_dir,
        save_images=True,
        save_matrix="zarr",
        n_clusters=10,
        sample_ratio=sample_ratio,
        transform_batch_size=transform_batch_size,
        random_state=random_state,
    )


def _validate_results_zarr(path1: str, path2: str) -> bool:
    "Validate saved result matrices."
    root_1 = zarr.open_group(f"{path1}/umap_kmeans_result.zarr", mode="r")
    root_2 = zarr.open_group(f"{path2}/umap_kmeans_result.zarr", mode="r")

    return (
        np.array_equal(root_1["umap_image"][:], root_2["umap_image"][:]) # type: ignore
        and np.array_equal(root_1["kmeans_image"][:], root_2["kmeans_image"][:]) # type: ignore
    )


def _image_pixels_equal(path1: str, path2: str) -> bool:
    "Check if two images have the same pixel values."
    with Image.open(path1) as raw_img1, Image.open(path2) as raw_img2:
        img1 = raw_img1.convert("RGBA")
        img2 = raw_img2.convert("RGBA")

        if img1.size != img2.size:
            return False

        arr1 = np.asarray(img1)
        arr2 = np.asarray(img2)

    return np.array_equal(arr1, arr2)


class TestRunUMAPKMeans:
    "test plot_umap_kmeans_image function."

    @pytest.fixture(scope="class")
    def data_manager(self) -> Iterator[MSDataManagerImzML]:
        "Fixture for data manager."
        _require_imzml_pair(Path(FILEPATH))
        data_manager = MSDataManagerImzML(filepath=FILEPATH)
        data_manager.load_head_data()
        yield data_manager

        data_manager.close()

    @pytest.mark.parametrize("sample_ratio, transform_batch_size, output_dir", [
        (0.2, 2048, "./ratio_0_2_benchmark"),
        (0.4, 2048, "./ratio_0_4_benchmark"),
        (1.0, 2048, "./ratio_1_0_benchmark"),
    ])
    @pytest.mark.parametrize("random_state", [None])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_feature_analysis(self, benchmark, data_manager, sample_ratio, transform_batch_size, output_dir, random_state) -> None:
        "Benchmark feature analysis."
        benchmark.pedantic(
            _run_feature_analysis,
            args=(data_manager, sample_ratio, transform_batch_size, output_dir, random_state),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

    @pytest.mark.parametrize("sample_ratio, transform_batch_size", [(0.2, 1024)])
    def test_feature_analysis_validation(self, data_manager, sample_ratio, transform_batch_size):
        "Validate feature analysis results."
        imzml_dm = data_manager
        zarr_dm = Preprocessor(imzml_dm, output_format="ion_image_zarr", temp_dir="./temp_zarr_dm").start()

        try:
            _run_feature_analysis(imzml_dm, sample_ratio, transform_batch_size, "./validation_imzml", 0)
            _run_feature_analysis(zarr_dm, sample_ratio, transform_batch_size, "./validation_zarr", 0)
        finally:
            zarr_dm.close()

        for image_name in ("umap_image.jpg", "kmeans_image.jpg"):
            if not _image_pixels_equal(f"./validation_imzml/{image_name}", f"./validation_zarr/{image_name}"):
                pytest.fail(f"{image_name} differs between imzML and Zarr results.")

        if not _validate_results_zarr("./validation_imzml", "./validation_zarr"):
            pytest.fail("Feature analysis results differ between imzML and Zarr results.")
