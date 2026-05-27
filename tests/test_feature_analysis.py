from __future__ import annotations

from typing import Any
import time
from collections.abc import Iterator

import pytest
import zarr
import numpy as np
from PIL import Image

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor
from massflow.tools.feature_analysis import run_analysis
from massflow.tools.logger import get_logger

logger = get_logger("massflow.tests.test_feature_analysis")


FILEPATH = "data/example.imzML"
ROUND = 1


@pytest.fixture(autouse=True, scope="module")
def cleanup_output_dirs():
    yield
    import shutil
    dirs = [
        "./all_benchmark",
        "./4096_benchmark",
        "./2000_benchmark",
        "./validation_imzml",
        "./validation_zarr",
        "./temp_zarr_dm"
    ]
    for d in dirs:
        shutil.rmtree(d, ignore_errors=True)


def _run_feature_analysis(
    data_manager: Any,
    umap_sample: int | None,
    transform_batch_size: int | None,
    output_dir: str
) -> None:
    "Run feature analysis."

    run_analysis(
        data_manager,
        output_dir=output_dir,
        save_images=True,
        save_matrix=True,
        n_clusters=10,
        random_state=0,
        umap_sample=umap_sample,
        transform_batch_size=transform_batch_size
    )


def _validate_results_zarr(path1: str, path2: str) -> bool:
    "Validate saved result matrices."
    root_1 = zarr.open_group(f"{path1}/result.zarr", mode="r")
    root_2 = zarr.open_group(f"{path2}/result.zarr", mode="r")

    return (
        np.array_equal(root_1["umap_image"][:], root_2["umap_image"][:])
        and np.array_equal(root_1["kmeans_image"][:], root_2["kmeans_image"][:])
        and np.allclose(root_1["tic_image"][:], root_2["tic_image"][:])
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
    "test run_analysis function."

    @pytest.fixture(scope="class")
    def data_manager(self) -> Iterator[MSDataManagerImzML]:
        "Fixture for data manager."
        data_manager = MSDataManagerImzML(filepath=FILEPATH)
        data_manager.load_head_data()
        yield data_manager

        data_manager.close()

    @pytest.mark.parametrize("umap_sample, transform_batch_size, output_dir", [
        (None, None, "./all_benchmark"),
        (4096, 2048, "./4096_benchmark"),
        (2000, 1024, "./2000_benchmark")
    ])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_feature_analysis(self, data_manager, benchmark, umap_sample, transform_batch_size, output_dir) -> None:
        "Benchmark feature analysis."
        benchmark.pedantic(
            _run_feature_analysis,
            args=(data_manager, umap_sample, transform_batch_size, output_dir),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
        )

    @pytest.mark.parametrize("umap_sample, transform_batch_size", [(2000, 1024)])
    def test_feature_analysis_validation(self, data_manager, umap_sample, transform_batch_size):
        "Validate feature analysis results."
        imzml_dm = data_manager
        zarr_dm = Preprocessor(imzml_dm, output_format="ion_image_zarr", temp_dir="./temp_zarr_dm").start()

        try:
            _run_feature_analysis(imzml_dm, umap_sample, transform_batch_size, "./validation_imzml")
            _run_feature_analysis(zarr_dm, umap_sample, transform_batch_size, "./validation_zarr")
        finally:
            zarr_dm.close()

        for image_name in ("umap_image.jpg", "kmeans_image.jpg", "tic_image.jpg"):
            if not _image_pixels_equal(f"./validation_imzml/{image_name}", f"./validation_zarr/{image_name}"):
                pytest.fail(f"{image_name} differs between imzML and Zarr results.")

        if not _validate_results_zarr("./validation_imzml", "./validation_zarr"):
            pytest.fail("Feature analysis results differ between imzML and Zarr results.")
