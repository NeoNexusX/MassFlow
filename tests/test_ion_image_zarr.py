"""Zarr 离子图像转写 + PyTorch 读取性能测试。"""

import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest
import torch
import zarr
from torch.utils.data import DataLoader

from massflow.data_manager import IonImageDataManagerZarr, MSDataManagerImzML
from massflow.module.ms_meta_data_zarr import MetaDataZarr
from massflow.morch import IonImageZarrDataset
from massflow.segmentation.region_mask import apply_region_mask
from massflow.tools.zarr_io import MSIZarrReader, prepare_coordinate_index

EXAMPLE_IMZML = Path("/Users/dre/Desktop/data/Example_read/example.imzML")
EXAMPLE_ZARR_DIR = Path("data")
ROUND = int(os.environ.get("MASSFLOW_BENCH_ROUNDS", "1"))
_PURGE_PASSWORD_ENV = "MASSFLOW_PURGE_PASSWORD"


def _purge_cache() -> None:
    """Optionally purge system disk cache when MASSFLOW_PURGE_CACHE_PASSWORD is set."""
    sudo_password = os.environ.get(_PURGE_PASSWORD_ENV)
    if not sudo_password:
        return
    try:
        subprocess.run(
            ["sudo", "-S", "purge"],
            input=sudo_password,
            text=True,
            check=True,
            capture_output=True,
            timeout=30,
        )
        time.sleep(5)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass


def _remove_zarr(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)


def _example_zarr_path(chunk_size: int) -> Path:
    return EXAMPLE_ZARR_DIR / f"example_chunk{int(chunk_size)}.zarr"


def _zarr_has_chunk_size(path: Path, chunk_size: int) -> bool:
    del chunk_size
    if not path.exists():
        return False
    try:
        reader = MSIZarrReader(path)
        try:
            return reader.row_axis == "ion"
        finally:
            reader.close()
    except Exception:
        return False


def _cache_example_flat_batches(*, batch_size: int = 256):
    """从 example.imzML 读取全部像素，返回内存中的 flat batches 缓存。"""
    if not EXAMPLE_IMZML.exists():
        pytest.skip(f"Example data not found at {EXAMPLE_IMZML}")

    with MSDataManagerImzML(filepath=str(EXAMPLE_IMZML)) as dm:
        dm.load_head_data()
        if len(dm.ms) == 0:
            pytest.skip("No spectra in example.imzML")

        batches = []
        for mz_axis, intensity_flat, lengths, coordinates in dm.flat_generator(
            batch_size=batch_size,
            include_mz=True,
            max_threads=2,
        ):
            batches.append(
                (
                    np.asarray(mz_axis).copy(),
                    np.asarray(intensity_flat).copy(),
                    np.asarray(lengths).copy(),
                    np.asarray(coordinates).copy(),
                )
            )

        assert batches
        return {
            "batches": batches,
            "width": int(dm.ms.meta.max_count_of_pixels_x or 0),
            "height": int(dm.ms.meta.max_count_of_pixels_y or 0),
            "mz_dtype": dm.mz_dtype,
            "intensity_dtype": dm.intensity_dtype,
        }


@pytest.fixture(scope="class")
def flat_cache():
    """class 级 fixture：一次性缓存全部 flat 数据到内存。"""
    return _cache_example_flat_batches()


@pytest.fixture
def zarr_path(request, chunk_size):
    """Return a reusable example Zarr store with the requested ion chunk size."""
    chunk_size = int(chunk_size)
    output_path = _example_zarr_path(chunk_size)
    if _zarr_has_chunk_size(output_path, chunk_size):
        return output_path
    _remove_zarr(output_path)

    flat_cache = request.getfixturevalue("flat_cache")
    zarr_dm = IonImageDataManagerZarr(
        str(output_path),
        mode="w",
        width=flat_cache["width"],
        height=flat_cache["height"],
        mz_dtype=flat_cache["mz_dtype"],
        intensity_dtype=flat_cache["intensity_dtype"],
        temp_dir="./temp",
        temp_workers=4,
        merge_workers=4,
        ion_block=512,
        chunk_size=chunk_size,
    )
    try:
        for mz_axis, intensity_flat, lengths, coordinates in flat_cache["batches"]:
            zarr_dm.write_flat_batch(
                intensity_flat=intensity_flat,
                lengths=lengths,
                coordinates=coordinates,
                mz_axis=mz_axis,
            )
        zarr_dm.finalize()
    finally:
        zarr_dm.close()
    return output_path


def _write_zarr_from_cache(
    flat_cache: dict,
    output_path: Path,
    tmp_path: Path,
    *,
    chunk_size: int = 1,
) -> None:
    """从内存缓存写出 Zarr store（不含 imzML 读取开销）。"""
    _remove_zarr(output_path)
    zarr_dm = IonImageDataManagerZarr(
        str(output_path),
        mode="w",
        width=flat_cache["width"],
        height=flat_cache["height"],
        mz_dtype=flat_cache["mz_dtype"],
        intensity_dtype=flat_cache["intensity_dtype"],
        temp_dir=str(tmp_path),
        temp_workers=4,
        merge_workers=4,
        ion_block=512,
        chunk_size=chunk_size,
    )
    try:
        for mz_axis, intensity_flat, lengths, coordinates in flat_cache["batches"]:
            zarr_dm.write_flat_batch(
                intensity_flat=intensity_flat,
                lengths=lengths,
                coordinates=coordinates,
                mz_axis=mz_axis,
            )
        zarr_dm.finalize()
    finally:
        zarr_dm.close()


def _write_tiny_zarr(
    path: Path,
    tmp_path: Path,
    coordinates,
    *,
    width: int | None = 2,
    height: int | None = 2,
    chunk_size: int = 1,
    compute_mean_spectrum: bool = False,
    metadata=None,
) -> IonImageDataManagerZarr:
    """Write a tiny two-ion Zarr store for correctness tests."""
    coords = np.asarray(coordinates, dtype=np.int64)
    n_pixels = coords.shape[0]
    pixel_values = np.arange(1, n_pixels + 1, dtype=np.float32)
    intensity = np.column_stack((pixel_values, pixel_values * 10)).ravel()

    dm = IonImageDataManagerZarr(
        str(path),
        mode="w",
        width=width,
        height=height,
        temp_dir=str(tmp_path),
        intensity_dtype=np.float32,
        ion_block=16,
        chunk_size=chunk_size,
        compute_mean_spectrum=compute_mean_spectrum,
        metadata=metadata,
        shared_mz_axis=(
            np.array([100.0, 200.0], dtype=np.float64) if metadata is not None else None
        ),
    )
    try:
        dm.write_flat_batch(
            intensity_flat=intensity,
            lengths=np.full(n_pixels, 2, dtype=np.uint32),
            coordinates=coords,
            mz_axis=np.array([100.0, 200.0], dtype=np.float64),
        )
        dm.finalize()
        return dm
    except BaseException:
        dm.close()
        raise


class TestZarrConvert:
    """Zarr 转写纯净计时（数据已在内存中）。"""

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_zarr_write_benchmark(self, benchmark, flat_cache, tmp_path):
        """benchmark：内存 → Zarr 完整转写时间（Phase 1 + Phase 2）。"""
        output = tmp_path / "bench.zarr"

        def convert():
            _write_zarr_from_cache(flat_cache, output, tmp_path)

        benchmark.pedantic(
            convert,
            rounds=ROUND,
            iterations=1,
            warmup_rounds=0,
            setup=_purge_cache,
        )
        assert output.exists()
        _remove_zarr(output)


class TestZarrDataLoader:
    """Zarr PyTorch DataLoader 全量读取 benchmark。"""

    BATCH_SIZE = 64

    @pytest.fixture
    def zarr_dataset(self, zarr_path):
        """创建 Zarr 数据集。"""
        dm = IonImageDataManagerZarr(str(zarr_path), mode="r")
        ds = IonImageZarrDataset(dm)
        yield ds
        dm.close()

    @staticmethod
    def _run_epoch(dataset, num_workers: int, shuffle: bool) -> int:
        loader = DataLoader(
            dataset.index_dataset(),
            batch_size=TestZarrDataLoader.BATCH_SIZE,
            shuffle=shuffle,
            num_workers=num_workers,
            pin_memory=False,
            collate_fn=dataset.get_batch,
        )
        total = 0
        for batch in loader:
            total += batch.size(0)
        return total

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("num_workers", [0])
    @pytest.mark.parametrize("chunk_size", [16])
    @pytest.mark.parametrize("shuffle", [False])
    def test_dataloader_epoch(self, benchmark, zarr_dataset, num_workers, chunk_size, shuffle):
        """全量 epoch 读取时间。"""
        total = benchmark.pedantic(
            lambda: self._run_epoch(zarr_dataset, num_workers, shuffle),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
            setup=_purge_cache,
        )
        assert total == len(zarr_dataset)


class TestZarrCorrectness:
    """Small synthetic checks for Zarr coordinate semantics."""

    def test_msi_zarr_v1_layout_and_metadata(self, tmp_path):
        path = tmp_path / "layout.zarr"
        meta = MetaDataZarr(
            name="layout",
            max_count_of_pixels_x=2,
            max_count_of_pixels_y=2,
            pixel_size_x=50.0,
            pixel_size_y=50.0,
        )
        meta.continuous = True
        meta.centroid_spectrum = True

        dm = _write_tiny_zarr(
            path,
            tmp_path,
            [(1, 1), (2, 1)],
            compute_mean_spectrum=True,
            metadata=meta,
        )
        dm.close()

        root = zarr.open_group(str(path), mode="r")
        assert dict(root.attrs) == {
            "format": "massflow.msi_zarr",
            "format_version": "1.0",
            "write_state": "complete",
            "coordinate_order": ["x", "y", "z"],
            "coordinate_base": 0,
            "spatial_shape": [2, 2],
        }
        assert dict(root["data"].attrs) == {
            "row_axis": "ion",
            "encoding": "continuous",
        }
        assert "axes/coordinates" in root
        assert "axes/mz" in root
        assert "data/intensity" in root
        assert "data/offsets" in root
        assert "stats/mean_spectrum" in root
        assert "metadata" in root
        np.testing.assert_array_equal(
            root["axes/coordinates"][:], [[0, 0, 0], [1, 0, 0]]
        )
        np.testing.assert_array_equal(root["data/offsets"][:], [0, 2, 4])
        np.testing.assert_array_equal(root["data/intensity"][:], [1, 2, 10, 20])
        assert root["metadata"].attrs["name"] == "layout"
        assert root["metadata"].attrs["pixel_size_x"] == 50.0

    def test_explicit_zero_based_subset_is_not_shifted(self):
        result = prepare_coordinate_index(
            np.array([[2, 1, 0]], dtype=np.uint32),
            width=3,
            height=2,
            coordinate_base=0,
        )
        np.testing.assert_array_equal(result.coordinates, [[2, 1, 0]])

    def test_reader_rejects_float_offsets(self, tmp_path):
        path = tmp_path / "float_offsets.zarr"
        dm = _write_tiny_zarr(path, tmp_path, [(1, 1), (2, 1)])
        dm.close()
        root = zarr.open_group(str(path), mode="a")
        del root["data/offsets"]
        root["data"].create_array("offsets", data=np.array([0.0, 2.0, 4.0]))
        with pytest.raises(ValueError, match="integer dtype"):
            MSIZarrReader(path)

    def test_mean_spectrum_uses_valid_pixels_and_roundtrips(self, tmp_path):
        path = tmp_path / "mean.zarr"
        dm = _write_tiny_zarr(
            path,
            tmp_path,
            [(1, 1), (2, 1)],
            compute_mean_spectrum=True,
        )
        try:
            np.testing.assert_allclose(dm.mean_spectrum, np.array([1.5, 15.0])) # type: ignore
        finally:
            dm.close()

        with IonImageDataManagerZarr(str(path), mode="r") as reopened:
            np.testing.assert_allclose(reopened.mean_spectrum, np.array([1.5, 15.0])) # type: ignore

    def test_mean_spectrum_is_absent_by_default(self, tmp_path):
        dm = _write_tiny_zarr(tmp_path / "no_mean.zarr", tmp_path, [(1, 1), (2, 1)])
        try:
            assert dm.mean_spectrum is None
        finally:
            dm.close()

    def test_one_based_coordinates_roundtrip(self, tmp_path):
        path = tmp_path / "one.zarr"
        dm = _write_tiny_zarr(path, tmp_path, [(1, 1), (2, 1)])
        try:
            assert dm.width == 2
            assert dm.height == 2
            assert dm.coordinate_base == 0
            np.testing.assert_array_equal(
                dm.coordinates, np.array([[0, 0, 0], [1, 0, 0]])
            )
            np.testing.assert_array_equal(
                dm.get_ion_image(0),
                np.array([[1, 2], [0, 0]], dtype=np.float32),
            )
            np.testing.assert_array_equal(dm.get_ion_vector(0), np.array([1, 2], dtype=np.float32))
        finally:
            dm.close()

        with IonImageDataManagerZarr(str(path), mode="r") as reopened:
            np.testing.assert_array_equal(
                reopened.get_ion_image(0),
                np.array([[1, 2], [0, 0]], dtype=np.float32),
            )
            np.testing.assert_array_equal(reopened.get_ion_vector(0), np.array([1, 2], dtype=np.float32))

    def test_zero_based_coordinates_infer_shape(self, tmp_path):
        dm = _write_tiny_zarr(
            tmp_path / "zero.zarr",
            tmp_path,
            [(0, 0), (1, 0)],
            width=None,
            height=None,
        )
        try:
            assert dm.width == 2
            assert dm.height == 1
            assert dm.coordinate_base == 0
            np.testing.assert_array_equal(dm.get_ion_image(0), np.array([[1, 2]], dtype=np.float32))
            np.testing.assert_array_equal(dm.get_ion_vector(0), np.array([1, 2], dtype=np.float32))
        finally:
            dm.close()

    def test_auto_detects_zero_based_y_axis(self, tmp_path):
        dm = _write_tiny_zarr(
            tmp_path / "y_zero.zarr",
            tmp_path,
            [(1, 0), (2, 0)],
            width=3,
            height=1,
        )
        try:
            assert dm.coordinate_base == 0
            np.testing.assert_array_equal(dm.get_ion_image(0), np.array([[0, 1, 2]], dtype=np.float32))
            np.testing.assert_array_equal(dm.get_ion_vector(0), np.array([1, 2], dtype=np.float32))
        finally:
            dm.close()

    def test_out_of_bounds_coordinates_fail_fast(self, tmp_path):
        with pytest.raises(ValueError, match="out-of-bounds"):
            _write_tiny_zarr(tmp_path / "bad.zarr", tmp_path, [(1, 1), (3, 1)], width=2, height=2)

    def test_duplicate_coordinates_keep_last_pixel(self, tmp_path):
        path = tmp_path / "duplicate.zarr"
        dm = _write_tiny_zarr(
            path,
            tmp_path,
            [(1, 1), (1, 1)],
            compute_mean_spectrum=True,
        )
        try:
            assert dm.n_pixels == 1
            np.testing.assert_array_equal(dm.coordinates, [[0, 0, 0]])
            np.testing.assert_array_equal(
                dm.get_ion_image(0),
                np.array([[2, 0], [0, 0]], dtype=np.float32),
            )
            np.testing.assert_array_equal(dm.get_ion_vector(0), [2])
            np.testing.assert_allclose(dm.mean_spectrum, [2.0, 20.0])
        finally:
            dm.close()

        with IonImageDataManagerZarr(str(path), mode="r") as reopened:
            assert reopened.n_pixels == 1
            np.testing.assert_array_equal(reopened.coordinates, [[0, 0, 0]])
            np.testing.assert_array_equal(reopened.get_ion_vector(0), [2])
            np.testing.assert_allclose(reopened.mean_spectrum, [2.0, 20.0])

    def test_dataset_batch_matches_single_item_stack(self, tmp_path):
        dm = _write_tiny_zarr(
            tmp_path / "batch.zarr",
            tmp_path,
            [(1, 1), (2, 1)],
            chunk_size=2,
        )
        try:
            dataset = IonImageZarrDataset(dm)
            indices = [1, 0, -1]
            batched = dataset.get_batch(indices)
            legacy = torch.stack([dataset[i] for i in indices], dim=0)
            assert torch.equal(batched, legacy)
        finally:
            dm.close()

    def test_region_mask_uses_stored_zero_based_coordinates(self, tmp_path):
        source = _write_tiny_zarr(
            tmp_path / "mask_source.zarr",
            tmp_path,
            [(2, 2), (3, 2)],
            width=3,
            height=2,
        )
        output = None
        try:
            mask = np.zeros((2, 3), dtype=bool)
            mask[1, 2] = True
            output = apply_region_mask(
                source,
                mask,
                tmp_path / "mask_output.zarr",
            )
            np.testing.assert_array_equal(output.coordinates, [[2, 1, 0]])
            np.testing.assert_array_equal(
                output.get_ion_image(0),
                np.array([[0, 0, 0], [0, 0, 2]], dtype=np.float32),
            )
        finally:
            if output is not None:
                output.close()
            source.close()
