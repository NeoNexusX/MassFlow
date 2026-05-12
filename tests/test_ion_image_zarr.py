"""Zarr 离子图像转写 + PyTorch 读取性能测试。"""

import json
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest
import torch
from torch.utils.data import DataLoader

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerZarr,
    MSDataManagerImzML,
)
from massflow.morch import IonImageZarrDataset

EXAMPLE_IMZML = Path("/Users/dre/Desktop/data/Example_read/example.imzML")
EXAMPLE_IIB = Path("data/example.iib")
EXAMPLE_ZARR = Path("data/example.zarr")
_SUDO_PASSWORD = "660789"
ROUND = 1


def _purge_cache() -> None:
    """清除系统磁盘缓存，随后暂停 5s 让系统稳定。"""
    try:
        subprocess.run(
            ["sudo", "-S", "purge"],
            input=_SUDO_PASSWORD,
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


@pytest.fixture(scope="class")
def zarr_path(request):
    """class 级 fixture：生成 example.zarr，后续测试复用。"""
    if EXAMPLE_ZARR.exists():
        return EXAMPLE_ZARR

    flat_cache = request.getfixturevalue("flat_cache")
    zarr_dm = IonImageDataManagerZarr(
        str(EXAMPLE_ZARR),
        mode="w",
        width=flat_cache["width"],
        height=flat_cache["height"],
        mz_dtype=flat_cache["mz_dtype"],
        intensity_dtype=flat_cache["intensity_dtype"],
        temp_dir="./temp",
        temp_workers=4,
        merge_workers=4,
        ion_block=512,
        chunk_size=1,
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
    return EXAMPLE_ZARR


def _write_zarr_from_cache(flat_cache: dict, output_path: Path, tmp_path: Path) -> None:
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
        chunk_size=16,
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

    def test_zarr_print_phase_timing(self, flat_cache, tmp_path):
        """打印 Phase 1 / Phase 2 分别耗时（不含 benchmark 框架）。"""
        output = tmp_path / "phase_timing.zarr"
        _remove_zarr(output)

        zarr_dm = IonImageDataManagerZarr(
            str(output),
            mode="w",
            width=flat_cache["width"],
            height=flat_cache["height"],
            mz_dtype=flat_cache["mz_dtype"],
            intensity_dtype=flat_cache["intensity_dtype"],
            temp_dir=str(tmp_path),
            temp_workers=4,
            merge_workers=4,
            ion_block=512,
            chunk_size=16,
        )

        t0 = time.perf_counter()
        for mz_axis, intensity_flat, lengths, coordinates in flat_cache["batches"]:
            zarr_dm.write_flat_batch(
                intensity_flat=intensity_flat,
                lengths=lengths,
                coordinates=coordinates,
                mz_axis=mz_axis,
            )
        t1 = time.perf_counter()
        zarr_dm.finalize()
        t2 = time.perf_counter()

        print(f"\n  像素数:     {zarr_dm.n_pixels}")
        print(f"  离子数:     {zarr_dm.n_ions}")
        print(f"  图像尺寸:   {zarr_dm.width} x {zarr_dm.height}")
        print(f"  Phase 1 (temp 写入): {t1 - t0:.3f}s")
        print(f"  Phase 2 (Zarr 合并): {t2 - t1:.3f}s")
        print(f"  Zarr 大小: {_store_size_mb(output):.0f} MB")
        print(f"  总转写:     {t2 - t0:.3f}s\n")
        zarr_dm.close()
        _remove_zarr(output)


class TestZarrRead:
    """Zarr 读取测试（基于已生成的 example.zarr fixture）。"""

    def test_zarr_read_performance(self, zarr_path):
        """随机 + 顺序各读取 1000 张离子图像。"""
        dm = IonImageDataManagerZarr(str(zarr_path), mode="r")
        n_ions = dm.n_ions
        n_reads = min(1000, n_ions)
        rng = np.random.default_rng(20260512)
        ion_indices = rng.integers(0, n_ions, size=n_reads, endpoint=False)

        _purge_cache()
        t0 = time.perf_counter()
        checksum = 0.0
        for idx in ion_indices:
            checksum += float(np.sum(dm.get_ion_image(int(idx))))
        t1 = time.perf_counter()
        print(f"\n  Zarr 随机 {n_reads} 张: {t1 - t0:.3f}s ({(t1 - t0) / n_reads * 1000:.2f}ms/张)")
        assert checksum > 0

        _purge_cache()
        t0 = time.perf_counter()
        checksum = 0.0
        for i in range(n_reads):
            checksum += float(np.sum(dm.get_ion_image(i)))
        t1 = time.perf_counter()
        print(f"  Zarr 顺序 {n_reads} 张: {t1 - t0:.3f}s ({(t1 - t0) / n_reads * 1000:.2f}ms/张)")
        assert checksum > 0
        dm.close()

    def test_correctness_vs_binary(self, zarr_path):
        """逐 ion 对比 Zarr 与 binary 文件的图像数据一致性。"""
        if not EXAMPLE_IIB.exists():
            pytest.skip(f"Binary reference not found at {EXAMPLE_IIB}")

        bin_dm = IonImageDataManagerBinary(str(EXAMPLE_IIB), mode="r")
        zarr_dm = IonImageDataManagerZarr(str(zarr_path), mode="r")

        assert bin_dm.n_ions == zarr_dm.n_ions
        assert bin_dm.width == zarr_dm.width and bin_dm.height == zarr_dm.height
        np.testing.assert_allclose(bin_dm.mz_axis, zarr_dm.mz_axis, rtol=1e-7, atol=1e-10)

        rng = np.random.default_rng(20260512)
        sample_ions = rng.integers(0, bin_dm.n_ions, size=min(100, bin_dm.n_ions), endpoint=False)

        max_diff = 0.0
        for ion_idx in sample_ions:
            img_bin = bin_dm.get_ion_image(int(ion_idx))
            img_zarr = zarr_dm.get_ion_image(int(ion_idx))
            diff = float(np.max(np.abs(img_bin.astype(np.float64) - img_zarr.astype(np.float64))))
            max_diff = max(max_diff, diff)

        print(f"\n  Zarr vs Binary 最大像素差异: {max_diff:.2e}")
        assert max_diff < 1e-4, f"差异过大: {max_diff}"

        for ion_idx in sample_ions[:10]:
            np.testing.assert_allclose(
                bin_dm.get_ion_vector(int(ion_idx)),
                zarr_dm.get_ion_vector(int(ion_idx)),
                rtol=1e-5,
                atol=1e-8,
            )

        bin_dm.close()
        zarr_dm.close()


class TestZarrDataLoader:
    """Zarr PyTorch DataLoader 全量读取 benchmark。"""

    BATCH_SIZE = 64

    @pytest.fixture(scope="class")
    def zarr_dataset(self, zarr_path):
        """创建 Zarr 数据集，class 级复用。"""
        dm = IonImageDataManagerZarr(str(zarr_path), mode="r")
        ds = IonImageZarrDataset(dm)
        yield ds
        dm.close()

    @staticmethod
    def _run_epoch(dataset, num_workers: int) -> int:
        loader = DataLoader(
            dataset.index_dataset(),
            batch_size=TestZarrDataLoader.BATCH_SIZE,
            shuffle=True,
            num_workers=num_workers,
            pin_memory=False,
            collate_fn=dataset.get_batch,
        )
        total = 0
        for batch in loader:
            total += batch.size(0)
        return total

    def test_get_batch_matches_single_item_stack(self, zarr_dataset):
        """batch-aware Zarr 读取与逐张读取结果一致。"""
        indices = [0, 2, 7, -1]
        indices = [idx for idx in indices if -len(zarr_dataset) <= idx < len(zarr_dataset)]
        batched = zarr_dataset.get_batch(indices)
        legacy = torch.stack([zarr_dataset[i] for i in indices], dim=0)
        assert torch.equal(batched, legacy)

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("num_workers", [0])
    def test_dataloader_epoch(self, benchmark, zarr_dataset, num_workers):
        """全量 epoch 读取时间。"""
        total = benchmark.pedantic(
            lambda: self._run_epoch(zarr_dataset, num_workers),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1,
            setup=_purge_cache,
        )
        assert total == len(zarr_dataset)


class TestZarrStorageLayout:
    """不依赖 example 数据的 Zarr 存储布局回归测试。"""

    def test_ion_images_use_blosc_lz4_without_sharding(self, tmp_path):
        output = tmp_path / "tiny.zarr"
        mz_axis = np.array([100.0, 101.0, 102.0], dtype=np.float64)
        lengths = np.array([3, 3, 3, 3], dtype=np.uint32)
        coordinates = np.array([[1, 1], [2, 1], [1, 2], [2, 2]], dtype=np.uint32)
        intensity = np.arange(1, 13, dtype=np.float32)

        with IonImageDataManagerZarr(
            output,
            mode="w",
            temp_dir=tmp_path,
            mz_dtype=np.float64,
            intensity_dtype=np.float32,
            ion_block=2,
            chunk_size=2,
        ) as dm:
            dm.write_flat_batch(intensity, lengths, coordinates, mz_axis)
            dm.finalize()

        metadata = json.loads((output / "ion_images" / "zarr.json").read_text())

        assert metadata["chunk_grid"]["configuration"]["chunk_shape"] == [2, 2, 2]
        codec_names = [codec["name"] for codec in metadata["codecs"]]
        assert codec_names == ["bytes", "blosc"]

        blosc_config = metadata["codecs"][1]["configuration"]
        assert blosc_config["cname"] == "lz4"
        assert blosc_config["clevel"] == 1
        assert blosc_config["shuffle"] == "bitshuffle"


def _store_size_mb(path: Path) -> float:
    return sum(file.stat().st_size for file in path.rglob("*") if file.is_file()) / (1024 * 1024)
