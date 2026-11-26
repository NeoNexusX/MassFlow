"""H5 离子图像转写 + PyTorch 读取性能测试。

先通过 fixture 把 example.imzML 的全部 flat 数据缓存到内存，
然后在此基础上测量 H5 转写和读取的纯净时间（排除 imzML 解析开销）。
"""

import subprocess
import time
from pathlib import Path

import numpy as np
import pytest
import torch
from torch.utils.data import DataLoader

from massflow.data_manager import (
    IonImageDataManagerBinary,
    IonImageDataManagerH5,
    MSDataManagerImzML,
)
from massflow.morch import IonImageH5Dataset

EXAMPLE_IMZML = Path("data/example.imzML")
EXAMPLE_IIB = Path("data/example.iib")
EXAMPLE_H5 = Path("data/example.h5")
_SUDO_PASSWORD = "660789"
ROUND = 1


def test_empty_ion_block_keeps_pixel_width(tmp_path):
    """空 ion 区间应返回 (0, n_pixels)，与 numpy 切片语义一致。"""
    with IonImageDataManagerH5(
        tmp_path / "tiny.h5", mode="w", width=2, height=2, temp_dir=tmp_path
    ) as dm:
        dm.write_flat_batch(
            intensity_flat=np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.float32),
            lengths=np.array([2, 2, 2, 2]),
            coordinates=np.array([[1, 1], [2, 1], [1, 2], [2, 2]]),
            mz_axis=np.array([100.0, 200.0]),
        )
        dm.finalize()
        assert dm.get_ion_block(1, 1).shape == (0, dm.n_pixels)


def _purge_cache() -> None:
    """清除系统磁盘缓存，随后暂停 5s 让系统稳定。"""
    try:
        subprocess.run(
            ["sudo", "-S", "purge"],
            input=_SUDO_PASSWORD, text=True, check=True,
            capture_output=True, timeout=30,
        )
        time.sleep(5)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass


def _cache_example_flat_batches(*, batch_size: int = 256):
    """从 example.imzML 读取全部像素，返回内存中的 flat batches 缓存。

    返回 dict: {"batches": [(mz, intensity, lengths, coords), ...],
                "width", "height", "mz_dtype", "intensity_dtype"}
    """
    if not EXAMPLE_IMZML.exists():
        pytest.skip(f"Example data not found at {EXAMPLE_IMZML}")

    with MSDataManagerImzML(filepath=str(EXAMPLE_IMZML)) as dm:
        dm.load_head_data()
        if len(dm.ms) == 0:
            pytest.skip("No spectra in example.imzML")

        batches = []
        for mz_axis, intensity_flat, lengths, coordinates in dm.flat_generator(
            batch_size=batch_size, include_mz=True,
        ):
            batches.append((
                np.asarray(mz_axis).copy(),
                np.asarray(intensity_flat).copy(),
                np.asarray(lengths).copy(),
                np.asarray(coordinates).copy(),
            ))

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
def h5_path(request):
    """class 级 fixture：生成 example.h5，后续测试复用。"""
    if EXAMPLE_H5.exists():
        return EXAMPLE_H5

    flat_cache = request.getfixturevalue("flat_cache")
    h5_dm = IonImageDataManagerH5(
        str(EXAMPLE_H5), mode="w",
        width=flat_cache["width"], height=flat_cache["height"],
        mz_dtype=flat_cache["mz_dtype"],
        intensity_dtype=flat_cache["intensity_dtype"],
        temp_dir="./temp", temp_workers=4, merge_workers=4, ion_block=512,
        chunk_size=1,
    )
    try:
        for mz_axis, intensity_flat, lengths, coordinates in flat_cache["batches"]:
            h5_dm.write_flat_batch(
                intensity_flat=intensity_flat, lengths=lengths,
                coordinates=coordinates, mz_axis=mz_axis,
            )
        h5_dm.finalize()
    finally:
        h5_dm.close()
    return EXAMPLE_H5

def _write_h5_from_cache(flat_cache: dict, output_path: Path, tmp_path: Path) -> None:
    """从内存缓存写出 H5 文件（不含 imzML 读取开销）。"""
    h5_dm = IonImageDataManagerH5(
        str(output_path), mode="w",
        width=flat_cache["width"], height=flat_cache["height"],
        mz_dtype=flat_cache["mz_dtype"],
        intensity_dtype=flat_cache["intensity_dtype"],
        temp_dir=str(tmp_path), temp_workers=4, merge_workers=4, ion_block=512,
        chunk_size=16,
    )
    try:
        for mz_axis, intensity_flat, lengths, coordinates in flat_cache["batches"]:
            h5_dm.write_flat_batch(
                intensity_flat=intensity_flat, lengths=lengths,
                coordinates=coordinates, mz_axis=mz_axis,
            )
        h5_dm.finalize()
    finally:
        h5_dm.close()


class TestH5Convert:
    """H5 转写纯净计时（数据已在内存中）。"""

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_h5_write_benchmark(self, benchmark, flat_cache, tmp_path):
        """benchmark：内存 → H5 完整转写时间（Phase 1 + Phase 2）。"""
        output = tmp_path / "bench.h5"

        def convert():
            _write_h5_from_cache(flat_cache, output, tmp_path)

        benchmark.pedantic(
            convert, rounds=ROUND, iterations=1, warmup_rounds=0, setup=_purge_cache,
        )
        assert output.exists()
        output.unlink(missing_ok=True)

    def test_h5_print_phase_timing(self, flat_cache, tmp_path):
        """打印 Phase 1 / Phase 2 分别耗时（不含 benchmark 框架）。"""
        output = tmp_path / "phase_timing.h5"

        h5_dm = IonImageDataManagerH5(
            str(output), mode="w",
            width=flat_cache["width"], height=flat_cache["height"],
            mz_dtype=flat_cache["mz_dtype"],
            intensity_dtype=flat_cache["intensity_dtype"],
            temp_dir=str(tmp_path), temp_workers=4, merge_workers=4, ion_block=512,
            chunk_size=16,
        )

        t0 = time.perf_counter()
        for mz_axis, intensity_flat, lengths, coordinates in flat_cache["batches"]:
            h5_dm.write_flat_batch(
                intensity_flat=intensity_flat, lengths=lengths,
                coordinates=coordinates, mz_axis=mz_axis,
            )
        t1 = time.perf_counter()
        h5_dm.finalize()
        t2 = time.perf_counter()

        print(f"\n  像素数:     {h5_dm.n_pixels}")
        print(f"  离子数:     {h5_dm.n_ions}")
        print(f"  图像尺寸:   {h5_dm.width} x {h5_dm.height}")
        print(f"  Phase 1 (temp 写入): {t1 - t0:.3f}s")
        print(f"  Phase 2 (H5 合并):   {t2 - t1:.3f}s")
        mb = output.stat().st_size / (1024 * 1024)
        iib_mb = EXAMPLE_IIB.stat().st_size / (1024 * 1024) if EXAMPLE_IIB.exists() else 0
        print(f"  H5 大小: {mb:.0f} MB" + (f"  .iib: {iib_mb:.0f} MB (比 {iib_mb/mb:.1f}x)" if iib_mb else ""))
        print(f"  总转写:     {t2 - t0:.3f}s\n")
        h5_dm.close()
        output.unlink(missing_ok=True)


class TestH5Read:
    """H5 读取性能测试（基于已生成的 example.h5 fixture）。"""

    def test_h5_read_performance(self, h5_path):
        """随机 + 顺序各读取 1000 张离子图像。"""
        dm = IonImageDataManagerH5(str(h5_path), mode="r")
        n_ions = dm.n_ions
        rng = np.random.default_rng(20260430)
        ion_indices = rng.integers(0, n_ions, size=1000, endpoint=False)

        _purge_cache()
        t0 = time.perf_counter()
        checksum = 0.0
        for idx in ion_indices:
            checksum += float(np.sum(dm.get_ion_image(int(idx))))
        t1 = time.perf_counter()
        print(f"\n  H5 随机 1000 张: {t1 - t0:.3f}s ({(t1 - t0) / 1000 * 1000:.2f}ms/张)")
        assert checksum > 0

        _purge_cache()
        t0 = time.perf_counter()
        checksum = 0.0
        for i in range(1000):
            checksum += float(np.sum(dm.get_ion_image(i)))
        t1 = time.perf_counter()
        print(f"  H5 顺序 1000 张: {t1 - t0:.3f}s ({(t1 - t0) / 1000 * 1000:.2f}ms/张)")

        dm.close()

    def test_correctness_vs_binary(self, h5_path):
        """逐 ion 对比 H5 与 binary 文件的图像数据一致性。"""
        if not EXAMPLE_IIB.exists():
            pytest.skip(f"Binary reference not found at {EXAMPLE_IIB}")

        bin_dm = IonImageDataManagerBinary(str(EXAMPLE_IIB), mode="r")
        h5_dm = IonImageDataManagerH5(str(h5_path), mode="r")

        assert bin_dm.n_ions == h5_dm.n_ions
        assert bin_dm.width == h5_dm.width and bin_dm.height == h5_dm.height
        np.testing.assert_allclose(bin_dm.mz_axis, h5_dm.mz_axis, rtol=1e-7, atol=1e-10)

        rng = np.random.default_rng(20260430)
        sample_ions = rng.integers(0, bin_dm.n_ions, size=min(100, bin_dm.n_ions), endpoint=False)

        max_diff = 0.0
        for ion_idx in sample_ions:
            img_bin = bin_dm.get_ion_image(int(ion_idx))
            img_h5 = h5_dm.get_ion_image(int(ion_idx))
            diff = float(np.max(np.abs(
                img_bin.astype(np.float64) - img_h5.astype(np.float64)
            )))
            max_diff = max(max_diff, diff)

        print(f"\n  H5 vs Binary 最大像素差异: {max_diff:.2e}")
        assert max_diff < 1e-4, f"差异过大: {max_diff}"

        for ion_idx in sample_ions[:10]:
            np.testing.assert_allclose(
                bin_dm.get_ion_vector(int(ion_idx)),
                h5_dm.get_ion_vector(int(ion_idx)),
                rtol=1e-5, atol=1e-8,
            )

        bin_dm.close()
        h5_dm.close()


class TestH5DataLoader:
    """H5 PyTorch DataLoader 全量读取 benchmark。"""

    BATCH_SIZE = 64

    @pytest.fixture(scope="class")
    def h5_dataset(self, h5_path):
        """创建 H5 数据集，class 级复用。"""
        dm = IonImageDataManagerH5(str(h5_path), mode="r")
        ds = IonImageH5Dataset(dm)
        yield ds
        dm.close()

    @staticmethod
    def _run_epoch(dataset, num_workers: int) -> int:
        loader = DataLoader(
            dataset.index_dataset(),
            batch_size=TestH5DataLoader.BATCH_SIZE,
            shuffle=True,
            num_workers=num_workers, 
            pin_memory=False,
            collate_fn=dataset.get_batch,
        )
        total = 0
        for batch in loader:
            total += batch.size(0)
        return total

    def test_get_batch_matches_single_item_stack(self, h5_dataset):
        """batch-aware H5 读取与逐张读取结果一致。"""
        indices = [0, 7, 100, -1]
        batched = h5_dataset.get_batch(indices)
        legacy = torch.stack([h5_dataset[i] for i in indices], dim=0)
        assert torch.equal(batched, legacy)

    @pytest.mark.benchmark(timer=time.perf_counter)
    @pytest.mark.parametrize("num_workers", [0])
    def test_dataloader_epoch(self, benchmark, h5_dataset, num_workers):
        """全量 epoch 读取时间。"""
        total = benchmark.pedantic(
            lambda: self._run_epoch(h5_dataset, num_workers),
            rounds=ROUND, iterations=1, warmup_rounds=1, setup=_purge_cache,
        )
        assert total == len(h5_dataset)
