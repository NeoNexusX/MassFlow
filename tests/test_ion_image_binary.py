import subprocess
import time
from pathlib import Path

import numpy as np
import pytest

from massflow.data_manager import IonImageDataManagerBinary, MSDataManagerImzML
from massflow.preprocess import Preprocessor

EXAMPLE_IMZML = Path("data/example.imzML")
ROUND = 10

_SUDO_PASSWORD = "660789"


def _purge_cache() -> None:
    """通过 sudo purge 清除磁盘缓存，然后暂停 10s 让系统稳定。"""
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

def _cache_example_flat_batches(
    *,
    target_locs=None,
    batch_size: int = 512,
):
    """从真实 example.imzML 读取全部像素，并缓存 flat batches。

    测试刻意使用项目自带 example 文件，而不是人工构造数据；默认不设置
    target_locs，因此会覆盖 example.imzML 中的全部像素。
    """
    if not EXAMPLE_IMZML.exists():
        pytest.skip(f"Example data not found at {EXAMPLE_IMZML}")

    with MSDataManagerImzML(filepath=str(EXAMPLE_IMZML), target_locs=target_locs) as dm:
        dm.load_head_data()
        if len(dm.ms) == 0:
            pytest.skip(f"No spectra found in target_locs={target_locs}")

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
def full_example_flat_cache():
    return _cache_example_flat_batches()


@pytest.fixture(scope="class")
def full_example_binary_path(tmp_path_factory, full_example_flat_cache):
    """预生成一份完整 example binary，供正确性测试和随机读取 benchmark 共用。

    这一步不是 benchmark；它的目的只是准备可读的完整 binary 文件，避免
    roundtrip 测试和随机读取测试各自重复写出大文件。
    """
    output_dir = tmp_path_factory.mktemp("ion_image_binary_full")
    output_path = output_dir / "example_full.iib"
    return _write_cached_binary(full_example_flat_cache, output_path, output_dir)


def _write_cached_binary(cache, output_path: Path, temp_dir: Path) -> Path:
    """写出 binary。

    这个函数本身不做手动计时；耗时统一交给 pytest-benchmark 插件统计。
    """
    output_path.unlink(missing_ok=True)
    manager = IonImageDataManagerBinary(
        output_path,
        mode="w",
        width=cache["width"],
        height=cache["height"],
        temp_dir=temp_dir,
        mz_dtype=cache["mz_dtype"],
        intensity_dtype=cache["intensity_dtype"],
        ion_block=16,
        temp_workers=1,
        merge_workers=1,
    )
    try:
        for mz_axis, intensity_flat, lengths, coordinates in cache["batches"]:
            manager.write_flat_batch(
                intensity_flat=intensity_flat,
                lengths=lengths,
                coordinates=coordinates,
                mz_axis=mz_axis,
            )
        manager.finalize()
    finally:
        manager.close()

    return output_path


def _expected_pixel_spectrum(cache, pixel_idx: int) -> np.ndarray:
    """从原始 flat cache 中还原某个像素的期望强度谱。"""
    offset = 0
    for _, intensity_flat, lengths, _coordinates in cache["batches"]:
        n_pixels = len(lengths)
        if offset <= pixel_idx < offset + n_pixels:
            local_idx = pixel_idx - offset
            n_ions = int(lengths[local_idx])
            start = local_idx * n_ions
            return intensity_flat[start : start + n_ions]
        offset += n_pixels
    raise IndexError(pixel_idx)


def _expected_ion_vector(cache, ion_idx: int) -> np.ndarray:
    """从原始 flat cache 中还原某个 ion 的期望 pixel 向量。"""
    values = []
    for _, intensity_flat, lengths, _coordinates in cache["batches"]:
        n_ions = int(lengths[0])
        values.append(intensity_flat.reshape(len(lengths), n_ions)[:, ion_idx])
    return np.concatenate(values)


class TestIonImageBinary:
    def test_binary_roundtrip_with_example_flat_batches(self, full_example_flat_cache, full_example_binary_path):
        """校验完整 example 数据写成 binary 后仍能按 pixel/ion 正确读回。"""
        assert full_example_binary_path.exists()

        with IonImageDataManagerBinary(full_example_binary_path, mode="r") as binary_dm:
            first_mz = full_example_flat_cache["batches"][0][0]
            assert binary_dm.mz_axis is not None
            assert binary_dm.n_ions == len(first_mz)
            assert binary_dm.n_pixels == sum(len(batch[2]) for batch in full_example_flat_cache["batches"])
            np.testing.assert_allclose(binary_dm.mz_axis, first_mz, rtol=1e-7, atol=1e-10)

            np.testing.assert_allclose(
                binary_dm.get_pixel_spectrum(0),
                _expected_pixel_spectrum(full_example_flat_cache, 0),
                rtol=1e-5,
                atol=1e-8,
            )

            ion_idx = min(1000, binary_dm.n_ions - 1)
            np.testing.assert_allclose(
                binary_dm.get_ion_vector(ion_idx),
                _expected_ion_vector(full_example_flat_cache, ion_idx),
                rtol=1e-5,
                atol=1e-8,
            )

            image = binary_dm.get_ion_image(ion_idx)
            assert image.shape == (full_example_flat_cache["height"], full_example_flat_cache["width"])
            assert np.count_nonzero(image) > 0

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_random_read_1000_ions_total_time(self, benchmark, full_example_binary_path):
        """benchmark：随机读取 1000 个离子图像的总时间。

        读取用的完整 binary 来自 fixture；pytest-benchmark 只统计下面
        read_1000_ion_images 的随机读取循环，不包含任何写出准备时间。
        """
        try:
            with IonImageDataManagerBinary(full_example_binary_path, mode="r") as binary_dm:
                rng = np.random.default_rng(20260428)
                ion_indices = rng.integers(0, binary_dm.n_ions, size=1000, endpoint=False)

                def read_1000_ion_images():
                    # pytest-benchmark 会统计整个循环的耗时；checksum 用于防止读取
                    # 结果被测试逻辑完全忽略。
                    checksum = 0.0
                    for ion_idx in ion_indices:
                        checksum += float(np.sum(binary_dm.get_ion_image(int(ion_idx))))
                    return checksum

                checksum = benchmark.pedantic(
                    read_1000_ion_images, rounds=ROUND, iterations=1, warmup_rounds=0, setup=_purge_cache
                )
                assert checksum >= 0.0
        finally:
            # 随机读 benchmark 完成后释放预生成的完整 binary，避免后续写出
            # benchmark 同时持有多份 5GB 级文件。
            full_example_binary_path.unlink(missing_ok=True)

    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_full_example_write_and_merge_total_time(self, benchmark, tmp_path, full_example_flat_cache):
        """benchmark：完整 example 文件写出和合并的总时间。

        pytest-benchmark 会统计 run_write 的耗时；run_write 覆盖了临时 batch
        写出和 finalize 合并两个阶段。
        """
        def run_write():
            print("Starting write benchmark...")
            return _write_cached_binary(full_example_flat_cache, tmp_path / "write_benchmark.iib", tmp_path)

        output_path = benchmark.pedantic(
            run_write, rounds=ROUND, iterations=1, warmup_rounds=0, setup=_purge_cache
        )
        assert output_path.exists()
        output_path.unlink(missing_ok=True)

    def test_pipeline_noise_reduction_to_binary(self, tmp_path):
        """测试使用 pipeline 和 noise reduction，输出为 binary 格式，并检查 head 部分。"""
        if not EXAMPLE_IMZML.exists():
            pytest.skip(f"Example data not found at {EXAMPLE_IMZML}")

        output_path = tmp_path / "test_pipeline.iib"
        with MSDataManagerImzML(filepath=str(EXAMPLE_IMZML)) as dm:
            dm.load_head_data()
            
            processed_dm = (
                Preprocessor(
                    data_manager=dm,
                    output_format="ion_image_binary",
                    output_path=str(output_path),
                )
                .noise_reduction(method="ma_numba", window=5)
                .start()
            )
            
            try:
                assert isinstance(processed_dm, IonImageDataManagerBinary)
                assert processed_dm.n_ions > 0
                assert processed_dm.n_pixels > 0
                assert processed_dm.mz_axis is not None
                assert len(processed_dm.mz_axis) == processed_dm.n_ions
                assert processed_dm.header is not None
                assert processed_dm.header.width > 0
                assert processed_dm.header.height > 0
                if output_path.exists():
                    print(f"Merged binary size: {output_path.stat().st_size} bytes")
            finally:
                processed_dm.close()

