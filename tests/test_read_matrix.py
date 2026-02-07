import time
import os
import pytest
import numpy as np
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger

logger = get_logger("test_matrix_read_swap")

ROUND = 1

def data_read2matrix(dm, batch_size, include_mz, max_threads, swap):
    """
    Test: Disk data -> Read to matrix, with/without swap.
    """
    swap_dm = MSDataManagerImzML(MassSpectrumSet(), temp_dir="./tmp")
    swap_dm.copy_meta(dm)

    for mz_data, intensity_matrix, batch_lengths, batch_coordinates in dm.get_matrix_generator(
        batch_size=batch_size,
        include_mz=include_mz,
        max_threads=max_threads,
    ):

        if swap:
            swap_dm.swap_matrix_data_out2disk(mz_data, intensity_matrix, batch_lengths, batch_coordinates)
        else:
            _ = intensity_matrix.sum()

    if swap:
        swap_dm.close_writer()
        return swap_dm

def data_read2batch(dm, batch_size, max_threads, swap):
    """
    Test: Disk data -> Read to batch, with/without swap.
    """
    swap_dm = MSDataManagerImzML(MassSpectrumSet(), temp_dir="./tmp")
    swap_dm.copy_meta(dm)
    writer = swap_dm.writer

    for batch in dm.get_batch_generator(
        batch_size=batch_size,
        max_threads=max_threads,
    ):

        if swap:
            for spec in batch:
                writer.addSpectrum(spec.mz_list, spec.intensity, spec.coordinate.get_tuple())
        else:
            _ = batch

        dm.clear_batch_data_memory(batch)

    if swap:
        swap_dm.close_writer()
        return swap_dm

def _collect_batch_data(dm):
    """
    collect all data via get_batch_generator.
    """
    results = []
    for batch in dm.get_batch_generator(batch_size=1024, max_threads=8):
        for spec in batch:
            results.append({
                "mz": np.array(spec.mz_list, copy=True),
                "intensity": np.array(spec.intensity, copy=True),
                "coord": spec.coordinate.get_tuple(),
            })
        dm.clear_batch_data_memory(batch)
    return results

def _collect_matrix_data(dm):
    """
    collect all data via get_matrix_generator (include_mz=True).
    """
    results = []
    for mz_data, intensity_matrix, lengths, coordinates in dm.get_matrix_generator(
        batch_size=1024, include_mz=True, max_threads=8,
    ):
        n = intensity_matrix.shape[0]
        is_continuous = (mz_data is not None and mz_data.ndim == 1)
        for i in range(n):
            length = int(lengths[i])
            intensity = intensity_matrix[i, :length].copy()
            mz = mz_data[:length].copy() if is_continuous else mz_data[i, :length].copy()
            coord = tuple(int(c) for c in coordinates[i])

            results.append({
                "mz": mz,
                "intensity": intensity,
                "coord": coord,
            })
    return results

def _compare_spectrum_lists(list_a, list_b):
    """
    严格比较两组谱图数据列表的一致性。

    检查项：
        1. 谱图总数一致
        2. 每条谱图的坐标一致
        3. 每条谱图的 m/z 数组长度一致
        4. 每条谱图的 intensity 数组长度一致
        5. m/z 值逐元素一致 (浮点容差 atol=1e-6, rtol=1e-6)
        6. intensity 值逐元素一致 (浮点容差 atol=1e-4, rtol=1e-5)

    Args:
        list_a: 第一组数据
        list_b: 第二组数据

    Raises:
        AssertionError: 任意检查项失败时
    """
    # 总数一致
    assert len(list_a) == len(list_b), (f"谱图数量不一致: batch_read={len(list_a)}, matrix_read={len(list_b)}")

    # 对按顺序排列的谱图逐一比较
    for idx, (a, b) in enumerate(zip(list_a, list_b)):

        # 判断坐标一致性；m/z、intensity长度一致
        assert a["coord"] == b["coord"], (f"spectrum #{idx} (batch_read vs matrix_read): 坐标不一致  {a['coord']} != {b['coord']}")
        assert len(a["mz"]) == len(b["mz"]), (f"spectrum #{idx} (batch_read vs matrix_read): m/z 长度不一致  {len(a['mz'])} != {len(b['mz'])}")
        assert len(a["intensity"]) == len(b["intensity"]), (
            f"spectrum #{idx} (batch_read vs matrix_read): intensity 长度不一致  {len(a['intensity'])} != {len(b['intensity'])}"
        )

        # m/z 值比较
        np.testing.assert_allclose(
            a["mz"], b["mz"],
            atol=1e-6, rtol=1e-6,
            err_msg=f"spectrum #{idx} (batch_read vs matrix_read): m/z 值不一致"
        )

        # intensity 值比较
        np.testing.assert_allclose(
            a["intensity"], b["intensity"],
            atol=1e-4, rtol=1e-5,
            err_msg=f"spectrum #{idx} (batch_read vs matrix_read): intensity 值不一致"
        )

    logger.info(f"数据一致性验证通过 (batch_read vs matrix_read): 共 {len(list_a)} 条谱图全部匹配。")

def _compare_by_coordinate_map(list_a, list_b):
    """
    基于坐标映射比较两组谱图数据（顺序无关）。

    检查项：
        1. 谱图总数一致
        2. 坐标集合完全一致（无遗漏、无多余）
        3. 相同坐标下的 m/z 值一致
        4. 相同坐标下的 intensity 值一致

    Raises:
        AssertionError: 任意检查项失败时
    """
    assert len(list_a) == len(list_b), (
        f"谱图数量不一致: original={len(list_a)}, matrix_swap_readback={len(list_b)}"
    )

    # 构建坐标 -> 数据映射
    map_a = {item["coord"]: item for item in list_a}
    map_b = {item["coord"]: item for item in list_b}

    # 检测坐标集合是否一致
    coords_a = set(map_a.keys())
    coords_b = set(map_b.keys())
    only_in_a = coords_a - coords_b
    only_in_b = coords_b - coords_a
    assert coords_a == coords_b, (
        f"坐标集合不一致: "
        f"仅在 original 中: {list(only_in_a)[:5]}{'...' if len(only_in_a) > 5 else ''}, "
        f"仅在 matrix_swap_readback 中: {list(only_in_b)[:5]}{'...' if len(only_in_b) > 5 else ''}"
    )

    # 逐坐标比较
    for coord in sorted(coords_a):
        a = map_a[coord]
        b = map_b[coord]

        # 长度
        assert len(a["mz"]) == len(b["mz"]), (
            f"coord {coord} (original vs matrix_swap_readback): m/z 长度不一致  {len(a['mz'])} != {len(b['mz'])}"
        )
        assert len(a["intensity"]) == len(b["intensity"]), (
            f"coord {coord} (original vs matrix_swap_readback): intensity 长度不一致  {len(a['intensity'])} != {len(b['intensity'])}"
        )

        # m/z、intensity 值
        np.testing.assert_allclose(
            a["mz"], b["mz"],
            atol=1e-6, rtol=1e-6,
            err_msg=f"coord {coord} (original vs matrix_swap_readback): m/z 值不一致"
        )
        np.testing.assert_allclose(
            a["intensity"], b["intensity"],
            atol=1e-4, rtol=1e-5,
            err_msg=f"coord {coord} (original vs matrix_swap_readback): intensity 值不一致"
        )

    logger.info(
        f"基于坐标映射的数据一致性验证通过 (original vs matrix_swap_readback): "
        f"共 {len(list_a)} 条谱图，{len(coords_a)} 个唯一坐标全部匹配。"
    )

def verify_read_data_consistency(dm):
    """
    验证 batch 读取模式与 matrix 读取模式的数据一致性。

    核心逻辑：
        1. 通过 get_batch_generator 逐谱图收集所有 mz, intensity, coordinate
        2. 通过 get_matrix_generator (include_mz=True) 逐行收集同样的数据
        3. 严格比较两者在谱图数量、顺序、坐标、m/z、intensity 上的完全一致性

    这是数据正确性的基础验证：确保两种读取通道对同一原始文件产出完全一致的结果。
    """
    # 收集 batch 数据
    batch_data = _collect_batch_data(dm)
    matrix_data = _collect_matrix_data(dm)

    _compare_spectrum_lists(batch_data, matrix_data)
    logger.info("READ 一致性验证完成: batch、matrix 读取数据完全一致。")

def verify_swap_data_consistency(dm):
    """
    验证 swap（写入磁盘再读回）后的数据正确性。

    核心逻辑：
        1. 从原始数据通过 batch 路径收集基准数据 (ground truth)
        2. 通过 matrix 路径 swap 出到磁盘 → 重新读回 → 通过 batch 路径收集数据
        3. 比较: 原始 vs matrix-swap 读回 → 验证 matrix swap 写入/读回的数据保真性

    使用坐标映射比较，因为 swap 后重新读取时谱图顺序可能不同。
    """
    logger.info("开始验证 SWAP 数据一致性 (原始 vs swap 读回)")

    # 收集原始基准数据
    original_data = _collect_batch_data(dm)

    # matrix swap → 读回
    matrix_swap_dm = data_read2matrix(dm, batch_size=1024, include_mz=True, max_threads=8, swap=True)
    matrix_swap_dm.load_full_data_from_file()
    matrix_swap_readback = _collect_batch_data(matrix_swap_dm)

    # 原始 vs matrix-swap 读回
    _compare_by_coordinate_map(original_data, matrix_swap_readback)

    # 清理临时资源
    matrix_swap_dm.close()

    logger.info("SWAP 一致性验证完成: matrix-swap 读回与原始数据完全一致，测试通过。")

class TestMatrixGenerator:
    """
    测试 MSDataManagerImzML 的 get_matrix_generator 方法.
    """

    @pytest.fixture(scope="class")
    def ms_raw_data(self):

        filepath = r"E:\data\data\example\example.imzML"
        if not os.path.exists(filepath):
            pytest.skip(f"Original data not found at {filepath}")
        # loda meta data
        dm_raw = MSDataManagerImzML(filepath=filepath)
        dm_raw.load_full_data_from_file()

        return dm_raw

    @pytest.mark.parametrize("batch_size", [1024])
    @pytest.mark.parametrize("max_threads", [8])
    @pytest.mark.parametrize(("include_mz", "swap"), [(True, False), (True, True)]) # swap=True need include_mz=True
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_matrix_read_swap_performance(self, benchmark, ms_raw_data, batch_size, max_threads, include_mz, swap):
        """
        Performance test: Data read to matrix with/without swap.
        """

        dm_result = benchmark.pedantic(
            data_read2matrix,
            args=(ms_raw_data, batch_size, include_mz, max_threads, swap),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1
        )

        if swap:
            dm_result.close()

    @pytest.mark.parametrize("batch_size", [1024])
    @pytest.mark.parametrize("max_threads", [8])
    @pytest.mark.parametrize("swap", [True, False])
    @pytest.mark.benchmark(timer=time.perf_counter)
    def test_batch_read_swap_performance(self, benchmark, ms_raw_data, batch_size, max_threads, swap):
        """
        Performance test: Data read to batch with/without swap.
        """

        dm_result = benchmark.pedantic(
            data_read2batch,
            args=(ms_raw_data, batch_size, max_threads, swap),
            rounds=ROUND,
            iterations=1,
            warmup_rounds=1
        )

        if swap:
            dm_result.close()

    def test_verify_read_data_consistency(self, ms_raw_data):
        """
        Verify data consistency between batch and matrix read (without swap).
        """

        verify_read_data_consistency(ms_raw_data)

    def test_verify_swap_data_consistency(self, ms_raw_data):
        """
        Verify data correctness after swap (write to disk and read back).
        """

        verify_swap_data_consistency(ms_raw_data)
