"""用于共享 m/z 特征轴上加合物过滤的 Numba 内核。"""

import numpy as np
from numba import njit, prange

from massflow.tools.funs import lengths_to_offsets


# ---------------------------------------------------------------------------
# 分批累积候选峰对的 Pearson 统计量
# ---------------------------------------------------------------------------


@njit(parallel=True, cache=True)
def update_pair_statistics_flat_parallel_jit(
    intensity_flat: np.ndarray,
    lengths: np.ndarray,
    left_peak_indices: np.ndarray,
    right_peak_indices: np.ndarray,
    counts: np.ndarray,
    mean_left: np.ndarray,
    mean_right: np.ndarray,
    m2_left: np.ndarray,
    m2_right: np.ndarray,
    comoment: np.ndarray,
) -> None:
    """用一个批次原地更新候选峰对的 Pearson 统计量。

    每个峰对在批内使用 Welford 算法计算均值、离差平方和与
    co-moment，再使用 Chan 合并公式更新跨批次状态。若某个峰索引超出
    当前谱长，该谱上的强度按 0.0 处理。

    ``counts``、``mean_left``、``mean_right``、``m2_left``、
    ``m2_right`` 和 ``comoment`` 的每个元素对应一个峰对。
    最终可以通过 ``comoment / sqrt(m2_left * m2_right)`` 计算
    Pearson 相关系数。
    """
    offsets = lengths_to_offsets(lengths)

    for pair_index in prange(left_peak_indices.size):  # pylint: disable=not-an-iterable
        left_peak_index = int(left_peak_indices[pair_index])
        right_peak_index = int(right_peak_indices[pair_index])

        batch_count = 0
        batch_mean_left = 0.0
        batch_mean_right = 0.0
        batch_m2_left = 0.0
        batch_m2_right = 0.0
        batch_comoment = 0.0

        for spectrum_index in range(lengths.size):
            left_value = 0.0
            right_value = 0.0
            spectrum_length = lengths[spectrum_index]
            spectrum_offset = offsets[spectrum_index]

            if left_peak_index < spectrum_length:
                left_value = intensity_flat[spectrum_offset + left_peak_index]
            if right_peak_index < spectrum_length:
                right_value = intensity_flat[spectrum_offset + right_peak_index]

            batch_count += 1
            delta_left = left_value - batch_mean_left
            batch_mean_left += delta_left / batch_count
            delta_right = right_value - batch_mean_right
            batch_mean_right += delta_right / batch_count
            batch_m2_left += delta_left * (left_value - batch_mean_left)
            batch_m2_right += delta_right * (right_value - batch_mean_right)
            batch_comoment += delta_left * (right_value - batch_mean_right)

        if batch_count == 0:
            continue

        previous_count = int(counts[pair_index])
        if previous_count == 0:
            counts[pair_index] = batch_count
            mean_left[pair_index] = batch_mean_left
            mean_right[pair_index] = batch_mean_right
            m2_left[pair_index] = batch_m2_left
            m2_right[pair_index] = batch_m2_right
            comoment[pair_index] = batch_comoment
            continue

        total_count = previous_count + batch_count
        delta_left = batch_mean_left - mean_left[pair_index]
        delta_right = batch_mean_right - mean_right[pair_index]
        merge_weight = float(previous_count) * float(batch_count) / float(total_count)

        comoment[pair_index] += batch_comoment + delta_left * delta_right * merge_weight
        m2_left[pair_index] += batch_m2_left + delta_left * delta_left * merge_weight
        m2_right[pair_index] += (
            batch_m2_right + delta_right * delta_right * merge_weight
        )
        mean_left[pair_index] += delta_left * float(batch_count) / float(total_count)
        mean_right[pair_index] += delta_right * float(batch_count) / float(total_count)
        counts[pair_index] = total_count
