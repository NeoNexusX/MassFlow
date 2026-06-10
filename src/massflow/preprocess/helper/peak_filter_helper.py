"""Peak filtering helper functions for mass spectrometry data preprocessing."""

from typing import Literal

import numpy as np

from massflow.data_manager import MSDataManagerImzML


# ---------------------------------------------------------------
# Helper functions for peak filtering by intensity
# ---------------------------------------------------------------
def _input_validation(
    data_manager: MSDataManagerImzML,
    method: Literal["mean", "sum"],
    quantile: float,
) -> None:
    """Validate the inputs for peak filtering."""
    if data_manager.ms.meta.continuous is not True:
        raise ValueError("Peak filtering is only applicable to continuous data.")
    if method not in {"mean", "sum"}:
        raise ValueError("Method must be either 'mean' or 'sum'.")
    if not 0 < quantile < 1:
        raise ValueError("Quantile must be between 0 and 1.")


def compute_sum_intensity(data_manager: MSDataManagerImzML) -> np.ndarray:
    """Compute the sum spectrum across all spectra in the data manager."""
    # 计算的时候同时统计一下数量，方便计算平均谱
    # 使用 flat_generator 使用numba加速批量并行计算，具体逻辑（相关到numba）下方到 peak_filter_numba.py 中实现
    pass


def compute_mean_intensity(data_manager: MSDataManagerImzML) -> np.ndarray:
    """Compute the mean spectrum across all spectra in the data manager."""
    # 这里可以调用 compute_sum_spectrum 来计算总谱，然后除以光谱数量来得到平均谱
    pass


def filter_peaks_by_intensity(
    data_manager: MSDataManagerImzML,
    *,
    method: Literal["mean", "sum"] = "mean",
    quantile: float = 0.05,
) -> MSDataManagerImzML:
    """Filter peaks based on intensity."""

    _input_validation(data_manager, method, quantile)

    if method == "sum":
        spectrum = compute_sum_intensity(data_manager)
    else:
        spectrum = compute_mean_intensity(data_manager)

    # 根据计算得到的谱来过滤峰，保留强度在指定分位数以上的峰
    # 之后创建新 的 MSDataManagerImzML 对象，批量写入 ibd 文件中

# 注意： 请参考 massflow.preprocess 下的代码，涉及 numba 并行加速/批量处理逻辑
# 以及创建新

