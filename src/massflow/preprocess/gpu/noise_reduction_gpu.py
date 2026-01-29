from typing import Optional, Sequence
import numpy as np

import cupy as cp
from cupyx.scipy.ndimage import convolve1d



def _prepare_kernel(coef: Optional[np.ndarray], window: int) -> np.ndarray:
    """准备卷积核：如果未提供coef，则生成全1的移动平均核，并归一化为float32。"""
    if coef is None:
        if not isinstance(window, int) or window <= 0:
            raise ValueError("window must be a positive integer when coef is None")
        window = window + 1 - (window % 2)
        coef = np.ones(window, dtype=np.float32)
    else:
        if not isinstance(coef, np.ndarray) or coef.ndim != 1 or coef.size == 0:
            raise TypeError("coef must be a non-empty 1D numpy array")

    coef = coef.astype(np.float32)
    coef = coef / np.sum(coef)
    return coef


def _convolve1d_gpu(data: np.ndarray, coef: np.ndarray, axis: int) -> np.ndarray:
    """GPU卷积核心：将数据上传至GPU，执行并行卷积，再将结果拉回CPU。"""
    if cp is None or convolve1d is None:
        raise ImportError(
            "cupy is required for ma_gpu. Please install a compatible cupy package."
        )

    # Ensure input is float32
    # User confirmed input is always float32, so we can skip checks and astype
    x = cp.asarray(data)  # Host -> Device (内存 -> 显存)
    w = cp.asarray(coef)

    # 在GPU上并行执行一维卷积
    y = convolve1d(x, w, axis=axis, mode="nearest")
    return cp.asnumpy(y)  # Device -> Host (显存 -> 内存)


def smooth_signal_ma_gpu(
    intensity: np.ndarray,
    coef: Optional[np.ndarray] = None,
    window: int = 5,
) -> np.ndarray:
    """单谱图GPU平滑入口。"""
    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1 or intensity.size == 0:
        raise ValueError("intensity must be a non-empty 1D numpy array")

    coef = _prepare_kernel(coef, window)
    return _convolve1d_gpu(intensity, coef, axis=0)


def smooth_signal_ma_gpu_batch(
    intensities: Sequence[np.ndarray],
    coef: Optional[np.ndarray] = None,
    window: int = 5,
) -> list[np.ndarray]:
    """批量GPU平滑入口：通过矩阵堆叠实现并行加速。"""
    if not intensities:
        return []

    lengths: list[int] = []
    # for arr in intensities:
    #     if not isinstance(arr, np.ndarray) or arr.ndim != 1 or arr.size == 0:
    #         raise ValueError("each intensity must be a non-empty 1D numpy array")
    #     lengths.append(arr.size)

    coef = _prepare_kernel(coef, window)

    n = len(intensities)
    max_len = max(lengths)

    # 创建全零大矩阵进行对齐
    padded = np.zeros((n, max_len), dtype=np.float32)
    for i, arr in enumerate(intensities):
        # Use simple slicing instead of padding for better performance if possible
        # but here we need to align lengths.
        # Direct copy is faster than np.pad for simple edge cases
        valid_len = lengths[i]
        padded[i, :valid_len] = arr.astype(np.float32, copy=False)
        # # Handle edge padding simulation: replicate last value
        # if valid_len < max_len:
        #     padded[i, valid_len:] = padded[i, valid_len - 1]

    # 沿 axis=1 (谱图长度方向) 进行批量卷积
    smoothed = _convolve1d_gpu(padded, coef, axis=1)

    # 结果裁剪：仅保留每条谱图的有效长度
    result: list[np.ndarray] = []
    for i, length in enumerate(lengths):
        # Result is already float32 from GPU
        result.append(smoothed[i, :length])

    return result