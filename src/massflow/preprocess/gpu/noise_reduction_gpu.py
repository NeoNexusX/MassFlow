from typing import Optional, Sequence
import numpy as np
import cupy as cp
from cupyx.scipy.ndimage import convolve1d
from cupy.cuda.compiler import CompileException

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

    try:
        # 在GPU上并行执行一维卷积
        y = convolve1d(x, w, axis=axis, mode="nearest")
    except CompileException as e:
        raise RuntimeError(
            "CuPy compilation failed. This typically indicates a mismatch between "
            "the installed CuPy version and the CUDA Toolkit version, or an invalid "
            "GPU architecture configuration. Please verify your environment.\n"
            f"Original error: {e}"
        ) from e
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

    lengths: list[int] = [arr.size for arr in intensities]

    coef = _prepare_kernel(coef, window)

    n = len(intensities)
    max_len = max(lengths)

    # 创建全零大矩阵进行对齐
    padded = np.zeros((n, max_len), dtype=np.float32)
    for i, arr in enumerate(intensities):
        # 优化：直接赋值比 np.pad 更快
        valid_len = lengths[i]
        padded[i, :valid_len] = arr.astype(np.float32, copy=False)
        #如果你想测试 np.pad 的性能，可用下面的替代实现（取消注释即可）：
        # padded[i] = np.pad(
        #     arr.astype(np.float32, copy=False),
        #     (0, max_len - valid_len),
        #     mode="edge",
        # )
        
        # 边缘处理：如果有填充区域，用最后一个有效值填充，减少卷积边缘伪影
        if valid_len < max_len:
             padded[i, valid_len:] = padded[i, valid_len - 1]

    # 沿 axis=1 (谱图长度方向) 进行批量卷积
    smoothed = _convolve1d_gpu(padded, coef, axis=1)

    # 结果裁剪：仅保留每条谱图的有效长度
    result: list[np.ndarray] = []
    for i, length in enumerate(lengths):
        # Result is already float32 from GPU
        result.append(smoothed[i, :length])

    return result