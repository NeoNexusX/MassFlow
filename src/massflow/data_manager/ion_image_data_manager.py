"""Ion-image data manager interfaces."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Self

import numpy as np


class IonImageDataManager(ABC):
    """离子图像存储后端的统一接口。

    这里不绑定具体文件格式；binary、HDF5、Zarr 等后端都可以实现同一组
    读写方法，方便上层 pipeline 只关心“离子图像”语义。
    """

    @abstractmethod
    def write_flat_batch(
        self,
        intensity_flat: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
        mz_axis: np.ndarray,
    ) -> None:
        """
        写入一个 flat batch。
        输入来自现有 flat_generator,布局是 pixel-major:
        intensity_flat.reshape(n_pixels, n_ions)。
        具体后端可以在内部转换成自己的落盘布局。
        """

    @abstractmethod
    def finalize(self) -> None:
        """完成所有延迟写入，让输出文件进入可读取状态。"""

    @abstractmethod
    def close(self) -> None:
        """释放后端资源。"""

    @abstractmethod
    def get_ion_image(self, ion_idx: int) -> np.ndarray:
        """按 ion index 返回一张稠密二维离子图像。"""

    @abstractmethod
    def get_pixel_spectrum(self, pixel_idx: int) -> np.ndarray:
        """按写入顺序的 pixel index 返回该像素的完整强度谱。"""

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type: type | None, exc_val: BaseException | None, exc_tb: Any) -> None:
        self.close()
