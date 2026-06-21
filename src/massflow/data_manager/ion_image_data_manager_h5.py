"""HDF5 离子图像后端 —— bitshuffle+lz4 压缩，每 chunk 一张离子图像。

Phase 1（write_flat_batch）继承自 IonImageDataManagerBinary 的 temp 文件写入。
Phase 2（finalize）从 temp 文件按 ion block 读取 → 向量化散射 → H5 压缩写入。
"""

from __future__ import annotations

import os
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from typing import cast

import h5py
import hdf5plugin
import numpy as np

from massflow.data_manager.ion_image_data_manager_binary import IonImageDataManagerBinary
from massflow.tools.ion_major_temp_store import BatchInfo, IonMajorTempStore
from massflow.tools.logger import get_logger


# bitshuffle 预处理 + LZ4 最快压缩级别
_H5_ARGS = hdf5plugin.Bitshuffle(nelems=0, cname="lz4", clevel=1)  # type: ignore
logger = get_logger("massflow.data_manager.ion_image_data_manager_h5")


class IonImageDataManagerH5(IonImageDataManagerBinary):
    """HDF5 离子图像后端。

    继承 binary 版本的 Phase 1（temp 文件转置写入），
    重写 finalize 实现 H5 合并写入，重写读路径直接解压 chunk。codex: 注释不应该这样写，不要提重写之类的，别人又没有看过以前的代码
    """

    def __init__(
        self,
        filepath: str | os.PathLike[str] | None = None,
        *,
        mode: str = "r",
        width: int | None = None,
        height: int | None = None,
        temp_dir: str | os.PathLike[str] | None = None,
        mz_dtype: np.dtype | type = np.float64,
        intensity_dtype: np.dtype | type = np.float32,
        ion_block: int = 256,
        temp_workers: int = 4,
        merge_workers: int = 4,
        chunk_size: int = 1,
        metadata=None,
        shared_mz_axis: np.ndarray | None = None,
        compute_mean_spectrum: bool = False,
    ):
        import uuid

        # 必须在 super().__init__ 前初始化 —— 父类读模式会立即调用 _open_for_read()
        self._h5: h5py.File | None = None
        self._dset: h5py.Dataset | None = None
        self._flat_idx: np.ndarray | None = None
        self._chunk_size = max(1, int(chunk_size))

        if filepath is None and mode == "w":
            filepath = Path("output") / f"{uuid.uuid4().hex}.h5"
        super().__init__(
            filepath=filepath, mode=mode, width=width, height=height,
            temp_dir=temp_dir, mz_dtype=mz_dtype, intensity_dtype=intensity_dtype,
            ion_block=ion_block, temp_workers=temp_workers, merge_workers=merge_workers,
            metadata=metadata, shared_mz_axis=shared_mz_axis,
            compute_mean_spectrum=compute_mean_spectrum,
        )

    # ---- properties ---------------------------------------------------------

    @property
    def dset(self) -> "h5py.Dataset":
        if self._dset is None:
            self._open_for_read()
        if self._dset is None:
            raise ValueError("H5 dataset has not been initialized.")
        return self._dset

    @dset.setter
    def dset(self, value: "h5py.Dataset | None"):
        self._dset = value

    @property
    def flat_idx(self) -> np.ndarray:
        if self._flat_idx is None:
            self._open_for_read()
        if self._flat_idx is None:
            raise ValueError("flat_idx has not been initialized.")
        return self._flat_idx

    @flat_idx.setter
    def flat_idx(self, value: np.ndarray | None):
        self._flat_idx = value

    # ---- Phase 2 ------------------------------------------------------------

    def finalize(self) -> None:
        """合并 temp 文件 → H5（bitshuffle+lz4 压缩）。"""
        if self.mode != "w" or self._finalized:
            return
        if (
            self._temp_store is None
            or self._n_ions is None
            or self._mz_axis is None
            or not self._coordinate_chunks
        ):
            raise ValueError("No ion-image batches were written.")
        batch_infos = self._temp_store.finish()
        self.coordinates = np.vstack(self._coordinate_chunks)
        if self.width is None:
            self.width = int(np.max(self.coordinates[:, 0]))
        if self.height is None:
            self.height = int(np.max(self.coordinates[:, 1]))

        self.filepath.parent.mkdir(parents=True, exist_ok=True)
        self._merge_to_h5(batch_infos)
        self._finalized = True
        self._temp_store.close()
        self._temp_store = None
        self._open_for_read()

    def _merge_to_h5(self, batch_infos: list[BatchInfo]) -> None:
        """从 temp 文件读取并散射写入 H5。"""
        n, h, w = self._n_ions, self.height, self.width
        npix = self._n_pixels
        dtype = self.intensity_dtype
        assert n is not None and self._temp_store is not None

        # 散射索引（所有 ion 共用）
        xs = self.coordinates[:, 0].astype(np.int64, copy=False)
        ys = self.coordinates[:, 1].astype(np.int64, copy=False)
        if xs.size > 0 and int(xs.min()) != 0:
            xs, ys = xs - 1, ys - 1
        valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
        fidx = (ys[valid] * w + xs[valid])
        self.flat_idx = fidx

        with h5py.File(self.filepath, "w") as f:
            f.attrs["width"], f.attrs["height"] = w, h
            f.attrs["n_ions"], f.attrs["n_pixels"] = n, npix
            f.create_dataset("mz_axis", data=self.mz_axis)
            f.create_dataset("coordinates", data=self.coordinates, dtype="<u4")
            dset = f.create_dataset(
                "ion_images", shape=(n, h, w),
                chunks=(self._chunk_size, h, w), dtype=dtype, **_H5_ARGS,
            )
            block = max(16, self.ion_block)
            starts = list(range(0, n, block))
            with ThreadPoolExecutor(max_workers=1) as executor:
                future: Future[np.ndarray] | None = None
                for index, start in enumerate(starts):
                    end = min(n, start + block)
                    if future is None:
                        current = _read_scatter(
                            self._temp_store,
                            batch_infos,
                            fidx,
                            start,
                            end,
                            npix,
                            h,
                            w,
                            dtype,
                        )
                    else:
                        current = future.result()
                    future = None
                    if index + 1 < len(starts):
                        next_start = starts[index + 1]
                        next_end = min(n, next_start + block)
                        future = executor.submit(
                            _read_scatter,
                            self._temp_store,
                            batch_infos,
                            fidx,
                            next_start,
                            next_end,
                            npix,
                            h,
                            w,
                            dtype,
                        )
                    dset[start:end] = current

    # ---- 读路径 --------------------------------------------------------------

    def _open_for_read(self) -> None:
        if self._h5 is not None:
            return
        self._h5 = h5py.File(self.filepath, "r")
        self.dset = cast("h5py.Dataset", self._h5["ion_images"])
        self.mz_axis = cast("h5py.Dataset", self._h5["mz_axis"])[:]
        self.coordinates = np.array(self._h5["coordinates"])
        self.width = int(cast(int, self._h5.attrs["width"]))
        self.height = int(cast(int, self._h5.attrs["height"]))
        self._n_ions = int(cast(int, self._h5.attrs["n_ions"]))
        self._n_pixels = int(cast(int, self._h5.attrs["n_pixels"]))
        # 重建散射索引（供 get_ion_vector 使用）
        xs = self.coordinates[:, 0].astype(np.int64, copy=False)
        ys = self.coordinates[:, 1].astype(np.int64, copy=False)
        if xs.size > 0 and int(xs.min()) != 0:
            xs, ys = xs - 1, ys - 1
        w, h = self.width, self.height
        valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
        self._flat_idx = (ys[valid] * w + xs[valid])

    def get_ion_image(self, ion_idx: int) -> np.ndarray:
        """读取单 chunk → 解压 → 返回 (H, W)。"""
        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx {ion_idx} 越界 [0, {self.n_ions})")
        return np.asarray(self.dset[ion_idx])

    def get_ion_vector(self, ion_idx: int) -> np.ndarray:
        """从稠密图像提取有效像素强度向量。"""
        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx {ion_idx} 越界 [0, {self.n_ions})")
        return np.asarray(self.dset[ion_idx]).ravel()[self.flat_idx]

    def get_pixel_spectrum(self, pixel_idx: int) -> np.ndarray:
        """H5 格式不支持高效像素谱读取。"""
        raise NotImplementedError("H5 稠密存储不支持 get_pixel_spectrum，请使用 binary 后端")

    def close(self) -> None:
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None
            self.dset = None
        self.mz_axis = None
        self._coordinates = None
        if self._temp_store is not None:
            self._temp_store.close()
            self._temp_store = None

    @classmethod
    def from_data_manager(cls, data_manager, filepath, *, temp_dir=None,
                          mz_dtype=None, intensity_dtype=None, ion_block=256,
                          temp_workers=2, merge_workers=2, chunk_size=1,
                          compute_mean_spectrum=False,
                          ) -> "IonImageDataManagerH5":
        meta = data_manager.ms.meta

        if (getattr(meta, "continuous", None) is not True or 
            data_manager.ms.shared_mz_list is None or 
            len(data_manager.ms.shared_mz_list) == 0):
            raise ValueError("from_data_manager requires continuous=True metadata and a shared mz_axis.")

        return cls(
            filepath=filepath, mode="w",
            width=int(meta.max_count_of_pixels_x),
            height=int(meta.max_count_of_pixels_y),
            temp_dir=temp_dir,
            mz_dtype=data_manager.mz_dtype if mz_dtype is None else mz_dtype,
            intensity_dtype=data_manager.intensity_dtype if intensity_dtype is None else intensity_dtype,
            ion_block=ion_block, temp_workers=temp_workers, merge_workers=merge_workers,
            chunk_size=chunk_size,
            metadata=meta, shared_mz_axis=data_manager.ms.shared_mz_list,
            compute_mean_spectrum=compute_mean_spectrum,
        )


# ---- 辅助函数（模块级，便于 pickle）-------------------------------------------

def _read_scatter(
    temp_store: IonMajorTempStore,
    batch_infos: list[BatchInfo],
    flat_idx: np.ndarray,
    ion_start: int,
    ion_end: int,
    n_pixels: int,
    height: int,
    width: int,
    dtype: np.dtype,
) -> np.ndarray:
    """从所有 temp 文件读取一个 ion block，向量化散射为 (block, H, W)。"""
    b = ion_end - ion_start
    sparse = np.empty((b, n_pixels), dtype=dtype)
    p = 0
    for bi in batch_infos:
        sparse[:, p:p + bi.n_pixels] = temp_store.read_ion_block(
            bi, ion_start, ion_end
        )
        p += bi.n_pixels
    dense = np.zeros((b, height * width), dtype=dtype)
    dense[:, flat_idx] = sparse
    return dense.reshape(b, height, width) # pylint: disable=E1121
