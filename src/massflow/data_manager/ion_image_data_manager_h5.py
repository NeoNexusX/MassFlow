"""HDF5 离子图像后端 —— bitshuffle+lz4 压缩，每 chunk 一张离子图像。

Phase 1（write_flat_batch）继承自 IonImageDataManagerBinary 的 temp 文件写入。
Phase 2（finalize）从 temp 文件按 ion block 读取 → 向量化散射 → H5 压缩写入。
"""

from __future__ import annotations

import os
import uuid
from concurrent.futures import Future
from contextlib import nullcontext
from copy import deepcopy
from pathlib import Path
from typing import cast

import h5py
import hdf5plugin
import numpy as np

from massflow.concurrency import get_effective_budget
from massflow.concurrency.runtime import io_executor
from massflow.data_manager.ion_image_data_manager import IonImageDataManager
from massflow.tools.infer_spectrum_type import SpectrumType, set_spectrum_type_metadata
from massflow.tools.ion_major_temp_store import BatchInfo, IonMajorTempStore
from massflow.tools.logger import get_logger


# bitshuffle 预处理 + LZ4 最快压缩级别
_H5_ARGS = hdf5plugin.Bitshuffle(nelems=0, cname="lz4", clevel=1)  # type: ignore
logger = get_logger("massflow.data_manager.ion_image_data_manager_h5")


class IonImageDataManagerH5(IonImageDataManager):
    """使用 HDF5 稠密离子图像存储的 continuous MSI 后端。"""

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
        if mode not in {"r", "w"}:
            raise ValueError("mode must be 'r' or 'w'")
        if filepath is None:
            if mode == "r":
                raise ValueError("filepath is required in read mode")
            filepath = Path("output") / f"{uuid.uuid4().hex}.h5"
        self.filepath = Path(filepath)
        self.mode = mode
        self._width = None if width is None else int(width)
        self._height = None if height is None else int(height)
        self.mz_dtype = np.dtype(mz_dtype).newbyteorder("<")
        self.intensity_dtype = np.dtype(intensity_dtype).newbyteorder("<")
        self.ion_block = max(1, int(ion_block))
        self._temp_workers = max(1, int(temp_workers))
        self._merge_workers = max(1, int(merge_workers))
        self._chunk_size = max(1, int(chunk_size))
        self._h5: h5py.File | None = None
        self._dset: h5py.Dataset | None = None
        self._flat_idx: np.ndarray | None = None
        self._mz_axis: np.ndarray | None = None
        self._coordinates: np.ndarray | None = None
        self._coordinate_chunks: list[np.ndarray] = []
        self._n_ions: int | None = None
        self._n_pixels = 0
        self._finalized = False
        self._temp_store: IonMajorTempStore | None = None

        if compute_mean_spectrum:
            logger.warning(
                "IonImageDataManagerH5 does not support mean spectrum "
                "computation; ignoring option."
            )
        if mode == "r":
            self._open_for_read()
            return
        self._validate_metadata(metadata, shared_mz_axis)
        self._temp_store = IonMajorTempStore(
            temp_dir=temp_dir,
            dtype=self.intensity_dtype,
            workers=self._temp_workers,
        )

    @property
    def width(self) -> int:
        if self._width is None:
            raise ValueError("width has not been initialized.")
        return self._width

    @property
    def height(self) -> int:
        if self._height is None:
            raise ValueError("height has not been initialized.")
        return self._height

    @property
    def n_ions(self) -> int:
        if self._n_ions is None:
            raise ValueError("n_ions has not been initialized.")
        return self._n_ions

    @property
    def n_pixels(self) -> int:
        return self._n_pixels

    @property
    def mz_axis(self) -> np.ndarray:
        if self._mz_axis is None:
            raise ValueError("mz_axis has not been initialized.")
        return self._mz_axis

    @property
    def coordinates(self) -> np.ndarray:
        if self._coordinates is None:
            raise ValueError("coordinates have not been initialized.")
        return self._coordinates

    def write_flat_batch(
        self,
        intensity_flat: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
        mz_axis: np.ndarray,
    ) -> None:
        if self.mode != "w":
            raise ValueError("write_flat_batch is only available in write mode.")
        if self._finalized:
            raise ValueError("Cannot write after finalize().")
        if self._temp_store is None:
            raise RuntimeError("Temporary store is not initialized.")

        lengths = np.asarray(lengths)
        coordinates = np.asarray(coordinates)
        intensity = np.asarray(intensity_flat, dtype=self.intensity_dtype)
        mz = np.asarray(mz_axis, dtype=self.mz_dtype)
        if lengths.ndim != 1 or intensity.ndim != 1 or mz.ndim != 1:
            raise ValueError("intensity_flat, lengths, and mz_axis must be 1-D.")
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError("coordinates must have shape (n_pixels, >=2).")
        if lengths.size == 0:
            return
        if self._n_ions is None:
            self._n_ions = int(mz.size)
            self._mz_axis = mz.copy()
        elif self._mz_axis is None or not np.array_equal(mz, self._mz_axis):
            raise ValueError("All batches must use the same shared mz_axis.")
        if np.any(lengths != self._n_ions):
            raise ValueError("Every H5 ion-image spectrum must match mz_axis.")
        if coordinates.shape[0] != lengths.size:
            raise ValueError("coordinates row count must equal lengths size.")
        if intensity.size != lengths.size * self._n_ions:
            raise ValueError("intensity_flat size does not match the batch shape.")

        if np.any(coordinates[:, :2] < 0):
            raise ValueError("coordinates must be non-negative.")
        xy = coordinates[:, :2].astype("<u4", copy=False)
        self._coordinate_chunks.append(xy.copy())
        self._n_pixels += int(lengths.size)
        self._temp_store.append_batch(
            intensity,
            n_pixels=int(lengths.size),
            n_ions=self._n_ions,
        )

    def _build_flat_index(self) -> np.ndarray:
        """Build a checked pixel index from zero- or one-based coordinates."""
        xs = self.coordinates[:, 0].astype(np.int64, copy=False)
        ys = self.coordinates[:, 1].astype(np.int64, copy=False)
        zero_based = (xs.size > 0 and int(xs.min()) == 0) or (
            ys.size > 0 and int(ys.min()) == 0
        )
        if not zero_based:
            xs, ys = xs - 1, ys - 1
        if self._width is None:
            self._width = int(xs.max()) + 1
        if self._height is None:
            self._height = int(ys.max()) + 1
        valid = (xs >= 0) & (xs < self.width) & (ys >= 0) & (ys < self.height)
        if not np.all(valid):
            raise ValueError(
                "coordinates are outside the H5 image bounds "
                f"(width={self.width}, height={self.height})."
            )
        return ys * self.width + xs

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
        self._coordinates = np.vstack(self._coordinate_chunks)
        self._flat_idx = self._build_flat_index()

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

        assert self._flat_idx is not None
        fidx = self._flat_idx

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
            prefetch = get_effective_budget() != 1
            with io_executor(max_workers=1) if prefetch else nullcontext() as executor:
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
                    if executor is not None and index + 1 < len(starts):
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
        self._dset = cast("h5py.Dataset", self._h5["ion_images"])
        self._mz_axis = cast("h5py.Dataset", self._h5["mz_axis"])[:]
        self._coordinates = np.array(self._h5["coordinates"])
        self._width = int(cast(int, self._h5.attrs["width"]))
        self._height = int(cast(int, self._h5.attrs["height"]))
        self._n_ions = int(cast(int, self._h5.attrs["n_ions"]))
        self._n_pixels = int(cast(int, self._h5.attrs["n_pixels"]))
        # 重建散射索引（供 get_ion_vector 使用）
        self._flat_idx = self._build_flat_index()

    def get_ion_image(self, ion_idx: int) -> np.ndarray:
        """读取单 chunk → 解压 → 返回 (H, W)。"""
        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx {ion_idx} 越界 [0, {self.n_ions})")
        assert self._dset is not None
        return np.asarray(self._dset[ion_idx])

    def get_ion_vector(self, ion_idx: int) -> np.ndarray:
        """从稠密图像提取有效像素强度向量。"""
        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx {ion_idx} 越界 [0, {self.n_ions})")
        assert self._dset is not None and self._flat_idx is not None
        return np.asarray(self._dset[ion_idx]).ravel()[self._flat_idx]

    def get_ion_block(self, start: int, end: int) -> np.ndarray:
        if start < 0 or end < start or end > self.n_ions:
            raise IndexError(f"ion range out of bounds: [{start}, {end})")
        if start == end:
            return np.empty((0, self.n_pixels), dtype=self.intensity_dtype)
        assert self._dset is not None and self._flat_idx is not None
        dense = np.asarray(self._dset[start:end])
        return dense.reshape(end - start, -1)[:, self._flat_idx]

    def get_pixel_spectrum(self, pixel_idx: int) -> np.ndarray:
        """H5 格式不支持高效像素谱读取。"""
        raise NotImplementedError("H5 稠密存储不支持 get_pixel_spectrum，请使用 binary 后端")

    def close(self) -> None:
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None
            self._dset = None
            self._flat_idx = None
        self._mz_axis = None
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

    @classmethod
    def from_flat_output_contract(
        cls,
        data_manager,
        filepath,
        *,
        shared_mz_axis: np.ndarray,
        output_spectrum_type: SpectrumType = "centroid",
        temp_dir=None,
        mz_dtype=None,
        intensity_dtype=None,
        ion_block=256,
        temp_workers=2,
        merge_workers=2,
        chunk_size=1,
        compute_mean_spectrum=False,
    ) -> "IonImageDataManagerH5":
        mz_axis = np.asarray(shared_mz_axis)
        if mz_axis.ndim != 1 or mz_axis.size == 0:
            raise ValueError("shared_mz_axis must be a non-empty 1-D array.")
        output_meta = deepcopy(data_manager.ms.meta)
        output_meta.continuous = True
        output_meta.processed = None
        set_spectrum_type_metadata(output_meta, output_spectrum_type)
        return cls(
            filepath=filepath,
            mode="w",
            width=int(output_meta.max_count_of_pixels_x),
            height=int(output_meta.max_count_of_pixels_y),
            temp_dir=temp_dir,
            mz_dtype=data_manager.mz_dtype if mz_dtype is None else mz_dtype,
            intensity_dtype=(
                data_manager.intensity_dtype
                if intensity_dtype is None
                else intensity_dtype
            ),
            ion_block=ion_block,
            temp_workers=temp_workers,
            merge_workers=merge_workers,
            chunk_size=chunk_size,
            metadata=output_meta,
            shared_mz_axis=mz_axis,
            compute_mean_spectrum=compute_mean_spectrum,
        )

    @staticmethod
    def _validate_metadata(metadata, shared_mz_axis: np.ndarray | None) -> None:
        if metadata is None:
            return
        if getattr(metadata, "continuous", None) is not True:
            raise ValueError("H5 ion-image output requires continuous=True.")
        if getattr(metadata, "centroid_spectrum", None) is False:
            raise ValueError("H5 ion-image output requires centroid data.")
        if shared_mz_axis is None or len(shared_mz_axis) == 0:
            raise ValueError("H5 ion-image output requires a shared m/z axis.")


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
