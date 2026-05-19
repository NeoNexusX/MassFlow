"""
Zarr 离子图像后端。

Phase 1（write_flat_batch）继承自 IonImageDataManagerBinary 的 temp 文件写入。
Phase 2（finalize）从 temp 文件按 ion block 读取，散射成稠密 ion image 后写入 Zarr。
"""

from __future__ import annotations

import os
from concurrent.futures import Future
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import zarr

from massflow.data_manager.ion_image_data_manager_binary import (
    _BatchInfo,
    IonImageDataManagerBinary,
)
from massflow.tools.logger import get_logger

logger = get_logger("massflow.data_manager.ion_image_data_manager_zarr")


class IonImageDataManagerZarr(IonImageDataManagerBinary):
    """
    Zarr 离子图像后端。

    继承 binary 版本的 Phase 1（temp 文件转置写入），重写 finalize 和读路径。
    文件内主数组为 ``ion_images``，shape 为 ``(n_ions, height, width)``。
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
        chunk_size: int = 16,
        metadata=None,
        shared_mz_axis: np.ndarray | None = None,
    ):
        import uuid

        self._group: Any | None = None
        self._dset: Any | None = None
        self._flat_idx: np.ndarray | None = None
        self._source_idx: np.ndarray | None = None
        self._coordinate_base: str | None = None
        self._chunk_size = max(1, int(chunk_size))

        if filepath is None and mode == "w":
            filepath = Path("output") / f"{uuid.uuid4().hex}.zarr"
        super().__init__(
            filepath=filepath,
            mode=mode,
            width=width,
            height=height,
            temp_dir=temp_dir,
            mz_dtype=mz_dtype,
            intensity_dtype=intensity_dtype,
            ion_block=ion_block,
            temp_workers=temp_workers,
            merge_workers=merge_workers,
            metadata=metadata,
            shared_mz_axis=shared_mz_axis,
        )

    @property
    def dset(self):
        if self._dset is None:
            self._open_for_read()
        if self._dset is None:
            raise ValueError("Zarr dataset has not been initialized.")
        return self._dset

    @dset.setter
    def dset(self, value):
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

    def finalize(self) -> None:
        """合并 temp 文件到 Zarr。"""
        if self.mode != "w" or self._finalized:
            return
        self._finish_temp_writes()
        if self._n_ions is None or self._mz_axis is None or not self._coordinate_chunks:
            raise ValueError("No ion-image batches were written.")

        self.coordinates = np.vstack(self._coordinate_chunks)
        coordinate_index = _prepare_coordinate_index(self.coordinates, self._width, self._height)
        if coordinate_index.duplicate_count:
            logger.warning(
                "Zarr ion-image output contains duplicate 2D pixels after coordinate normalization; "
                "keeping the last spectrum for each duplicated pixel. "
                f"duplicate_count={coordinate_index.duplicate_count}, "
                f"coordinate_base={coordinate_index.coordinate_base}"
            )
        self.coordinates = coordinate_index.coordinates.astype("<u4", copy=False)
        self._n_pixels = int(self.coordinates.shape[0])
        self._width = coordinate_index.width
        self._height = coordinate_index.height
        self._flat_idx = coordinate_index.flat_idx
        self._source_idx = coordinate_index.source_idx
        self._coordinate_base = coordinate_index.coordinate_base

        self.filepath.parent.mkdir(parents=True, exist_ok=True)
        self._merge_to_zarr()
        self._finalized = True
        self._shutdown_io_executor()
        self._cleanup_temp()
        self._open_for_read()

    def _merge_to_zarr(self) -> None:
        n, h, w = self._n_ions, self.height, self.width
        npix = self._n_pixels
        dtype = self.intensity_dtype
        executor = self._io_executor
        batch_infos = self._batch_infos
        input_npix = sum(batch_info.n_pixels for batch_info in batch_infos)
        assert n is not None and executor is not None

        if self._flat_idx is None or self._source_idx is None:
            coordinate_index = _prepare_coordinate_index(self.coordinates, w, h)
            self._flat_idx = coordinate_index.flat_idx
            self._source_idx = coordinate_index.source_idx
            self._coordinate_base = coordinate_index.coordinate_base
        flat_idx = self._flat_idx
        source_idx = self._source_idx

        mmaps: list[np.ndarray] = []
        try:
            for bi in batch_infos:
                mmaps.append(np.memmap(bi.path, dtype=dtype, mode="r", shape=(n, bi.n_pixels), order="C"))

            root = zarr.open_group(str(self.filepath), mode="w")
            root.attrs.update(
                {
                    "width": int(w),
                    "height": int(h),
                    "n_ions": int(n),
                    "n_pixels": int(npix),
                    "intensity_dtype": str(np.dtype(dtype)),
                    "mz_dtype": str(np.dtype(self.mz_dtype)),
                    "coordinate_base": self._coordinate_base or "auto",
                }
            )
            _create_array(root, "mz_axis", data=np.asarray(self.mz_axis, dtype=self.mz_dtype))
            _create_array(root, "coordinates", data=np.asarray(self.coordinates, dtype="<u4"))
            dset = _create_array(
                root,
                "ion_images",
                shape=(n, h, w),
                chunks=(self._chunk_size, h, w),
                dtype=dtype,
            )

            block = max(16, self.ion_block)
            starts = list(range(0, n, block))
            if not starts:
                return

            curr = _read_scatter(
                mmaps,
                batch_infos,
                flat_idx,
                source_idx,
                starts[0],
                min(n, starts[0] + block),
                input_npix,
                h,
                w,
                dtype,
            )
            fut: Future | None = None
            if len(starts) > 1:
                s1 = starts[1]
                fut = executor.submit(
                    _read_scatter,
                    mmaps,
                    batch_infos,
                    flat_idx,
                    source_idx,
                    s1,
                    min(n, s1 + block),
                    input_npix,
                    h,
                    w,
                    dtype,
                )

            for i, start in enumerate(starts):
                end = min(n, start + block)
                dset[start:end] = curr

                nxt: Future | None = None
                if i + 2 < len(starts):
                    s2 = starts[i + 2]
                    nxt = executor.submit(
                        _read_scatter,
                        mmaps,
                        batch_infos,
                        flat_idx,
                        source_idx,
                        s2,
                        min(n, s2 + block),
                        input_npix,
                        h,
                        w,
                        dtype,
                    )
                if fut is not None:
                    curr = fut.result()
                fut = nxt
        finally:
            for mmap in mmaps:
                mmap_ext = getattr(mmap, "_mmap", None)
                if mmap_ext is not None and hasattr(mmap_ext, "close"):
                    mmap_ext.close()

    def _open_for_read(self) -> None:
        if self._group is not None:
            return
        self._group = zarr.open_group(str(self.filepath), mode="r")
        self.dset = self._group["ion_images"]
        self.mz_axis = np.asarray(self._group["mz_axis"][:])
        self.coordinates = np.asarray(self._group["coordinates"][:])
        attrs = self._group.attrs
        self.width = int(attrs["width"])
        self.height = int(attrs["height"])
        self._n_ions = int(attrs["n_ions"])
        self._n_pixels = int(attrs["n_pixels"])
        self.mz_dtype = np.dtype(attrs.get("mz_dtype", self.mz_axis.dtype))
        self.intensity_dtype = np.dtype(attrs.get("intensity_dtype", self.dset.dtype))

        coordinate_base = attrs.get("coordinate_base", None)
        coordinate_index = _prepare_coordinate_index(
            self.coordinates,
            self.width,
            self.height,
            coordinate_base=None if coordinate_base == "auto" else coordinate_base,
        )
        self._flat_idx = coordinate_index.flat_idx
        self._source_idx = coordinate_index.source_idx
        self._coordinate_base = coordinate_index.coordinate_base

    def get_ion_image(self, ion_idx: int) -> np.ndarray:
        """读取单张 ion image，返回 shape=(H, W)。"""
        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx out of range: {ion_idx}")
        return np.asarray(self.dset[ion_idx])

    def get_ion_vector(self, ion_idx: int) -> np.ndarray:
        """从稠密图像提取有效像素强度向量。"""
        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx out of range: {ion_idx}")
        return np.asarray(self.dset[ion_idx]).ravel()[self.flat_idx]

    def get_pixel_spectrum(self, pixel_idx: int) -> np.ndarray:
        """Zarr 稠密存储不支持高效像素谱读取。"""
        raise NotImplementedError("Zarr 稠密存储不支持 get_pixel_spectrum，请使用 binary 后端")

    def close(self) -> None:
        self._group = None
        self.dset = None
        self.mz_axis = None
        self._coordinates = None
        self._flat_idx = None
        self._source_idx = None
        self._coordinate_base = None
        self._shutdown_io_executor()
        if self.mode == "w" and not self._finalized:
            self._cleanup_temp()

    @classmethod
    def from_data_manager(
        cls,
        data_manager,
        filepath,
        *,
        temp_dir=None,
        mz_dtype=None,
        intensity_dtype=None,
        ion_block=256,
        temp_workers=2,
        merge_workers=2,
        chunk_size=16,
    ) -> "IonImageDataManagerZarr":
        meta = data_manager.ms.meta
        if (
            getattr(meta, "continuous", None) is not True
            or data_manager.ms.shared_mz_list is None
            or len(data_manager.ms.shared_mz_list) == 0
        ):
            raise ValueError("from_data_manager requires continuous=True metadata and a shared mz_axis.")

        return cls(
            filepath=filepath,
            mode="w",
            width=int(meta.max_count_of_pixels_x),
            height=int(meta.max_count_of_pixels_y),
            temp_dir=temp_dir,
            mz_dtype=data_manager.mz_dtype if mz_dtype is None else mz_dtype,
            intensity_dtype=data_manager.intensity_dtype if intensity_dtype is None else intensity_dtype,
            ion_block=ion_block,
            temp_workers=temp_workers,
            merge_workers=merge_workers,
            chunk_size=chunk_size,
            metadata=meta,
            shared_mz_axis=data_manager.ms.shared_mz_list,
        )


@dataclass(frozen=True, slots=True)
class _CoordinateIndex:
    width: int
    height: int
    flat_idx: np.ndarray
    source_idx: np.ndarray
    coordinates: np.ndarray
    coordinate_base: str
    duplicate_count: int


def _coord_samples(coords: np.ndarray, *, limit: int = 5) -> list[tuple[int, int]]:
    return [tuple(map(int, row[:2])) for row in coords[:limit]] # type: ignore[return-value]


def _infer_coordinate_base(xy: np.ndarray, width: int | None, height: int | None) -> str:
    if width is not None and height is not None:
        zero_valid = bool(
            np.all((xy[:, 0] >= 0) & (xy[:, 0] < width) & (xy[:, 1] >= 0) & (xy[:, 1] < height))
        )
        one_valid = bool(
            np.all((xy[:, 0] >= 1) & (xy[:, 0] <= width) & (xy[:, 1] >= 1) & (xy[:, 1] <= height))
        )
        if zero_valid and not one_valid:
            return "zero"
        if one_valid and not zero_valid:
            return "one"

    return "zero" if bool(np.any(xy == 0)) else "one"


def _prepare_coordinate_index(
    coordinates: np.ndarray,
    width: int | None,
    height: int | None,
    coordinate_base: str | None = None,
) -> _CoordinateIndex:
    coords = np.asarray(coordinates)
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("coordinates must have shape (n_pixels, >=2).")
    if coords.shape[0] == 0:
        raise ValueError("coordinates must contain at least one pixel.")

    xy = coords[:, :2].astype(np.int64, copy=False)
    base = _infer_coordinate_base(xy, width, height) if coordinate_base is None else str(coordinate_base)
    if base not in {"zero", "one"}:
        raise ValueError(f"coordinate_base must be 'zero' or 'one', got {base!r}.")

    xy0 = xy if base == "zero" else xy - 1
    if width is None:
        width = int(np.max(xy0[:, 0])) + 1
    if height is None:
        height = int(np.max(xy0[:, 1])) + 1
    width = int(width)
    height = int(height)
    if width <= 0 or height <= 0:
        raise ValueError(f"width and height must be positive, got width={width}, height={height}.")

    valid = (xy0[:, 0] >= 0) & (xy0[:, 0] < width) & (xy0[:, 1] >= 0) & (xy0[:, 1] < height)
    if not bool(np.all(valid)):
        bad = xy[~valid]
        raise ValueError(
            "coordinates contain out-of-bounds pixels: "
            f"count={bad.shape[0]}, width={width}, height={height}, "
            f"coordinate_base={base}, examples={_coord_samples(bad)}"
        )

    flat_idx = (xy0[:, 1] * width + xy0[:, 0]).astype(np.intp, copy=False)
    unique_flat_idx = np.unique(flat_idx)
    duplicate_count = int(flat_idx.size - unique_flat_idx.size)
    if duplicate_count:
        # Dense ion-image storage has one value per 2D pixel. Match the H5 scatter
        # behavior by keeping the last spectrum written for each duplicate pixel.
        _, reverse_idx = np.unique(flat_idx[::-1], return_index=True)
        source_idx = (flat_idx.size - 1 - reverse_idx).astype(np.intp, copy=False)
        source_idx.sort()
        flat_idx = flat_idx[source_idx]
        xy = xy[source_idx]
    else:
        source_idx = np.arange(flat_idx.size, dtype=np.intp)

    return _CoordinateIndex(
        width=width,
        height=height,
        flat_idx=flat_idx,
        source_idx=source_idx,
        coordinates=xy,
        coordinate_base=base,
        duplicate_count=duplicate_count,
    )


def _create_array(group, name: str, *, data=None, shape=None, chunks=None, dtype=None):
    """Create a Zarr v3 array."""
    kwargs: dict[str, Any] = {"overwrite": True}
    if data is not None:
        kwargs["data"] = data
    else:
        kwargs.update({"shape": shape, "dtype": dtype})
    if chunks is not None:
        kwargs["chunks"] = chunks
    return group.create_array(name, **kwargs)


def _read_scatter(
    mmaps: list[np.ndarray],
    batch_infos: list[_BatchInfo],
    flat_idx: np.ndarray,
    source_idx: np.ndarray | None,
    ion_start: int,
    ion_end: int,
    n_pixels: int,
    height: int,
    width: int,
    dtype: np.dtype,
) -> np.ndarray:
    """从所有 temp 文件读取一个 ion block，向量化散射为 (block, H, W)。"""
    n_block = ion_end - ion_start
    sparse = np.empty((n_block, n_pixels), dtype=dtype)
    pixel_start = 0
    for mmap, batch_info in zip(mmaps, batch_infos):
        pixel_end = pixel_start + batch_info.n_pixels
        sparse[:, pixel_start:pixel_end] = mmap[ion_start:ion_end, :]
        pixel_start = pixel_end

    if source_idx is not None and source_idx.size != n_pixels:
        sparse = sparse[:, source_idx]

    dense = np.zeros((n_block, height * width), dtype=dtype)
    dense[:, flat_idx] = sparse
    return dense.reshape((n_block, height, width))
