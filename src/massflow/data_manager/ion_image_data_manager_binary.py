"""Binary backend for ion-image data."""

from __future__ import annotations

import os
import shutil
import struct
import tempfile
import uuid
from copy import deepcopy
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import astuple, dataclass, fields
from pathlib import Path
from typing import Any, Callable

import numpy as np

from massflow.data_manager.ion_image_data_manager import IonImageDataManager
from massflow.tools.infer_spectrum_type import SpectrumType, set_spectrum_type_metadata
from massflow.tools.logger import get_logger

# Header 中每个字段都是 little-endian uint64；sentinel 作为终止标记，
# HEADER_SIZE / HEADER_STRUCT 在 IonImageBinaryHeader 定义之后自动计算。
logger = get_logger("massflow.data_manager.ion_image_data_manager_binary")
UINT32_SIZE = np.dtype("<u4").itemsize
_HEADER_SENTINEL = 0xFFFFFFFFFFFFFFFF
_TEMP_WRITE_ION_BLOCK = 2048


def _float_dtype_from_bits(bits: int, *, field_name: str) -> np.dtype:
    if bits == 16:
        return np.dtype("<f2")
    if bits == 32:
        return np.dtype("<f4")
    if bits == 64:
        return np.dtype("<f8")
    raise ValueError(f"Unsupported {field_name} bit width: {bits}")


def _float_bits(dtype: np.dtype, *, field_name: str, allowed: set[int]) -> int:
    dtype = np.dtype(dtype)
    if dtype.kind != "f":
        raise TypeError(f"{field_name} must be a floating dtype, got {dtype}")
    bits = dtype.itemsize * 8
    if bits not in allowed:
        raise ValueError(f"{field_name} must be one of {sorted(allowed)} bits, got {bits}")
    return bits


@dataclass(slots=True)
class _BatchInfo:
    path: Path
    n_pixels: int


@dataclass(frozen=True, slots=True)
class IonImageBinaryHeader:
    """
    binary 文件头。包含以下内容：
    width, height, n_ions, n_pixels, mz_dtype, intensity_dtype,
    mz_offset, coordinate_offset, intensity_offset, _sentinel。
    """

    width: int
    height: int
    n_ions: int
    n_pixels: int
    mz_dtype_bits: int
    intensity_dtype_bits: int
    mz_offset: int
    coordinate_offset: int
    intensity_offset: int
    _sentinel: int = _HEADER_SENTINEL

    @classmethod
    def from_bytes(cls, data: bytes) -> "IonImageBinaryHeader":
        """
        从 bytes 反序列化，以 sentinel 作为 header 终止标记做前后向兼容。
        """
        n_current = len(fields(cls))
        data_size = len(data)

        if data_size < 8 or data_size % 8 != 0:
            raise ValueError(f"Binary header must be a multiple of 8 bytes, got {data_size}")

        # 按 8 bytes 一组解析全部 uint64
        n_available = data_size // 8
        all_values = struct.unpack(f"<{n_available}Q", data[:data_size])

        # 扫描 sentinel 截止位
        end = n_available
        try:
            end = all_values.index(_HEADER_SENTINEL)
        except ValueError:
            pass  # 老格式无 sentinel，全部作为数据字段

        # 截止位之前的部分，不足当前字段数则补 0
        values = list(all_values[:end])
        if len(values) < n_current - 1:  # -1 因为 _sentinel 不算数据字段
            values.extend([0] * (n_current - 1 - len(values)))

        # sentinel 本身和之后的新字段
        values.append(_HEADER_SENTINEL)
        if len(values) < n_current:
            values.extend([0] * (n_current - len(values)))

        return cls(*values)

    def to_bytes(self) -> bytes:
        return HEADER_STRUCT.pack(*astuple(self))


HEADER_STRUCT = struct.Struct(f"<{len(fields(IonImageBinaryHeader))}Q")
HEADER_SIZE = HEADER_STRUCT.size


class IonImageDataManagerBinary(IonImageDataManager):
    """continuous MSI 数据的 ion-major binary 后端。

    写入时先把每个 flat batch 转置成临时 batch 文件，再在 finalize 阶段按
    ion block 合并为最终的 D[ion, pixel] 连续矩阵；读取时使用 memmap 做
    随机离子图像访问。
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
        metadata=None,
        shared_mz_axis: np.ndarray | None = None,
        compute_mean_spectrum: bool = False,
    ):
        if filepath is None:
            if mode != "w":
                raise ValueError("filepath is required in read mode")
            filepath = Path("output") / f"{uuid.uuid4().hex}.bin"
        self.filepath = Path(filepath)
        self.mode = mode
        self._width = None if width is None else int(width)
        self._height = None if height is None else int(height)
        self.mz_dtype = np.dtype(mz_dtype).newbyteorder("<")
        self.intensity_dtype = np.dtype(intensity_dtype).newbyteorder("<")
        self.ion_block = int(ion_block)
        self.temp_workers = max(1, int(temp_workers))
        self.merge_workers = max(1, int(merge_workers))
        self.header: IonImageBinaryHeader | None = None
        self._mz_axis: np.ndarray | None = None
        self._coordinates: np.ndarray | None = None
        self._intensity_matrix: np.memmap | None = None

        self._batch_infos: list[_BatchInfo] = []
        self._coordinate_chunks: list[np.ndarray] = []
        self._pending_temp_futures: set[Future[tuple[int, _BatchInfo]]] = set()
        self._io_executor: ThreadPoolExecutor | None = None
        self._batch_counter = 0
        self._n_ions: int | None = None
        self._n_pixels = 0
        self._finalized = False
        self._owns_temp_dir = False
        self._temp_root: Path | None = None

        if mode not in {"r", "w"}:
            raise ValueError("mode must be 'r' or 'w'")
        if compute_mean_spectrum:
            logger.warning(f"{type(self).__name__} does not support mean spectrum computation; ignoring option.")

        if mode == "r":
            # 读模式只解析 header，并为三个数据区建立 memmap 视图。
            self._open_for_read()
            return

        _float_bits(self.mz_dtype, field_name="mz_dtype", allowed={32, 64})
        _float_bits(self.intensity_dtype, field_name="intensity_dtype", allowed={16, 32, 64})
        self._validate_metadata(metadata, shared_mz_axis)

        base_temp = Path(temp_dir) if temp_dir is not None else Path(tempfile.gettempdir())
        base_temp.mkdir(parents=True, exist_ok=True)
        self._temp_root = base_temp / f"ion_image_binary_{uuid.uuid4().hex}"
        self._temp_root.mkdir(parents=True, exist_ok=False)
        self._owns_temp_dir = True
        self._io_executor = ThreadPoolExecutor(max_workers=max(self.temp_workers, self.merge_workers))

    @classmethod
    def from_data_manager(
        cls,
        data_manager,
        filepath: str | os.PathLike[str],
        *,
        temp_dir: str | os.PathLike[str] | None = None,
        mz_dtype: np.dtype | type | None = None,
        intensity_dtype: np.dtype | type | None = None,
        ion_block: int = 256,
        temp_workers: int = 2,
        merge_workers: int = 2,
        compute_mean_spectrum: bool = False,
    ) -> "IonImageDataManagerBinary":
        """
        从已有 MSDataManager 创建 binary writer。
        这里使用原数据的图像尺寸、dtype 和 shared m/z axis 约束，确保输出
        binary 与当前 preprocessing pipeline 的最终数据一致。
        
        仅接受 continuous=True 的元数据转换。若数据类型不符，立即报错。
        """
        meta = data_manager.ms.meta
        
        # 在创建对象前进行验证，确保元数据符合要求
        if (getattr(meta, "continuous", None) is not True or 
            data_manager.ms.shared_mz_list is None or 
            len(data_manager.ms.shared_mz_list) == 0):
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
            metadata=meta,
            shared_mz_axis=data_manager.ms.shared_mz_list,
            compute_mean_spectrum=compute_mean_spectrum,
        )

    @classmethod
    def from_flat_output_contract(
        cls,
        data_manager,
        filepath: str | os.PathLike[str],
        *,
        shared_mz_axis: np.ndarray,
        output_spectrum_type: SpectrumType = "centroid",
        temp_dir: str | os.PathLike[str] | None = None,
        mz_dtype: np.dtype | type | None = None,
        intensity_dtype: np.dtype | type | None = None,
        ion_block: int = 256,
        temp_workers: int = 2,
        merge_workers: int = 2,
        compute_mean_spectrum: bool = False,
    ) -> "IonImageDataManagerBinary":
        """
        Create an ion-image writer from the final flat pipeline output contract.

        Unlike `from_data_manager`, this is intended for flat pipelines where the
        input data manager may still be processed/variable-length, but the output
        chunk already carries the final shared m/z axis produced by alignment.
        """
        mz_axis = np.asarray(shared_mz_axis)
        if mz_axis.ndim != 1 or mz_axis.size == 0:
            raise ValueError("shared_mz_axis must be a non-empty 1-D array.")

        source_meta = data_manager.ms.meta
        output_meta = deepcopy(source_meta)
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
            intensity_dtype=data_manager.intensity_dtype if intensity_dtype is None else intensity_dtype,
            ion_block=ion_block,
            temp_workers=temp_workers,
            merge_workers=merge_workers,
            metadata=output_meta,
            shared_mz_axis=mz_axis,
            compute_mean_spectrum=compute_mean_spectrum,
        )

    def _validate_metadata(self, metadata, shared_mz_axis: np.ndarray | None) -> None:
        """写出前做 metadata 级别的快速失败检查。"""
        if metadata is None:
            return
        if getattr(metadata, "continuous", None) is not True:
            raise ValueError("ion_image_binary output requires continuous=True metadata.")
        if getattr(metadata, "centroid_spectrum", None) is False:
            raise ValueError("ion_image_binary output requires centroid_spectrum not to be False.")
        if shared_mz_axis is None or len(shared_mz_axis) == 0:
            raise ValueError("ion_image_binary output requires a shared mz_axis.")

    @property
    def n_ions(self) -> int:
        if self._n_ions is not None:
            return self._n_ions
        if self.header is None:
            raise ValueError("n_ions is unknown before the first batch or read header.")
        return self.header.n_ions

    @property
    def n_pixels(self) -> int:
        if self.header is not None:
            return self.header.n_pixels
        return self._n_pixels

    @property
    def width(self) -> int:
        if self._width is None:
            if self.mode == "r":
                self._open_for_read()
            if self._width is None:
                raise ValueError("width has not been initialized.")
        return self._width

    @width.setter
    def width(self, value: int | None):
        self._width = value

    @property
    def height(self) -> int:
        if self._height is None:
            if self.mode == "r":
                self._open_for_read()
            if self._height is None:
                raise ValueError("height has not been initialized.")
        return self._height

    @height.setter
    def height(self, value: int | None):
        self._height = value

    @property
    def coordinates(self) -> np.ndarray:
        if self._coordinates is None:
            if self.mode == "r":
                self._open_for_read()
            if self._coordinates is None:
                raise ValueError("coordinates have not been initialized.")
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value: np.ndarray):
        self._coordinates = value

    @property
    def intensity_matrix(self) -> np.memmap:
        if self._intensity_matrix is None:
            if self.mode == "r":
                self._open_for_read()
            if self._intensity_matrix is None:
                raise ValueError("intensity_matrix has not been initialized.")
        return self._intensity_matrix

    @intensity_matrix.setter
    def intensity_matrix(self, value: np.memmap | None):
        self._intensity_matrix = value

    @property
    def mz_axis(self) -> np.ndarray:
        if self._mz_axis is None:
            if self.mode == "r":
                self._open_for_read()
            if self._mz_axis is None:
                raise ValueError("mz_axis has not been initialized.")
        return self._mz_axis

    @mz_axis.setter
    def mz_axis(self, value: np.ndarray | None):
        self._mz_axis = value

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
        if self._temp_root is None:
            raise RuntimeError("Temporary directory is not initialized.")

        lengths = np.asarray(lengths)
        coordinates = np.asarray(coordinates)
        intensity_flat = np.asarray(intensity_flat, dtype=self.intensity_dtype)
        mz_axis = np.asarray(mz_axis, dtype=self.mz_dtype)

        if lengths.ndim != 1 or intensity_flat.ndim != 1:
            raise ValueError("lengths and intensity_flat must be 1-D arrays.")
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError("coordinates must have shape (n_pixels, >=2).")

        batch_pixels = int(lengths.size)
        if batch_pixels == 0:
            return

        if self._n_ions is None:
            # 第一个 batch 决定全局 ion 数量和 shared m/z axis。
            self._n_ions = int(mz_axis.size)
            self.mz_axis = mz_axis.copy()
        elif int(mz_axis.size) != self._n_ions or (
            (self._mz_axis is None) or not np.array_equal(mz_axis, self._mz_axis)
        ):
            raise ValueError("All batches must use the same shared mz_axis.")

        # binary ion image 要求所有像素共享同一个 m/z axis，因此每条谱长度
        # 必须等于 n_ions；这也是后续 reshape 和转置成立的前提。
        if int(lengths.min()) != self._n_ions or int(lengths.max()) != self._n_ions:
            raise ValueError("ion_image_binary requires every spectrum length to equal n_ions.")
        if coordinates.shape[0] != batch_pixels:
            raise ValueError("coordinates row count must equal the batch pixel count.")
        if intensity_flat.size != batch_pixels * self._n_ions:
            raise ValueError("intensity_flat size must equal batch_pixels * n_ions.")

        batch_index = self._batch_counter
        self._batch_counter += 1
        batch_path = self._temp_root / f"batch_{batch_index:05d}.bin"
        if self._io_executor is None:
            raise RuntimeError("I/O executor is not initialized.")

        # 方案只保存 xy，不保存 z；坐标等到 finalize 时一次性写入最终文件。
        xy = coordinates[:, :2].astype("<u4", copy=False)
        self._coordinate_chunks.append(xy.copy())
        self._n_pixels += batch_pixels

        future = self._io_executor.submit(
            self._write_temp_batch,
            batch_index,
            batch_path,
            intensity_flat,
            batch_pixels,
            self._n_ions,
            self.intensity_dtype,
        )
        self._pending_temp_futures.add(future)
        if len(self._pending_temp_futures) >= self.temp_workers:
            self._collect_temp_futures(return_when=FIRST_COMPLETED)

    @staticmethod
    def _write_temp_batch(
        batch_index: int,
        batch_path: Path,
        intensity_flat: np.ndarray,
        batch_pixels: int,
        n_ions: int,
        dtype: np.dtype,
    ) -> tuple[int, _BatchInfo]:
        """把一个 pixel-major batch 分块转置为 [n_ions, batch_pixels] 写出。"""
        matrix = np.asarray(intensity_flat, dtype=dtype).reshape(batch_pixels, n_ions)
        with batch_path.open("wb") as temp_file:
            for ion_start in range(0, n_ions, _TEMP_WRITE_ION_BLOCK):
                ion_end = min(n_ions, ion_start + _TEMP_WRITE_ION_BLOCK)
                block = np.ascontiguousarray(matrix[:, ion_start:ion_end].T)
                temp_file.write(memoryview(block))

        return batch_index, _BatchInfo(path=batch_path, n_pixels=batch_pixels)

    def _collect_temp_futures(self, *, return_when) -> None:
        """收集已完成的临时写出任务，并保持异常能在主线程暴露。"""
        def add_batch_info(result: tuple[int, _BatchInfo]) -> None:
            _batch_index, batch_info = result
            self._batch_infos.append(batch_info)

        self._pending_temp_futures = self._collect_futures(
            self._pending_temp_futures,
            return_when=return_when,
            on_result=add_batch_info,
        )

    @staticmethod
    def _collect_futures(
        pending: set[Future[Any]],
        *,
        return_when: str,
        on_result: Callable[[Any], None] | None = None,
    ) -> set[Future[Any]]:
        """等待一批 future，抛出异常并可消费已完成结果。"""
        if not pending:
            return set()
        done, remaining = wait(pending, return_when=return_when)
        for future in done:
            result = future.result()
            if on_result is not None:
                on_result(result)
        return set(remaining)

    def _finish_temp_writes(self) -> None:
        """等待所有后台临时 batch 写出完成，并按 batch index 排序。"""
        while self._pending_temp_futures:
            self._collect_temp_futures(return_when=FIRST_COMPLETED)
        self._batch_infos.sort(key=lambda info: info.path.name)

    def finalize(self) -> None:
        """合并所有临时 batch，写出最终 binary 文件。"""
        if self.mode != "w":
            return
        if self._finalized:
            return
        self._finish_temp_writes()
        if self._n_ions is None or self._mz_axis is None or not self._coordinate_chunks:
            raise ValueError("No ion-image batches were written.")
        self._coordinates = np.vstack(self._coordinate_chunks)
        if self._width is None:
            self._width = int(np.max(self._coordinates[:, 0]))
        if self._height is None:
            self._height = int(np.max(self._coordinates[:, 1]))

        mz_bits = _float_bits(self.mz_dtype, field_name="mz_dtype", allowed={32, 64})
        intensity_bits = _float_bits(self.intensity_dtype, field_name="intensity_dtype", allowed={16, 32, 64})
        # 数据区 offset 完全按文档公式计算：header -> mz_axis -> coordinates
        # -> intensity_matrix。
        mz_offset = HEADER_SIZE
        coordinate_offset = mz_offset + self._n_ions * self.mz_dtype.itemsize
        intensity_offset = coordinate_offset + self._n_pixels * 2 * UINT32_SIZE
        self.header = IonImageBinaryHeader(
            width=int(self._width),
            height=int(self._height),
            n_ions=self._n_ions,
            n_pixels=self._n_pixels,
            mz_dtype_bits=mz_bits,
            intensity_dtype_bits=intensity_bits,
            mz_offset=mz_offset,
            coordinate_offset=coordinate_offset,
            intensity_offset=intensity_offset,
        )

        self.filepath.parent.mkdir(parents=True, exist_ok=True)
        final_size = intensity_offset + self._n_ions * self._n_pixels * self.intensity_dtype.itemsize

        fd = os.open(self.filepath, os.O_RDWR | os.O_CREAT | os.O_TRUNC, 0o644)
        batch_fds: list[int] = []
        try:
            # header、mz_axis、coordinates 都是一次性写入；大矩阵区域先
            # ftruncate 预分配，再由 _merge_intensity 分块 pwrite。
            os.write(fd, self.header.to_bytes())
            os.pwrite(fd, self.mz_axis.astype(self.mz_dtype, copy=False).tobytes(order="C"), mz_offset)
            os.pwrite(fd, self._coordinates.astype("<u4", copy=False).tobytes(order="C"), coordinate_offset)
            os.ftruncate(fd, final_size)

            batch_fds = [os.open(info.path, os.O_RDONLY) for info in self._batch_infos]
            self._merge_intensity(fd, batch_fds, intensity_offset)
        finally:
            for batch_fd in batch_fds:
                os.close(batch_fd)
            os.close(fd)

        self._finalized = True
        self._shutdown_io_executor()
        self._cleanup_temp()
        self._open_for_read()

    def _merge_intensity(self, output_fd: int, batch_fds: list[int], intensity_offset: int) -> None:
        """
        按 ion block 合并临时文件到最终 intensity_matrix。
        每个临时文件内部布局是 [ion, pixel_in_batch]。这里使用 os.pread
        显式指定读取 offset，最终写入使用 os.pwrite，避免共享文件指针带来
        的 seek/read 竞争，也为后续并发 merge 留出空间。
        """
        itemsize = self.intensity_dtype.itemsize
        block_size = max(16, self.ion_block)
        n_ions = self._n_ions
        assert n_ions is not None

        def merge_ion_block(ion_start: int) -> None:
            ion_end = min(n_ions, ion_start + block_size)
            n_block = ion_end - ion_start
            band = np.empty((n_block, self._n_pixels), dtype=self.intensity_dtype)

            pixel_start = 0
            for batch_fd, info in zip(batch_fds, self._batch_infos):
                # 对同一个 ion block，从每个 batch 文件读取对应 band，然后
                # 横向拼成 [ion_block, all_pixels]。
                nbytes = n_block * info.n_pixels * itemsize
                offset = ion_start * info.n_pixels * itemsize
                raw = os.pread(batch_fd, nbytes, offset)
                if len(raw) != nbytes:
                    raise OSError(f"Short read while merging {info.path}: expected {nbytes}, got {len(raw)}")
                pixel_end = pixel_start + info.n_pixels
                band[:, pixel_start:pixel_end] = np.frombuffer(raw, dtype=self.intensity_dtype).reshape(
                    n_block,
                    info.n_pixels,
                )
                pixel_start = pixel_end

            output_offset = intensity_offset + ion_start * self._n_pixels * itemsize
            # band 已经是 C-contiguous；memoryview 避免 band.tobytes() 再复制一份
            # 大块 bytes，和 benchmark 脚本的写法保持一致。
            os.pwrite(output_fd, memoryview(np.ascontiguousarray(band)), output_offset)

        ion_starts = range(0, n_ions, block_size)
        if self._io_executor is None:
            raise RuntimeError("I/O executor is not initialized.")

        pending_merge_futures: set[Future[None]] = set()
        for ion_start in ion_starts:
            pending_merge_futures.add(self._io_executor.submit(merge_ion_block, ion_start))
            if len(pending_merge_futures) >= self.merge_workers:
                pending_merge_futures = self._collect_futures(
                    pending_merge_futures,
                    return_when=FIRST_COMPLETED,
                )

        while pending_merge_futures:
            pending_merge_futures = self._collect_futures(
                pending_merge_futures,
                return_when=FIRST_COMPLETED,
            )

    def _open_for_read(self) -> None:
        """解析 header,并把三个数据区映射成 NumPy 视图。"""
        with self.filepath.open("rb") as f:
            self.header = IonImageBinaryHeader.from_bytes(f.read(HEADER_SIZE))

        self.width = self.header.width
        self.height = self.header.height
        self._n_ions = self.header.n_ions
        self._n_pixels = self.header.n_pixels
        self.mz_dtype = _float_dtype_from_bits(self.header.mz_dtype_bits, field_name="mz_dtype")
        self.intensity_dtype = _float_dtype_from_bits(self.header.intensity_dtype_bits, field_name="intensity_dtype")
        self.mz_axis = np.memmap(
            self.filepath,
            dtype=self.mz_dtype,
            mode="r",
            offset=self.header.mz_offset,
            shape=(self.header.n_ions,),
        )
        self._coordinates = np.memmap(
            self.filepath,
            dtype="<u4",
            mode="r",
            offset=self.header.coordinate_offset,
            shape=(self.header.n_pixels, 2),
        )
        self._intensity_matrix = np.memmap(
            self.filepath,
            dtype=self.intensity_dtype,
            mode="r",
            offset=self.header.intensity_offset,
            shape=(self.header.n_ions, self.header.n_pixels),
            order="C",
        )

    def get_ion_image(self, ion_idx: int) -> np.ndarray:
        """返回二维离子图像。

        文件内强度是按写入像素顺序保存的向量；这里根据 coordinates 把该
        ion 的强度散布回 height x width 的图像数组。imzML 常见坐标为
        1-based，因此在未检测到 0 坐标时会自动减 1。
        """

        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx out of range: {ion_idx}")

        image = np.zeros((self.height, self.width), dtype=self.intensity_dtype)
        coords = np.asarray(self.coordinates)
        xs = coords[:, 0].astype(np.int64, copy=False)
        ys = coords[:, 1].astype(np.int64, copy=False)

        zero_based = (xs.size > 0 and int(xs.min()) == 0) or (ys.size > 0 and int(ys.min()) == 0)
        if not zero_based:
            xs = xs - 1
            ys = ys - 1

        valid = (xs >= 0) & (xs < self.width) & (ys >= 0) & (ys < self.height)
        image[ys[valid], xs[valid]] = np.asarray(self.intensity_matrix[ion_idx])[valid]
        return image

    def get_ion_vector(self, ion_idx: int) -> np.ndarray:
        """返回某个 ion 在所有写入像素上的一维强度向量。"""
        if ion_idx < 0 or ion_idx >= self.n_ions:
            raise IndexError(f"ion_idx out of range: {ion_idx}")
        return np.asarray(self.intensity_matrix[ion_idx])

    def get_pixel_spectrum(self, pixel_idx: int) -> np.ndarray:
        """返回某个像素的完整强度谱，即 D[:, pixel_idx]。"""
        if pixel_idx < 0 or pixel_idx >= self.n_pixels:
            raise IndexError(f"pixel_idx out of range: {pixel_idx}")
        return np.asarray(self.intensity_matrix[:, pixel_idx])

    def close(self) -> None:
        self.mz_axis = None
        self._coordinates = None
        self._intensity_matrix = None
        self._shutdown_io_executor()
        if self.mode == "w" and not self._finalized:
            self._cleanup_temp()

    def _shutdown_io_executor(self) -> None:
        if self._io_executor is not None:
            self._io_executor.shutdown(wait=True)
            self._io_executor = None

    def _cleanup_temp(self) -> None:
        if self._owns_temp_dir and self._temp_root is not None:
            shutil.rmtree(self._temp_root, ignore_errors=True)
            self._temp_root = None

    def __del__(self):
        self.close()
