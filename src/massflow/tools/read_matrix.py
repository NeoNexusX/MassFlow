# 服务于 get_matrix_generator 的矩阵读取工具函数集

from __future__ import annotations
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence
import numpy as np

# 1. 文件布局描述符 (Immutable)
@dataclass(frozen=True, slots=True)
class IbdFileLayout:
    """
    .ibd 二进制文件的不可变布局描述。
    """

    ibd_path: str

    # 每条谱图的偏移量、长度 —— 来自 parser
    intensity_offsets: Sequence[int]
    intensity_lengths: Sequence[int]
    mz_offsets: Sequence[int]
    mz_lengths: Sequence[int]
    coordinates: Sequence[tuple[int, int, int]]

    # 磁盘上的 dtype
    file_intensity_dtype: np.dtype
    file_mz_dtype: np.dtype

    # 输出 dtype
    out_intensity_dtype: np.dtype
    out_mz_dtype: np.dtype

    # 每个元素的字节数
    intensity_element_bytes: int
    mz_element_bytes: int


# 2. 从 ImzMLParser 提取 layout
_PRECISION_BYTES: dict[str, int] = {"f": 4, "d": 8, "i": 4, "l": 8}
_PRECISION_DTYPE: dict[str, np.dtype] = {
    "f": np.float32,
    "d": np.float64,
    "i": np.int32,
    "l": np.int64,
}

def extract_ibd_layout(
    parser,
    ibd_path: str,
    out_intensity_dtype: np.dtype = np.float32,
    out_mz_dtype: np.dtype = np.float64,
) -> IbdFileLayout:
    """
    从 ImzMLParser 中一次性提取 .ibd 文件布局信息。
    """
    return IbdFileLayout(
        ibd_path=ibd_path,
        intensity_offsets=parser.intensityOffsets,
        intensity_lengths=parser.intensityLengths,
        mz_offsets=parser.mzOffsets,
        mz_lengths=parser.mzLengths,
        coordinates=parser.coordinates,
        file_intensity_dtype=_PRECISION_DTYPE.get(parser.intensityPrecision, np.float32),
        file_mz_dtype=_PRECISION_DTYPE.get(parser.mzPrecision, np.float64),
        out_intensity_dtype=out_intensity_dtype,
        out_mz_dtype=out_mz_dtype,
        intensity_element_bytes=_PRECISION_BYTES.get(parser.intensityPrecision, 4),
        mz_element_bytes=_PRECISION_BYTES.get(parser.mzPrecision, 8),
    )

# 3. 索引过滤
def filter_target_indices(
    coordinates: Sequence[tuple[int, int, int]],
    target_locs: Optional[tuple[list[int], list[int]]] = None,
) -> list[int]:
    """
    根据空间窗口过滤谱图索引。
    """
    if target_locs is None:
        return range(len(coordinates))

    (x1, y1), (x2, y2) = target_locs
    return [
        i for i, (x, y, *_) in enumerate(coordinates)
        if x1 <= x <= x2 and y1 <= y <= y2
    ]

# 4. 批次元数据
def batch_lengths_and_coords(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
) -> tuple[np.ndarray, int, np.ndarray]:
    """
    计算一个 batch 的长度数组、最大长度和坐标矩阵。
    """
    lengths = np.array(
        [layout.intensity_lengths[i] for i in batch_indices],
        dtype=np.int32,
    )
    max_len = int(np.max(lengths))
    coords = np.array(
        [layout.coordinates[i] for i in batch_indices],
        dtype=np.int32,
    )
    return lengths, max_len, coords

# 5. 连续性检测
def is_disk_contiguous(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
    max_length: int,
) -> bool:
    """
    判断 batch 中所有 intensity 数据是否在磁盘上物理连续。
    """
    n = len(batch_indices)
    if n <= 1:
        return False

    stride = max_length * layout.intensity_element_bytes
    expected = stride * (n - 1)
    actual = layout.intensity_offsets[batch_indices[-1]] - layout.intensity_offsets[batch_indices[0]]
    return expected == actual

# 6. 连续块读取 (Strategy A)
def read_contiguous_block(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
    batch_size: int,
    max_length: int,
) -> np.ndarray:
    """
    单次 seek + 单次 read 读取连续 intensity 数据并 reshape 为 2-D 矩阵。
    """
    first_offset = layout.intensity_offsets[batch_indices[0]]
    total_bytes = batch_size * max_length * layout.intensity_element_bytes

    with open(layout.ibd_path, "rb") as f:
        f.seek(first_offset)
        raw = f.read(total_bytes)

    flat = np.frombuffer(raw, dtype=layout.file_intensity_dtype)
    if layout.file_intensity_dtype != layout.out_intensity_dtype:
        flat = flat.astype(layout.out_intensity_dtype)

    return flat.reshape(batch_size, max_length)

# 7. 碎片化多线程读取 (Strategy B)
def read_fragmented_block(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
    intensity_out: np.ndarray,
    mz_out: Optional[np.ndarray],
    executor: ThreadPoolExecutor,
    max_threads: int,
) -> None:
    """
    多线程读取非连续 intensity（及可选 mz）数据，就地填充输出矩阵。

    每个线程独立打开文件句柄，负责填充矩阵的若干行。

    Args:
        layout: 文件布局描述符。
        batch_indices: 本批次谱图在文件中的全局索引。
        intensity_out: 预分配的 intensity 矩阵，形状为 (N, max_length)。
        mz_out: 预分配的 mz 矩阵（Processed + include_mz 时）；
                None 时跳过 mz 读取。
        executor: 线程池。
        max_threads: 线程数目（用于拆分 mini-batch）。
    """
    batch_size = len(batch_indices)
    chunk = max(1, (batch_size + max_threads - 1) // max_threads)

    def _read_chunk(info: tuple[Sequence[int], int]) -> None:
        indices, row_offset = info
        with open(layout.ibd_path, "rb") as f:
            for idx, spec_idx in enumerate(indices):
                row = row_offset + idx

                # intensity
                off = layout.intensity_offsets[spec_idx]
                length = layout.intensity_lengths[spec_idx]
                f.seek(off)
                raw = f.read(length * layout.intensity_element_bytes)
                arr = np.frombuffer(raw, dtype=layout.file_intensity_dtype)
                if layout.file_intensity_dtype != layout.out_intensity_dtype:
                    arr = arr.astype(layout.out_intensity_dtype)
                intensity_out[row, :length] = arr

                # mz (仅 Processed + include_mz)
                if mz_out is not None:
                    mz_off = layout.mz_offsets[spec_idx]
                    mz_len = layout.mz_lengths[spec_idx]
                    f.seek(mz_off)
                    mz_raw = f.read(mz_len * layout.mz_element_bytes)
                    mz_arr = np.frombuffer(mz_raw, dtype=layout.file_mz_dtype)
                    if layout.file_mz_dtype != layout.out_mz_dtype:
                        mz_arr = mz_arr.astype(layout.out_mz_dtype)
                    mz_out[row, :mz_len] = mz_arr

    # 拆分 mini-batches
    mini_batches: list[tuple[Sequence[int], int]] = []
    for i in range(0, batch_size, chunk):
        end = min(i + chunk, batch_size)
        mini_batches.append((batch_indices[i:end], i))

    list(executor.map(_read_chunk, mini_batches))
