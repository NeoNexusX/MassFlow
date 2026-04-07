# Matrix reading utility functions for get_matrix_generator

from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence
import numpy as np

# File layout descriptor (Immutable)
@dataclass(frozen=True, slots=True)
class IbdFileLayout:
    """Immutable layout description of a .ibd binary file."""

    ibd_path: str

    # Offset and length of each spectrum — from parser
    intensity_offsets: Sequence[int]
    intensity_lengths: Sequence[int]
    mz_offsets: Sequence[int]
    mz_lengths: Sequence[int]
    coordinates: Sequence[tuple[int, int, int]]

    # dtype on disk
    file_intensity_dtype: np.dtype
    file_mz_dtype: np.dtype

    # Output dtype
    out_intensity_dtype: np.dtype
    out_mz_dtype: np.dtype

    # Bytes per element
    intensity_element_bytes: int
    mz_element_bytes: int

# Extract layout from ImzMLParser
_PRECISION_BYTES: dict[str, int] = {"f": 4, "d": 8, "i": 4, "l": 8}
_PRECISION_DTYPE: dict[str, np.dtype] = {
    "f": np.dtype(np.float32),
    "d": np.dtype(np.float64),
    "i": np.dtype(np.int32),
    "l": np.dtype(np.int64),
}
_COMPRESSION_ACCESSIONS = {
    "MS:1000574",   # zlib compression
    "IMS:1005001",  # xz compression
    "IMS:1005002",  # lz4 compression
    "IMS:1005003",  # zstd compression
}

def extract_ibd_layout(
    parser,
    ibd_path: str,
    out_intensity_dtype: np.dtype = np.dtype(np.float32),
    out_mz_dtype: np.dtype = np.dtype(np.float64),
) -> IbdFileLayout:
    """Extract .ibd file layout information from ImzMLParser in one pass."""

    # Check compression， which is not supported for direct .ibd reading
    for group_id in (parser.mzGroupId, parser.intGroupId):
        group = parser.metadata.referenceable_param_groups.get(group_id)
        if group:
            for acc in _COMPRESSION_ACCESSIONS:
                if acc in group:
                    raise ValueError(
                        f"Compressed binary arrays detected (accession {acc}) in group '{group_id}'. "
                        f"Direct .ibd reading does not support compressed data."
                    )

    return IbdFileLayout(
        ibd_path=ibd_path,
        intensity_offsets=parser.intensityOffsets,
        intensity_lengths=parser.intensityLengths,
        mz_offsets=parser.mzOffsets,
        mz_lengths=parser.mzLengths,
        coordinates=parser.coordinates,
        file_intensity_dtype=_PRECISION_DTYPE.get(parser.intensityPrecision, np.dtype(np.float32)),
        file_mz_dtype=_PRECISION_DTYPE.get(parser.mzPrecision, np.dtype(np.float64)),
        out_intensity_dtype=out_intensity_dtype,
        out_mz_dtype=out_mz_dtype,
        intensity_element_bytes=_PRECISION_BYTES.get(parser.intensityPrecision, 4),
        mz_element_bytes=_PRECISION_BYTES.get(parser.mzPrecision, 8),
    )

# Index filtering
def filter_target_indices(
    coordinates: Sequence[tuple[int, int, int]],
    target_locs: Optional[tuple[list[int], list[int]]] = None,
) -> list[int]:
    """Filter spectrum indices by spatial window."""
    if target_locs is None:
        return list(range(len(coordinates)))

    (x1, y1), (x2, y2) = target_locs
    return [
        i for i, (x, y, *_) in enumerate(coordinates)
        if x1 <= x <= x2 and y1 <= y <= y2
    ]

# Batch metadata
def batch_lengths_and_coords(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
) -> tuple[np.ndarray, int, np.ndarray]:
    """Calculate length array, maximum length, and coordinate matrix for a batch."""
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

# Contiguity detection
def is_disk_contiguous(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
    max_length: int,
) -> bool:
    """Check if all intensity data in the batch is physically contiguous on disk."""
    n = len(batch_indices)
    if n <= 1:
        return False

    stride = max_length * layout.intensity_element_bytes
    expected = stride * (n - 1)
    actual = layout.intensity_offsets[batch_indices[-1]] - layout.intensity_offsets[batch_indices[0]]
    return expected == actual

# Contiguous block reading (Strategy A)
def read_contiguous_block(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
    batch_size: int,
    max_length: int,
) -> np.ndarray:
    """Read contiguous intensity data with single seek + single read and reshape to 2-D matrix."""
    first_offset = layout.intensity_offsets[batch_indices[0]]
    total_bytes = batch_size * max_length * layout.intensity_element_bytes

    with open(layout.ibd_path, "rb") as f:
        f.seek(first_offset)
        raw = f.read(total_bytes)

    flat = np.frombuffer(raw, dtype=layout.file_intensity_dtype)
    if layout.file_intensity_dtype != layout.out_intensity_dtype:
        flat = flat.astype(layout.out_intensity_dtype)

    return flat.reshape(batch_size, max_length)

# Fragmented multi-threaded reading (Strategy B)
def read_fragmented_block(
    layout: IbdFileLayout,
    batch_indices: Sequence[int],
    intensity_out: np.ndarray,
    mz_out: Optional[np.ndarray],
    executor: ThreadPoolExecutor,
    max_threads: int,
) -> None:
    """
    Multi-threaded reading of non-contiguous intensity (and optional mz) data, in-place filling of output matrices.

    Each thread independently opens a file handle and fills several rows of the matrix.

    Args:
        layout: File layout descriptor.
        batch_indices: Global indices of spectra in the file for this batch.
        intensity_out: Pre-allocated intensity matrix with shape (N, max_length).
        mz_out: Pre-allocated mz matrix (for Processed + include_mz);
                None to skip mz reading.
        executor: Thread pool.
        max_threads: Number of threads (for splitting mini-batches).
    """
    batch_size = len(batch_indices)
    chunk = max(1, (batch_size + max_threads - 1) // max_threads)

    def _read_chunk(info: tuple[Sequence[int], int]) -> None:
        indices, row_offset = info
        with open(layout.ibd_path, "rb") as f:
            for idx, spec_idx in enumerate(indices):
                row = row_offset + idx

                # Intensity
                off = layout.intensity_offsets[spec_idx]
                length = layout.intensity_lengths[spec_idx]
                f.seek(off)
                raw = f.read(length * layout.intensity_element_bytes)
                arr = np.frombuffer(raw, dtype=layout.file_intensity_dtype)
                if layout.file_intensity_dtype != layout.out_intensity_dtype:
                    arr = arr.astype(layout.out_intensity_dtype)
                intensity_out[row, :length] = arr

                # m/z (Processed + include_mz only)
                if mz_out is not None:
                    mz_off = layout.mz_offsets[spec_idx]
                    mz_len = layout.mz_lengths[spec_idx]
                    f.seek(mz_off)
                    mz_raw = f.read(mz_len * layout.mz_element_bytes)
                    mz_arr = np.frombuffer(mz_raw, dtype=layout.file_mz_dtype)
                    if layout.file_mz_dtype != layout.out_mz_dtype:
                        mz_arr = mz_arr.astype(layout.out_mz_dtype)
                    mz_out[row, :mz_len] = mz_arr

    # Split into mini-batches
    mini_batches: list[tuple[Sequence[int], int]] = []
    for i in range(0, batch_size, chunk):
        end = min(i + chunk, batch_size)
        mini_batches.append((batch_indices[i:end], i))

    list(executor.map(_read_chunk, mini_batches))
