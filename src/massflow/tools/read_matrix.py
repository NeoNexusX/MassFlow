# Matrix reading utility functions for get_matrix_generator

from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from typing import BinaryIO, Callable, Optional, Sequence
import numpy as np

# File layout descriptor (Immutable)
@dataclass(frozen=True, slots=True)
class IbdFileLayout:
    """Immutable layout description of a .ibd binary file."""

    ibd_path: str

    # Offset and length of each spectrum — from parser
    intensity_offsets: np.ndarray
    intensity_lengths: np.ndarray
    mz_offsets: np.ndarray
    mz_lengths: np.ndarray
    coordinates: np.ndarray

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
        intensity_offsets=np.asarray(parser.intensityOffsets, dtype=np.int64),
        intensity_lengths=np.asarray(parser.intensityLengths, dtype=np.uint32),
        mz_offsets=np.asarray(parser.mzOffsets, dtype=np.int64),
        mz_lengths=np.asarray(parser.mzLengths, dtype=np.uint32),
        coordinates=np.asarray(parser.coordinates, dtype=np.int32),
        file_intensity_dtype=_PRECISION_DTYPE.get(parser.intensityPrecision, np.dtype(np.float32)),
        file_mz_dtype=_PRECISION_DTYPE.get(parser.mzPrecision, np.dtype(np.float64)),
        out_intensity_dtype=out_intensity_dtype,
        out_mz_dtype=out_mz_dtype,
        intensity_element_bytes=_PRECISION_BYTES.get(parser.intensityPrecision, 4),
        mz_element_bytes=_PRECISION_BYTES.get(parser.mzPrecision, 8),
    )

# Index filtering
def filter_target_indices(
    coordinates: Sequence[tuple[int, int, int]] | np.ndarray,
    target_locs: Optional[tuple[list[int], list[int]]] = None,
) -> np.ndarray:
    """Filter spectrum indices by spatial window."""
    coords = np.asarray(coordinates)

    if target_locs is None:
        return np.arange(coords.shape[0], dtype=np.int64)

    (x1, y1), (x2, y2) = target_locs
    mask = (
        (coords[:, 0] >= x1)
        & (coords[:, 0] <= x2)
        & (coords[:, 1] >= y1)
        & (coords[:, 1] <= y2)
    )
    return np.nonzero(mask)[0].astype(np.int64, copy=False)

# Batch metadata
def batch_lengths(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the spectrum lengths and coordinates of a batch."""
    idx = np.asarray(batch_indices, dtype=np.intp)
    lengths = np.take(layout.intensity_lengths, idx).astype(np.uint32, copy=False)
    coords = np.take(layout.coordinates, idx, axis=0).astype(np.int32, copy=False)
    return lengths, coords


def batch_lengths_and_coords(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
) -> tuple[np.ndarray, int, np.ndarray]:
    """Calculate length array, maximum length, and coordinate matrix for a batch."""
    lengths, coords = batch_lengths(layout, batch_indices)
    return lengths, int(np.max(lengths)), coords

# Contiguity detection
def is_disk_contiguous(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
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

# Low-level reads
def read_into(f: BinaryIO, offset: int, out: np.ndarray, file_dtype: np.dtype) -> None:
    """Fill the contiguous 1-D ``out`` with values of ``file_dtype`` read at ``offset``.

    Reads straight into ``out`` when the dtypes match, so the result is always
    a writable array and no intermediate ``bytes`` copy is made.
    """
    f.seek(int(offset))
    target = out if out.dtype == file_dtype else np.empty(out.size, dtype=file_dtype)
    if f.readinto(target) != target.nbytes:
        raise OSError(f"Unexpected end of .ibd data: {out.size} values requested at offset {offset}.")
    if target is not out:
        out[...] = target


def read_fragmented(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
    executor: ThreadPoolExecutor,
    max_threads: int,
    read_one: Callable[[BinaryIO, int, int], None],
) -> None:
    """Call ``read_one(file, position_in_batch, spectrum_index)`` for every spectrum of a batch.

    The batch is split into ``max_threads`` mini-batches; each worker opens its
    own file handle so seeks never race.
    """
    if max_threads < 2:
        raise ValueError("Single-thread mode is not supported for fragmented reads. Set max_threads >= 2.")
    if executor is None:
        raise ValueError("A ThreadPoolExecutor is required for fragmented reads.")

    batch_size = len(batch_indices)
    chunk = max(1, (batch_size + max_threads - 1) // max_threads)

    def _read_chunk(start: int) -> None:
        with open(layout.ibd_path, "rb") as f:
            for position in range(start, min(start + chunk, batch_size)):
                read_one(f, position, int(batch_indices[position]))

    list(executor.map(_read_chunk, range(0, batch_size, chunk)))


# Contiguous block reading (Strategy A)
def read_contiguous_block(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
    batch_size: int,
    max_length: int,
) -> np.ndarray:
    """Read contiguous intensity data with single seek + single read and reshape to 2-D matrix."""
    out = np.empty((batch_size, max_length), dtype=layout.out_intensity_dtype)
    with open(layout.ibd_path, "rb") as f:
        read_into(f, layout.intensity_offsets[batch_indices[0]], out.reshape(-1), layout.file_intensity_dtype)
    return out

# Fragmented multi-threaded reading (Strategy B)
def read_fragmented_block(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
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
        max_threads: Number of threads (for splitting mini-batches), must be >= 2.
    """

    def _read_one(f: BinaryIO, row: int, spec_idx: int) -> None:
        length = int(layout.intensity_lengths[spec_idx])
        read_into(f, layout.intensity_offsets[spec_idx], intensity_out[row, :length], layout.file_intensity_dtype)
        if mz_out is not None:  # Processed + include_mz only
            mz_len = int(layout.mz_lengths[spec_idx])
            read_into(f, layout.mz_offsets[spec_idx], mz_out[row, :mz_len], layout.file_mz_dtype)

    read_fragmented(layout, batch_indices, executor, max_threads, _read_one)
