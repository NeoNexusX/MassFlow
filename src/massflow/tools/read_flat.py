from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import BinaryIO, Optional, Sequence
import numpy as np

from massflow.tools.read_matrix import IbdFileLayout, batch_lengths, read_fragmented, read_into

# Batch metadata for flat read
def batch_lengths_and_coords_flat(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Calculate length array, coordinate matrix, and total elements for a batch.
    """
    lengths, coords = batch_lengths(layout, batch_indices)
    return lengths, coords, int(np.sum(lengths, dtype=np.int64))

# Contiguous block reading (Strategy A) for flat
def read_contiguous_block_flat(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
    total_elements: int,
) -> np.ndarray:
    """Read contiguous intensity data with single seek + single read directly to 1D."""
    out = np.empty(total_elements, dtype=layout.out_intensity_dtype)
    with open(layout.ibd_path, "rb") as f:
        read_into(f, layout.intensity_offsets[batch_indices[0]], out, layout.file_intensity_dtype)
    return out

# Fragmented multi-threaded reading (Strategy B) for flat arrays
def read_fragmented_block_flat(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
    intensity_out: np.ndarray,
    mz_out: Optional[np.ndarray],
    offsets: Sequence[int] | np.ndarray,
    executor: ThreadPoolExecutor,
    max_threads: int,
) -> None:
    """
    Multi-threaded reading of non-contiguous intensity (and optional mz) data, in-place filling of 1D output flat arrays.

    Args:
        layout: File layout descriptor.
        batch_indices: Global indices of spectra in the file for this batch.
        intensity_out: Pre-allocated 1D intensity flat array of size total_elements.
        mz_out: Pre-allocated 1D mz flat array (for Processed + include_mz); None to skip mz reading.
        offsets: The starting index in the flat array for each spectrum.
        executor: Thread pool.
        max_threads: Number of threads (for splitting mini-batches), must be >= 2.
    """

    def _read_one(f: BinaryIO, position: int, spec_idx: int) -> None:
        start = int(offsets[position])
        length = int(layout.intensity_lengths[spec_idx])
        read_into(f, layout.intensity_offsets[spec_idx], intensity_out[start:start + length], layout.file_intensity_dtype)
        if mz_out is not None:  # Processed + include_mz only
            mz_len = int(layout.mz_lengths[spec_idx])
            read_into(f, layout.mz_offsets[spec_idx], mz_out[start:start + mz_len], layout.file_mz_dtype)

    read_fragmented(layout, batch_indices, executor, max_threads, _read_one)
