from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence
import numpy as np

from massflow.tools.read_matrix import IbdFileLayout

# Batch metadata for flat read
def batch_lengths_and_coords_flat(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Calculate length array, coordinate matrix, and total elements for a batch.
    """
    idx = np.asarray(batch_indices, dtype=np.intp)
    lengths = np.take(layout.intensity_lengths, idx).astype(np.int32, copy=False)
    coords = np.take(layout.coordinates, idx, axis=0).astype(np.int32, copy=False)
    total_elements = np.sum(lengths, dtype=np.int64)
    return lengths, coords, int(total_elements)

# Contiguous block reading (Strategy A) for flat
def read_contiguous_block_flat(
    layout: IbdFileLayout,
    batch_indices: Sequence[int] | np.ndarray,
    total_elements: int,
) -> np.ndarray:
    """Read contiguous intensity data with single seek + single read directly to 1D."""
    first_offset = layout.intensity_offsets[batch_indices[0]]
    total_bytes = total_elements * layout.intensity_element_bytes

    with open(layout.ibd_path, "rb") as f:
        f.seek(first_offset)
        raw = f.read(total_bytes)

    flat = np.frombuffer(raw, dtype=layout.file_intensity_dtype)
    if layout.file_intensity_dtype != layout.out_intensity_dtype:
        flat = flat.astype(layout.out_intensity_dtype)

    return flat

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
    if max_threads < 2:
        raise ValueError("Single-thread mode is not supported for fragmented flat reads. Set max_threads >= 2.")
    if executor is None:
        raise ValueError("A ThreadPoolExecutor is required for fragmented flat reads.")

    batch_size = len(batch_indices)
    chunk = max(1, (batch_size + max_threads - 1) // max_threads)

    def _read_chunk(info: tuple[Sequence[int] | np.ndarray, Sequence[int] | np.ndarray]) -> None:
        indices, chunk_offsets = info
        with open(layout.ibd_path, "rb") as f:
            for spec_idx, offset_start in zip(indices, chunk_offsets):
                length = layout.intensity_lengths[spec_idx]

                # Intensity
                off = layout.intensity_offsets[spec_idx]
                f.seek(off)
                raw = f.read(length * layout.intensity_element_bytes)
                arr = np.frombuffer(raw, dtype=layout.file_intensity_dtype)
                if layout.file_intensity_dtype != layout.out_intensity_dtype:
                    arr = arr.astype(layout.out_intensity_dtype)
                intensity_out[offset_start:offset_start + length] = arr

                # m/z (Processed + include_mz only)
                if mz_out is not None:
                    mz_off = layout.mz_offsets[spec_idx]
                    mz_len = layout.mz_lengths[spec_idx]
                    f.seek(mz_off)
                    mz_raw = f.read(mz_len * layout.mz_element_bytes)
                    mz_arr = np.frombuffer(mz_raw, dtype=layout.file_mz_dtype)
                    if layout.file_mz_dtype != layout.out_mz_dtype:
                        mz_arr = mz_arr.astype(layout.out_mz_dtype)
                    mz_out[offset_start:offset_start + mz_len] = mz_arr

    # Split into mini-batches
    mini_batches: list[tuple[Sequence[int] | np.ndarray, Sequence[int] | np.ndarray]] = []
    for i in range(0, batch_size, chunk):
        end = min(i + chunk, batch_size)
        mini_batches.append((batch_indices[i:end], offsets[i:end]))

    list(executor.map(_read_chunk, mini_batches))
