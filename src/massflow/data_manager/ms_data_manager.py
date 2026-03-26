"""
MS Data Management Module
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""

import os
import uuid
import random
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import Generator
import numpy as np
from massflow.tools.logger import get_logger
from massflow.module import MassSpectrumSet, Spectrum
from massflow.tools.stream_imzml_writer import ImzMLWriter

logger = get_logger("massflow.datamanager")


class MSDataManager(ABC):
    """
    Abstract base class for MS Data Manager.
    """

    def __init__(
        self,
        ms=None,
        target_mz_range=None,
        target_locs=None,
        filepath=None,
        temp_dir=None,
        max_threads: int = 4,
        mz_dtype=np.float64,
        intensity_dtype=np.float32,
    ):
        """Initialize the MS data manager."""
        # create MassSpectrumSet object if not provided
        self.ms = MassSpectrumSet() if ms is None else ms

        # Multi-threading control part
        self._max_threads = max_threads
        self._threads_executor = None

        # data type part
        self.mz_dtype = mz_dtype
        self.intensity_dtype = intensity_dtype

        # target_locs input verification
        self._init_target_locs(target_locs)
        self.target_mz_range = target_mz_range
        self.current_spectrum_num = 0

        # file swap and path control part
        self._init_file_base_path(filepath, temp_dir)
        self._writer = None


    def _init_file_base_path(self, filepath: str | None ,temp_dir: str | None):
        """Normalize external path input to suffix-free base path."""

        if filepath is None:
            base_dir = tempfile.gettempdir() if temp_dir is None else temp_dir
            os.makedirs(base_dir ,exist_ok=True)
            self.file_base_path = os.path.join(base_dir, f"temp_{uuid.uuid4().hex}")
            self.swap_filepath = self.file_base_path
        else:
            base_path = filepath
            while (next_base := os.path.splitext(base_path)[0]) != base_path:
                base_path = next_base
            self.file_base_path = base_path
            self.swap_filepath = None


    def _init_target_locs(self, target_locs):
        """Initialize target_locs with validation."""
        if target_locs is not None:
            if len(target_locs) <= 1 or (target_locs[0][0] > target_locs[1][0] or target_locs[0][1] > target_locs[1][1]):
                logger.error("target_locs must be non-empty or locs must x1<x2,y1<y2")
                raise ValueError
            self.target_locs = target_locs
        else:
            self.target_locs = ([0, 0], [float('inf'), float('inf')])


############# Data define part #############

    @property
    def imzml_filepath(self) -> str:
        """Return full .imzML path derived from the base path."""
        return f"{self.file_base_path}.imzML"

    @property
    def ibd_filepath(self) -> str:
        """
        Return full .ibd path derived from the base path.
        """
        return f"{self.file_base_path}.ibd"

    @property
    def max_threads(self) -> int:
        return self._max_threads

    @max_threads.setter
    def max_threads(self, value: int):
        if value != self._max_threads:
            self._max_threads = value
            if self._threads_executor is not None:
                logger.info(f"Recreating ThreadPoolExecutor with max_workers={self._max_threads}")
                self._threads_executor.shutdown(wait=True)
                self._threads_executor = ThreadPoolExecutor(max_workers=self._max_threads)

    @property
    def threads_executor(self) -> ThreadPoolExecutor:
        """lazy load for thread pool executor"""
        if self._threads_executor is None:
            logger.info(f"Creating ThreadPoolExecutor with max_workers={self.max_threads}")
            self._threads_executor = ThreadPoolExecutor(max_workers=self.max_threads)
        return self._threads_executor

    @threads_executor.setter
    def threads_executor(self, executor: ThreadPoolExecutor):
        """set thread pool executor"""
        self._threads_executor = executor

    @property
    def writer(self):
        """return writer."""
        if self.ms.meta is not None:
            if self._writer is None:
                logger.info("Creating ImzMLWriter from MS meta data")
                self._writer = ImzMLWriter.from_ms_meta(
                    output_filename=str(self.file_base_path),
                    meta=self.ms.meta,
                    mz_dtype=self.mz_dtype,
                    intensity_dtype=self.intensity_dtype
                )
            return self._writer
        else:
            logger.error("MS meta data is None. Please load or copy[use copy_meta method] meta data first.")
            raise ValueError("MS meta data is None. Please load or copy meta data first.")


############# Method define part #############
    @abstractmethod
    def copy_meta(self, source_dm: "MSDataManager"):
        pass

    @abstractmethod
    def preload_hook(self, *args, **kwargs):
        """Hook for pre-loading or preparing metadata before full data loading."""

    @abstractmethod
    def loading_hook(self, *args, **kwargs):
        """Hook invoked while spectra are being loaded."""

    @abstractmethod
    def loaded_hook(self):
        """Finalize metadata after loading spectra."""

    @abstractmethod
    def load_head_data(self):
        """Load data from the configured input file path."""

    @abstractmethod
    def batch_generator(self, batch_size: int = 256, max_threads: int = 0) -> Generator[list[Spectrum], None, None]:
        """Get a multi-threaded batch generator."""

    @abstractmethod
    def matrix_generator(self, batch_size: int = 256, include_mz: bool = True, max_threads: int = 0):
        """Get a multi-threaded matrix generator."""

############# Method part #############

    def inspect_data(self, inspect_target=None):
        """Inspect the data structure of the MSI object."""

        # get meta data info
        if self.ms.meta is None:
            raise ValueError("MS meta data is None. Please load meta data first.")

        meta_info = str(self.ms.meta) + "\n"
        logger.info(meta_info)

        # get base info
        base_info = "Spectrum information:\n"
        total_len = len(self.ms)
        if total_len > 0:
            # Determine index: use inspect_target if valid, else pick a random one
            is_valid_target = inspect_target is not None and 0 <= inspect_target < total_len
            idx = inspect_target if is_valid_target else random.randint(0, total_len - 1)
            label = "Target" if is_valid_target else "Random"

            base_info += f" [{label} Index {idx}]\n"
            base_info += f"{self.ms[idx]}\n" #type: ignore

        logger.info(base_info)


    def detect_sorted(self):
        """Check if the m/z arrays in the MS data are sorted."""

        total_length = len(self.ms)
        check_num = random.randint(0, total_length -1)
        check_rsult = self.ms[check_num].is_sorted()
        if not check_rsult:
            logger.info("Detected unsorted m/z data in the MS object.")
            for spectrum in self.ms:
                spectrum.sort_by_mz = False
        else:
            logger.info("random test m/z data in the MS object are sorted.")


    def create_ms_meta_mask(self):
        """Create a binary occupancy mask for all available spectra coordinates."""
        if self.ms.meta is None:
            logger.error("MS meta data is None. Please load meta data first.")
            raise ValueError("MS meta data is None. Please load meta data first.")

        if (self.ms.meta.max_count_of_pixels_x is None or self.ms.meta.max_count_of_pixels_y is None):
            logger.error("Image dimensions missing in meta data.")
            raise ValueError("Image dimensions missing in meta data.")

        # Get image dimensions
        # Get image dimensions
        width = int(self.ms.meta.max_count_of_pixels_x)
        height = int(self.ms.meta.max_count_of_pixels_y)

        # Prepare an empty occupancy mask
        mask = np.zeros((height, width), dtype=np.int8)

        # Fill the occupancy mask across all z-planes
        xs = np.fromiter((int(x) for z in self.ms.coordinate_index.values() for x, ydict in z.items() for _ in ydict),dtype=np.int16,)
        ys = np.fromiter((int(y) for z in self.ms.coordinate_index.values() for x, ydict in z.items() for y in ydict),dtype=np.int16,)

        valid = (0 <= xs) & (xs < width) & (0 <= ys) & (ys < height)
        mask[ys[valid], xs[valid]] = 1

        # Cache mask in metadata and return
        self.ms.meta.mask = mask


    def close(self):
        """
        Close any resources held by the data manager.
        """
        # Close writer first to release file locks
        self.close_writer()

        # Thread pool shutdown
        if self._threads_executor is not None:
            self._threads_executor.shutdown(wait=True)
            self._threads_executor = None

        # Temporary swap file cleanup
        swap_path_str = getattr(self, "swap_filepath", None)
        if swap_path_str:
            swap_base = Path(swap_path_str)
            suffixes = [".imzML", ".ibd", ".fdata", ".log", ".metadata", ".IBD"]
            #del all related files with the base name
            [swap_base.with_suffix(s).unlink(missing_ok=True) for s in suffixes]  #pylint: disable=expression-not-assigned


    def close_writer(self):
        """close ImzMLWriter."""
        if self._writer is not None:
            self._writer.close()


    def clear_batch_data_memory(self, batch):
        """Clear data in a batch of spectra to free memory."""
        for spec in batch:
            spec.clear_data()


    def clear_all_data_memory(self):
        """
        Clear data in all spectra to free memory.
        """
        for spec in self.ms:
            spec.clear_data()


    def swap_batch_data_out2disk(self, batch, writer):
        """Swap out a batch of spectra to disk to free memory."""
        for spec in batch:
            spec.swap_out2disk(writer=writer)

    def swap_matrix_data_out2disk(
        self,
        mz_data: np.ndarray | None,
        intensity_matrix: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray
    ):
        """swap out a matrix of spectra to disk to free memory."""
        if mz_data is None:
            logger.error("mz_data is None. Cannot swap matrix data to disk without m/z information.")
            raise ValueError("mz_data is None. Cannot swap matrix data to disk without m/z information.")

        if not np.issubdtype(coordinates.dtype, np.integer):
            coordinates = coordinates.astype(np.int32, copy=False)

        is_shared_mz = mz_data.ndim == 1

        for i,(intensity_row, length, coord) in enumerate(zip(intensity_matrix, lengths, coordinates)):
            length = int(length)

            mz = mz_data[:length] if is_shared_mz else mz_data[i, :length]
            intensity = intensity_row[:length]

            self.writer.add_spectrum(mz, intensity, tuple(coord.tolist()))

    def __enter__(self):
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


    def __del__(self):
        self.close()
