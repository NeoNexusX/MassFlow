"""
MS Data Management Module
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""

from __future__ import annotations

import os
import random
import tempfile
import uuid
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Generator

import numpy as np

from massflow.module import MassSpectrumSet
from massflow.tools.logger import get_logger

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


    def _init_file_base_path(self, filepath: str | None ,temp_dir: str | None):
        """Normalize external path input to suffix-free base path."""

        if filepath is None:
            base_dir = tempfile.gettempdir() if temp_dir is None else temp_dir
            os.makedirs(base_dir ,exist_ok=True)
            self.file_base_path = os.path.join(base_dir, f"temp_{uuid.uuid4().hex}")
        else:
            base_path = filepath
            while (next_base := os.path.splitext(base_path)[0]) != base_path:
                base_path = next_base
            self.file_base_path = base_path


    def _init_target_locs(self, target_locs):
        """Initialize target_locs with validation."""
        if target_locs is not None:
            if len(target_locs) <= 1 or (target_locs[0][0] > target_locs[1][0] or target_locs[0][1] > target_locs[1][1]):
                logger.error("target_locs must be non-empty or locs must x1<x2,y1<y2")
                raise ValueError("target_locs must be ([x1, y1], [x2, y2]) with x1 <= x2 and y1 <= y2.")
            self.target_locs = target_locs
        else:
            max_coord = int(np.iinfo(np.int32).max)
            self.target_locs = ([0, 0], [max_coord, max_coord])


############# Data define part #############

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

############# Method define part #############
    @abstractmethod
    def copy_meta(self, source_dm: "MSDataManager"):
        ...

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
    def matrix_generator(self, batch_size: int = 256, include_mz: bool = True, max_threads: int = 0) -> Generator[tuple, Any, None]:
        """Generate batches of spectra data in matrix format."""

    @abstractmethod
    def flat_generator(self, batch_size: int = 256, include_mz: bool = True, max_threads: int = 0) -> Generator[tuple, Any, None]:
        """Generate batches of spectra data in flat format."""

    @abstractmethod
    def write_flat_batch(
        self,
        intensity_flat: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
        mz_axis: np.ndarray | None,
    ) -> None:
        """Write one flat spectrum batch using the backend's storage format."""


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
        if total_length == 0:
            logger.warning("No spectra loaded; skipping m/z sort detection.")
            return
        check_num = random.randint(0, total_length - 1)
        if not self.ms[check_num].is_sorted():
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
        width = int(self.ms.meta.max_count_of_pixels_x)
        height = int(self.ms.meta.max_count_of_pixels_y)

        # Prepare an empty occupancy mask
        mask = np.zeros((height, width), dtype=np.int8)

        # Fill the occupancy mask across all z-planes
        xs = np.fromiter((int(x) for z in self.ms.coordinate_index.values() for x, ydict in z.items() for _ in ydict), dtype=np.int64)
        ys = np.fromiter((int(y) for z in self.ms.coordinate_index.values() for ydict in z.values() for y in ydict), dtype=np.int64)

        valid = (1 <= xs) & (xs <= width) & (1 <= ys) & (ys <= height)
        mask[ys[valid] - 1, xs[valid] - 1] = 1

        # Cache mask in metadata and return
        self.ms.meta.mask = mask


    def close(self):
        """Close resources shared by all data-manager implementations."""
        if self._threads_executor is not None:
            self._threads_executor.shutdown(wait=True)
            self._threads_executor = None


    def clear_all_data_memory(self):
        """
        Clear data in all spectra to free memory.
        """
        for spec in self.ms:
            spec.clear_data()


    def __enter__(self):
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


    def __del__(self):
        self.close()
