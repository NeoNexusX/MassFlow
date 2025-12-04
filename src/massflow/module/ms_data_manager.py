"""
MS Data Management Module

Provides abstract interfaces and common utilities for managing MS/MSI data,
including metadata inspection, batch spectrum loading, and occupancy mask
creation over spatial coordinates.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from abc import ABC, abstractmethod
import numpy as np
from massflow.logger import get_logger
from massflow.module.ms_module import MS
import random
from concurrent.futures import ThreadPoolExecutor

logger = get_logger("ms_data_manager")


class MSDataManager(ABC):
    """
    Abstract base class for MS Data Manager.
    """
    def __init__(self,
                 ms: MS,
                 target_mz_range=None,
                 target_locs=None,
                 filepath=None,
                 max_threads: int = 10):
        """
        Initialize the MS data manager.

        Args:
            ms (MS): Mass-spectrometry domain model instance.
            target_mz_range (tuple[float, float], optional): (min_mz, max_mz) to filter peaks.
            target_locs (list[tuple], optional): List of (x,y) or (x,y,z) coordinates to load.
            filepath (str, optional): Path to the input file.
        """
        self._ms = ms
        self.target_mz_range = target_mz_range

        # Multi-threading control part
        self.max_threads = max_threads
        self._threads_executor = None

        # target_locs input verification
        if target_locs is not None:
            # target locs input
            if len(target_locs) <= 1:
                logger.error("target_locs must be non-empty")
                raise ValueError("target_locs must be non-empty")
            elif target_locs[0][0] > target_locs[1][0] or target_locs[0][1] > target_locs[1][1]:
                logger.error("locs must x1<x2,y1<y2")
                raise ValueError("locs must x1<x2,y1<y2")

        # update target_locs
        self.target_locs = target_locs
        self.filepath = filepath
        self.current_spectrum_num = 0

    @property
    def threads_executor(self) -> ThreadPoolExecutor:
        """lazy load for thread pool executor"""
        if self._threads_executor is None:
            self._threads_executor = ThreadPoolExecutor(max_workers=self.max_threads)
        return self._threads_executor

    @threads_executor.setter
    def threads_executor(self, executor: ThreadPoolExecutor):
        """set thread pool executor"""
        self._threads_executor = executor
        

    @property
    def ms(self) -> MS:
        """
        Get the MS object.

        Returns:
            MS: The MS object.
        """
        return self._ms

    @ms.setter
    def ms(self, ms: MS):
        """
        Set the MS object.

        Args:
            ms (MS): Mass-spectrometry domain model instance.
        """
        self._ms = ms

    @abstractmethod
    def load_full_data_from_file(self):
        """Load data (and usually metadata) from the configured input file path.

        Concrete subclasses must implement the actual backend-specific loading
        logic (e.g., ImzML, HDF5) and use :attr:`self.filepath` and
        :attr:`self.target_locs` / :attr:`self.target_mz_range` as needed.
        """

    def inspect_data(self,inpect_num=2):
        """
        Inspect the data structure of the MSI object.

        Log metadata shapes and queue information, including max/min m/z values,
        queue length, and count of non-empty base masks.
        """
        meta_info = "MS meta data:\r\n"
        meta_info+=(f"  target_mz_range: {self.target_mz_range}\r\n")
        meta_info+=(f"  target_locs: {self.target_locs}\r\n")
        meta_info+=(f"  filepath: {self.filepath}\r\n")
        meta_info+=(f"  current_spectrum_num: {self.current_spectrum_num}\r\n")
        
        # get meta data info
        if self.ms.meta is not None:
            for attr, value in self.ms.meta.items():
                shape = getattr(value, 'shape', None)
                if shape is not None and len(shape) > 0:
                    meta_info+=(f"  meta_{attr}: {shape}\r\n")
                else:
                    meta_info+=(f"  meta_{attr}: {value}\r\n")
            logger.info(meta_info)

        # get base info
        base_info = "MS  information:\r\n"
        pointer4num = 0
        for spectrum in self._ms:
            if pointer4num >= inpect_num:
                break
            base_info+=(f"  MS len: {len(spectrum)}\r\n"
                        f"  MS range: {min(spectrum.mz_list)} - {max(spectrum.mz_list)}\r\n"
                        f"  MS coord: {spectrum.coordinates}\r\n"
                        f"  max and min mz_list: {max(spectrum.mz_list)} - {min(spectrum.mz_list)}\r\n"
                        f"  max intensity: {max(spectrum.intensity)}\r\n\r\n")
            pointer4num += 1
        logger.info(base_info)

    def get_batch_generator(self, batch_size: int = 256, max_threads: int = -1):

        # update max_threads if needed
        if max_threads != -1 and max_threads != self.max_threads:
            self.max_threads = max_threads
            # reload thread pool executor
            if self._threads_executor is not None:
                self._threads_executor.shutdown(wait=True)
            self._threads_executor = ThreadPoolExecutor(max_workers=self.max_threads)

        total_length = len(self.ms)

        # define the spectrum loading function
        def _trigger_for_load_spectrum(sp):
            _ = sp.intensity

        # 3. Batch processing loop
        for i in range(0, total_length, batch_size):
            batch = self.ms[i: i + batch_size]
            if self._threads_executor is not None:
                self._threads_executor.map(_trigger_for_load_spectrum, batch)
            yield batch

    def detect_sorted(self):

        total_length = len(self.ms)
        check_num = random.randint(0,total_length)
        check_rsult = self.ms[check_num].is_sorted()
        if not check_rsult:
            logger.info("Detected unsorted m/z data in the MS object.")
            for spectrum in self.ms:
                spectrum.sort_by_mz = False
        else:
            logger.info("random test m/z data in the MS object are sorted.")

    def create_ms_meta_mask(self):
        """
        Create a binary occupancy mask for all available spectra coordinates.

        Parameters
        - None

        Returns
        - np.ndarray: A 2D array of shape `(height, width)` where `1` marks a pixel
          with an available spectrum and `0` otherwise.

        Raises
        - ValueError: If required metadata (`max_count_of_pixels_x/y`) is missing.
        """
        if self.ms.meta is None:
            logger.error("MS meta data is None. Please load meta data first.")
            raise ValueError("MS meta data is None. Please load meta data first.")

        if self.ms.meta.max_count_of_pixels_x is None or self.ms.meta.max_count_of_pixels_y is None:
            logger.error("Image dimensions missing in meta data.")
            raise ValueError("Image dimensions missing in meta data.")

        # Get image dimensions
        width = int(self.ms.meta.max_count_of_pixels_x)
        height = int(self.ms.meta.max_count_of_pixels_y)

        # Prepare an empty occupancy mask
        mask = np.zeros((height, width), dtype=np.int8)

        # Fill the occupancy mask across all z-planes
        for z_dict in self.ms.coordinate_index.values():
            for x, y_dict in z_dict.items():
                for y in y_dict.keys():
                    ix = int(x)
                    iy = int(y)
                    if 0 <= ix < width and 0 <= iy < height:
                        mask[iy, ix] = 1

        # Cache mask in metadata and return
        self.ms.meta.mask = mask
    
    @abstractmethod
    def pre_load_meta(self,*args, **kwargs):
        """
        Hook for pre-loading or preparing metadata before full data loading.

        Typically used to populate `ms.meta` (e.g., image size, instrument
        information) prior to streaming spectra.
        """
        pass

    @abstractmethod
    def loading_meta(self,*args, **kwargs):
        """
        Hook invoked while spectra are being loaded, to update metadata
        incrementally (e.g., min/max coordinates, shared m/z list).
        """
        pass

    def loaded_meta(self,*args, **kwargs):
        """
        Finalize metadata after loading spectra.

        Default implementation computes and stores a spatial occupancy mask
        in `ms.meta.mask` based on the current `ms.coordinate_index`.
        """
        logger.info("creating ms mask.")
        self.create_ms_meta_mask()
        self.detect_sorted()

    def close(self):
        """
        Close any resources held by the data manager.
        """
        if self._threads_executor is not None:
            self._threads_executor.shutdown(wait=True)
            self._threads_executor = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()