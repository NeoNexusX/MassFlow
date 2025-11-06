"""
MS Data Management Module

Provides functions for reading/writing MS data, memory statistics, and visualization.
Supports .h5/.msi files and batch import from directories, filters by m/z range,
and generates merged or split outputs.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from abc import ABC, abstractmethod
import numpy as np
from logger import get_logger
from .ms_module import MS

logger = get_logger("ms_data_manager")


class MSDataManager(ABC):
    """
    Abstract base class for MS Data Manager.
    """
    def __init__(self,
                 ms: MS,
                 target_mz_range=None,
                 target_locs=None,
                 filepath=None):
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
        """
        Load metadata from a file.

        Args:
            filepath (str): Path to the input file.
        """

    def inspect_data(self,inpect_num=3):
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
        for attr, value in self.ms.meta.items():
            shape = getattr(value, 'shape', None)
            if shape is not None and len(shape) > 0:
                meta_info+=(f"  meta_{attr}: {shape}\r\n")
            else:
                meta_info+=(f"  meta_{attr}: {value}\r\n")
        logger.info(meta_info)

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

    def create_ms_mask(self):
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
