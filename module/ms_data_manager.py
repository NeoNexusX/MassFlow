"""
MS Data Management Module

Provides functions for reading/writing MS data, memory statistics, and visualization.
Supports .h5/.msi files and batch import from directories, filters by m/z range,
and generates merged or split outputs.
"""

import os
from abc import ABC, abstractmethod
import h5py
import numpy as np
from .ms_module import MS

class MSDataManager(ABC):
    """
    Abstract base class for MS Data Manager.
    """
    def __init__(self,
                 ms: MS,
                 target_mz_range=None,
                 target_locs=None,
                 filepath=None):
        self._ms = ms
        self.target_mz_range = target_mz_range
        self.target_locs = target_locs
        self.filepath = filepath
        self.current_spectrum_num = 0

    def get_ms(self) -> MS:
        """
        Get the MS object.

        Returns:
            MS: The MS object.
        """
        return self._ms

    @abstractmethod
    def load_full_data_from_file(self):
        """
        Load metadata from a file.

        Args:
            filepath (str): Path to the input file.
        """

    def inspect_data(self):
        """
        Inspect the data structure of the MSI object.

        Prints metadata shapes and queue information, including max/min m/z values,
        queue length, and count of non-empty base masks.
        """
        print("MS meta data:")
        #TODO: implement inspect meta data
        print("MS  information:")
        for spectrum in self._ms:
            print("ms len :", len(spectrum))
            print(f"mz range:{min(spectrum.mz_list)} - {max(spectrum.mz_list)}")
            print("max intensity:", max(spectrum.intensity_list))

    def _write_meta_data(self, output_path):
        pass

    def write2local(self, mode="merge", prefix="MS", output_fold=None, compression_opts=9):
        pass