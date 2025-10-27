"""
MS Data Management Module

Provides functions for reading/writing MS data, memory statistics, and visualization.
Supports .h5/.msi files and batch import from directories, filters by m/z range,
and generates merged or split outputs.
"""
from abc import ABC, abstractmethod
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
            assert len(target_locs) >= 1, "target_locs must be non-empty"
            if len(target_locs) ==2:
                assert target_locs[0][0] <= target_locs[1][0] and target_locs[0][1] <= target_locs[1][1], "locs must x1<x2,y1<y2"

        # update target_locs
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

    def inspect_data(self,inpect_num=10):
        """
        Inspect the data structure of the MSI object.

        Prints metadata shapes and queue information, including max/min m/z values,
        queue length, and count of non-empty base masks.
        """
        print("MS meta data:")
        #TODO: implement inspect meta data - dlq
        print(f"ms count is {self.current_spectrum_num}" )
        print("MS  information:")
        
        pointer4num = 0
        for spectrum in self._ms:
            if pointer4num >= inpect_num:
                break
            print("ms len :", len(spectrum))
            print(f"mz range:{min(spectrum.mz_list)} - {max(spectrum.mz_list)}")
            print("max intensity:", max(spectrum.intensity))
            pointer4num += 1
