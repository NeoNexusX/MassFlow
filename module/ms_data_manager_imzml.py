"""
MSI Data Manager for imzML format.

Handles reading .imzML files and loading data into the MSI domain model.
Supports filtering by m/z range and efficient ion image extraction.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
import os
from typing import Optional
import warnings
from pyimzml.ImzMLParser import ImzMLParser
from logger import get_logger
from .ms_data_manager import MSDataManager
from .ms_module import MS,MSImzMLBase
from .meta_data import MetaDataBase,MetaDataImzMl

logger = get_logger("ms_data_manager_imzml")


class MSDataManagerImzML(MSDataManager):
    """
    MSI Data Manager for .imzML files.

    Handles reading .imzML files, filtering by m/z range,
    and loading data into the MSI domain model.
    """

    def __init__(self,
                 ms: MS,
                 ms_meta: Optional[MetaDataBase] = None,
                 target_locs=None,
                 filepath=None):
        """
        Initialize the ImzML data manager.
        Args:
            ms (MS): Mass-spectrometry domain model instance.
            target_mz_range (tuple[float, float], optional): (min_mz, max_mz) to filter peaks.
            target_locs (list[tuple], optional): List of (x,y) or (x,y,z) coordinates to load.
            filepath (str, optional): Path to the .imzML file.
        """
        super().__init__(ms, None, target_locs, filepath)
        self.parser = ImzMLParser(self.filepath)

        # meta data protection and read meta data
        if ms_meta is None:
            self.meta  = MetaDataImzMl(parser=self.parser)
        elif ms_meta.parser is None:
            ms_meta.parser = self.parser
            self.meta = ms_meta
        else:
            self.meta = ms_meta


    def load_full_data_from_file(self):
        """
        Lazy-load spectra from .imzML:
        - Build coordinate->index map
        - Add MSImzML placeholders into MS
        """
        if not self.filepath or not os.path.exists(self.filepath):
            logger.error(f"Error: File {self.filepath} does not exist.")
            raise FileNotFoundError(f"Error: File {self.filepath} does not exist.")
        if not self.filepath.lower().endswith('.imzml'):
            logger.error(f"Error: {self.filepath} is not an .imzML file.")
            raise ValueError(f"Error: {self.filepath} is not an .imzML file.")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            if w:
                for warning in w:
                    if "pyimzml" in str(warning.filename):
                        logger.warning(f"pyimzml warn: {warning.message}")
                    else:
                        warnings.warn(warning.message, warning.category)

        # Build (x,y,z)->index mapping
        coords = self.parser.coordinates  # list of tuples
        for i, c in enumerate(coords):
            x, y, z = c
            c1, c2 = self.target_locs if self.target_locs is not None else [0,0],[999,999]
            if c1[0] <= x <= c2[0] and c1[1] <= y <= c2[1]:
                spectrum = MSImzMLBase(parser=self.parser, index=i, coordinates=[x-1, y-1, z-1])
                self._ms.add_spectrum(spectrum)
                self.current_spectrum_num += 1
