"""
MSI Data Manager for imzML format.

Handles reading .imzML files and loading data into the MSI domain model.
Supports filtering by m/z range and efficient ion image extraction.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
import os
import warnings
from pyimzml.metadata import ParamGroup
from pyimzml.ImzMLParser import ImzMLParser
from logger import get_logger
from .ms_data_manager import MSDataManager
from .ms_module import MS,SpectrumImzML
from .meta_data import ImzMlMetaData

logger = get_logger("ms_data_manager_imzml")

class MSDataManagerImzML(MSDataManager):
    """
    MSI Data Manager for .imzML files.

    Handles reading .imzML files, filtering by m/z range,
    and loading data into the MSI domain model.
    """

    def __init__(self,
                 ms: MS,
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

        if not self.filepath or not os.path.exists(self.filepath):
            logger.error(f"Error: File {self.filepath} does not exist.")
            raise FileNotFoundError(f"Error: File {self.filepath} does not exist.")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.parser = ImzMLParser(self.filepath)

            combined_message = "\r\n".join([ f"{wm.message}"for wm in w])
            if len(combined_message) > 0:
                logger.warning(f"{combined_message}")


        # meta data protection and read meta data
        if self.ms.meta is None:
            self.ms.meta  = ImzMlMetaData(parser=self.parser)
        elif self.ms.meta.parser is None:
            self.ms.meta.parser = self.parser


    def load_full_data_from_file(self):
        """
        Lazy-load spectra from .imzML:
        - Build coordinate->index map
        - Add MSImzML placeholders into MS
        """
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            if not self.filepath or not os.path.exists(self.filepath):
                logger.error(f"Error: File {self.filepath} does not exist.")
                raise FileNotFoundError(f"Error: File {self.filepath} does not exist.")

            if not self.filepath.lower().endswith('.imzml'):
                logger.error(f"Error: {self.filepath} is not an .imzML file.")
                raise ValueError(f"Error: {self.filepath} is not an .imzML file.")

            self.extract_metadata()

            # Build (x,y,z)->index mapping
            coords = self.parser.coordinates  # list of tuples
            for i, c in enumerate(coords):
                x, y, z = c
                #update min   pixel x,y
                self.ms.meta.min_pixel_x = min(self.ms.meta.min_pixel_x, x)
                self.ms.meta.min_pixel_y = min(self.ms.meta.min_pixel_y, y)

                c1, c2 = self.target_locs if self.target_locs is not None else [0,0],[999,999]
                if c1[0] <= x <= c2[0] and c1[1] <= y <= c2[1]:
                    spectrum = SpectrumImzML(parser=self.parser, index=i, coordinates=[x-1, y-1, z-1])
                    self._ms.add_spectrum(spectrum)
                    self.current_spectrum_num += 1

            combined_message = "\r\n".join([ f"{wm.message}"for wm in w])
            if len(combined_message) > 0:
                logger.warning(f"{combined_message}")

    def extract_metadata(self):
        """Iterate _meta_index and populate matching attributes from the parser."""

        logger.info("Extracting metadata...")

        if self.ms.meta.parser is None:
            logger.error("Parser is not initialized. Please set parser or filepath first.")
            raise ValueError("Parser is not initialized. Please set parser or filepath first.")

        for accession_id, prop_name in self.ms.meta.meta_index.items():
            param_value = self.find_meta_by_accession_id(accession_id)
            if param_value is not None:
                setattr(self.ms.meta, prop_name, param_value)

        logger.info("Metadata extraction completed.")

    def find_meta_by_accession_id(self, accession_id: str): # Use pyimzML.ImzMLParser to fetch metadata
        """Search the predefined metadata areas for the given accession identifier."""

        search_areas = [
            self.ms.meta.parser.metadata.file_description,  # File description (data type, creation time, etc.)
            self.ms.meta.parser.metadata.scan_settings,  # Scan settings (scan mode, m/z range, etc.)
            self.ms.meta.parser.metadata.instrument_configurations,  # Instrument configuration (model, ion source, etc.)
            self.ms.meta.parser.metadata.samples,  # Sample information (sample name, preparation, etc.)
            self.ms.meta.parser.metadata.softwares,  # Software information (parser, version, etc.)
            self.ms.meta.parser.metadata.data_processings,  # Data processing (peak picking, normalization, etc.)
            self.ms.meta.parser.metadata.referenceable_param_groups,  # Referenceable parameter groups (shared metadata)
        ]

        for area in search_areas:
            if area is None:
                continue

            result = self._search_in_area(area, accession_id)
            if result is not None:
                return result

        return None

    def _search_in_area(self, area, accession_id):
        """Search a single parameter area for the given accession identifier."""
        if isinstance(area, ParamGroup):
            if accession_id in area:
                return area[accession_id]

        elif isinstance(area, dict):
            for param_group in area.values():
                if accession_id in param_group:
                    return param_group[accession_id]

        return None
