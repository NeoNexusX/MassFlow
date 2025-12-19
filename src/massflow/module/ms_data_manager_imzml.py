"""
MSI Data Manager for imzML format.

Handles reading .imzML files and loading data into the MSI domain model.
Supports coordinate-based sub-region loading and a shared m/z list (shared m/z axis) for
continuous data.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""

import os
import warnings
import numpy as np
from pyimzml.metadata import ParamGroup
from pyimzml.ImzMLParser import ImzMLParser
from massflow.logger import get_logger
from massflow.module.ms_data_manager import MSDataManager
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.ms_meta_data import ImzMlMetaData
from massflow.module.mass_spectrum_set import MassSpectrumSet

logger = get_logger("ms_data_manager_imzml")


class MSDataManagerImzML(MSDataManager):
    """
    MSI Data Manager for .imzML files.

    Handles reading .imzML files, building MS objects with lazily-loaded
    spectra, optional spatial sub-region selection via `target_locs`, and
    metadata extraction into `ImzMlMetaData`.
    """

    def __init__(
        self,
        ms: MassSpectrumSet = None,
        target_locs=None,
        filepath=None,
        coordinates_zero_based: bool = True,
        max_threads: int = 10,
        temp_dir = None,
        mz_dtype = np.float64,
        intensity_dtype =  np.float32,
    ):
        """
        Initialize the ImzML data manager.

        Args:
            ms (MassSpectrumSet): Mass-spectrometry domain model instance.
            target_locs (optional): Spatial window for sub-region loading.
                Expected shape: ([x1, y1], [x2, y2]) (inclusive bounds). When None, load all.
            filepath (str, optional): Path to the .imzML file (XML header; paired .ibd is
                resolved automatically).
            coordinates_zero_based (bool): Whether coordinates exposed in metadata/MS
                are treated as zero-based.
            max_threads (int): Maximum worker threads for batch loading utilities.
        """
        super().__init__(ms, None, target_locs, filepath, max_threads=max_threads, mz_dtype=mz_dtype, intensity_dtype=intensity_dtype, temp_dir=temp_dir)

        self.coordinates_zero_based = coordinates_zero_based
        self.parser = None
        self.reader = None
        self.ibd_path = None

        if os.path.exists(self.filepath) and os.path.getsize(self.filepath) > 0:
            self._init_reader_and_parser()
        else:
            logger.info(f"File {self.filepath} is new or empty. Skipping initial parsing.")
            if self.ms.meta is None:
                self.ms.meta = ImzMlMetaData(parser=None, filepath=self.filepath, coordinates_zero_based=self.coordinates_zero_based)

    def _init_reader_and_parser(self):
        """Initialize the ImzML parser and spectrum reader."""
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.parser = ImzMLParser(self.filepath)
        
        if self.ms.meta is None:
            self.ms.meta = ImzMlMetaData(parser=self.parser, filepath=self.filepath, coordinates_zero_based=self.coordinates_zero_based)
        else:
            self.ms.meta.parser = self.parser

        ibd_path = self.filepath[:-5] + "ibd"
        if os.path.exists(ibd_path):
            self.ibd_path = ibd_path
            self.reader = self.parser.portable_spectrum_reader()
        else:
            logger.error(
                f"Error: Corresponding .ibd file for {self.filepath} does not exist."
            )
            raise FileNotFoundError(
                f"Error: Corresponding .ibd file for {self.filepath} does not exist."
            )

    def load_full_data_from_file(self):
        """
                Stream spectra from .imzML and add lazy SpectrumImzML placeholders into `ms`.

                Behavior:
                        - Iterates over `parser.coordinates` and creates a `SpectrumImzML` placeholder per pixel.
                        - If `target_locs` is provided, only pixels within the inclusive window are added.
                        - For continuous data, `ms.shared_mz_list` is populated once (on first seen spectrum)
                            and each placeholder gets `spectrum.mz_list = ms.shared_mz_list`.
        """
        if self.parser is None:
            if os.path.exists(self.filepath) and os.path.getsize(self.filepath) > 0:
                self._init_reader_and_parser()
            else:
                raise ValueError("Cannot load data: File is empty or does not exist.")
            
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # regular checks
            if not self.filepath or not os.path.exists(self.filepath):
                logger.error(f"Error: File {self.filepath} does not exist.")
                raise FileNotFoundError(f"Error: File {self.filepath} does not exist.")

            if not self.filepath.lower().endswith(".imzml"):
                logger.error(f"Error: {self.filepath} is not an .imzML file.")
                raise ValueError(f"Error: {self.filepath} is not an .imzML file.")

            # meta data load part
            self.pre_load_meta()
            # list of coords
            coords = self.parser.coordinates
            # spectrum  meta data load and create part
            for i, c in enumerate(coords):
                x, y, z = c
                c1, c2 = (
                    self.target_locs
                    if self.target_locs is not None
                    else ([0, 0], [999, 999])
                )
                if c1[0] <= x <= c2[0] and c1[1] <= y <= c2[1]:
                    # Build (x,y,z)->index mapping and add SpectrumImzML placeholders
                    spectrum = SpectrumImzML(
                        reader=self.reader,
                        ibd_path=self.ibd_path,
                        index=i,
                        coordinates=[x, y, z],
                    )

                    # judge mz_list is used or not
                    spectrum.mz_list = self.ms.shared_mz_list
                    self.ms.add_spectrum(spectrum)
                    self.loading_meta(
                        x=spectrum.coordinates.x, y=spectrum.coordinates.y, index=i
                    )
                    self.current_spectrum_num += 1

            combined_message = "\r\n".join([f"{wm.message}" for wm in w])
            if len(combined_message) > 0:
                logger.warning(f"{combined_message}")
            self.loaded_meta()

    def pre_load_meta(self, *args, **kwargs):
        """Pre-load and populate metadata from the ImzML parser before loading spectra."""
        self.extract_metadata()

    def extract_metadata(self):
        """Iterate `meta_index` and populate matching attributes from the parser."""

        logger.info("Extracting metadata...")

        if self.ms.meta.parser is None:  # type: ignore
            logger.error(
                "Parser is not initialized. Please set parser or filepath first."
            )
            raise ValueError(
                "Parser is not initialized. Please set parser or filepath first."
            )

        for accession_id, prop_name in self.ms.meta.meta_index.items():  # type: ignore
            param_value = self.find_meta_by_accession_id(accession_id)
            if param_value is not None:
                setattr(self.ms.meta, prop_name, param_value)

        logger.info("Metadata extraction completed.")

    def find_meta_by_accession_id(
        self, accession_id: str
    ):
        """Search the predefined metadata areas for the given accession identifier."""

        search_areas = [
            self.ms.meta.parser.metadata.file_description,  # File description (data type, creation time, etc.) #type: ignore
            self.ms.meta.parser.metadata.scan_settings,  # Scan settings (scan mode, m/z range, etc.)#type: ignore
            self.ms.meta.parser.metadata.instrument_configurations,  # Instrument configuration (model, ion source, etc.)#type: ignore
            self.ms.meta.parser.metadata.samples,  # Sample information (sample name, preparation, etc.) #type: ignore
            self.ms.meta.parser.metadata.softwares,  # Software information (parser, version, etc.) #type: ignore
            self.ms.meta.parser.metadata.data_processings,  # Data processing (peak picking, normalization, etc.) #type: ignore
            self.ms.meta.parser.metadata.referenceable_param_groups,  # Referenceable parameter groups (shared metadata) #type: ignore
        ]

        for area in search_areas:
            if area is None:
                continue

            result = self._search_in_area(area, accession_id)
            if result is not None:
                return result

        return None

    def _search_in_area(self, area, accession_id):
        """Search a single parameter area (ParamGroup or dict) for the given accession identifier."""
        if isinstance(area, ParamGroup):
            if accession_id in area:
                return area[accession_id]

        elif isinstance(area, dict):
            for param_group in area.values():
                if accession_id in param_group:
                    return param_group[accession_id]

        return None

    def loading_meta(self, *args, **kwargs):
        """Update metadata during spectrum streaming (continuous mode + min pixel tracking)."""

        if self.parser and self.ms.meta.continuous and self.ms.shared_mz_list is None:  # type: ignore
            logger.info(
                "Assigning shared mz_list for continuous data; subsequent spectra will reuse it."
            )
            self.ms.shared_mz_list, _ = self.parser.getspectrum(kwargs["index"])

        # update min pixel x,y
        if isinstance(self.ms.meta, ImzMlMetaData):
            self.ms.meta.min_pixel_x = min(self.ms.meta.min_pixel_x, kwargs["x"])
            self.ms.meta.min_pixel_y = min(self.ms.meta.min_pixel_y, kwargs["y"])

    def close(self):
        """
        Safely close underlying resources associated with the ImzML parser.

        Parameters:
            None

        Returns:
            None

        Raises:
            None directly. Any errors during closing are caught and logged as warnings.

        Notes:
            - This attempts to close common handles exposed by pyimzml.ImzMLParser,
              including the memory-mapped object and any file-like objects.
            - After closing, it clears local and metadata references to the parser to
              avoid accidental reuse of invalid handles.
        """

        if self.parser is not None:
            # Try closing memory-mapped handle if present
            m = getattr(self.parser, "m", None)
            if m is not None and hasattr(m, "close"):
                try:
                    m.close()
                except Exception as e:
                    logger.warning(f"Failed to close memory-mapped handle: {e}")
        super().close()

    def __enter__(self):
        return super().__enter__()  # father class

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
