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
        ms=None,
        target_locs=None,
        filepath=None,
        max_threads: int = 8,
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
            max_threads (int): Maximum worker threads for batch loading utilities.
        """
        super().__init__(ms=ms,
                         target_mz_range=None, 
                         target_locs=target_locs, 
                         filepath=filepath, 
                         max_threads=max_threads, 
                         mz_dtype=mz_dtype, 
                         intensity_dtype=intensity_dtype, 
                         temp_dir=temp_dir)

        # Initialize ImzML parser and reader       
        if os.path.exists(self.filepath) and self.filepath.endswith(".imzML") and os.path.exists(self.filepath[:-5] + "ibd"):
            self.lazy_init()
        else:
            self.parser = None
            self.reader = None
            self.ibd_path = None
            self.ms.meta  = ImzMlMetaData(filepath=None)
            logger.warning(f"Warning: ImzML file or corresponding ibd file not found at {self.filepath} or {self.filepath[:-5] + 'ibd'}. Deferred initialization until load.")

    def lazy_init(self):
        """
        Lazy initialization placeholder for ImzML data manager.
        Currently no additional lazy init steps are required beyond base class.
        """
        # catch warnings during parser initialization
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            self.parser = ImzMLParser(self.filepath)

            if self.parser is not None :
                self.reader = self.parser.portable_spectrum_reader()
                self.ms.meta.filepath = self.filepath # type: ignore
                self.ibd_path = self.filepath[:-5] + "ibd"

            # log warnings part
            combined_message = "\r\n".join([f"{wm.message}" for wm in w])
            if len(combined_message) > 0:
                logger.warning(f"{combined_message}")


    def load_full_data_from_file(self):
        """
            Stream spectra from .imzML and add lazy SpectrumImzML placeholders into `ms`.

            Behavior:
                    - Iterates over `parser.coordinates` and creates a `SpectrumImzML` placeholder per pixel.
                    - If `target_locs` is provided, only pixels within the inclusive window are added.
                    - For continuous data, `ms.shared_mz_list` is populated once (on first seen spectrum)
                        and each placeholder gets `spectrum.mz_list = ms.shared_mz_list`.
         """

        # Ensure parser is initialized before loading defender to avoid repeated checks
        if self.parser is None:
            if  os.path.exists(self.filepath) and os.path.getsize(self.filepath) > 0 and self.filepath.lower().endswith(".imzml"):
                self.lazy_init()
                if self.parser is None:
                    logger.error(f"load data faild : Unable to initialize ImzML parser for file {self.filepath}.")
                    raise ValueError(f"load data faild : Unable to initialize ImzML parser for file {self.filepath}.")
            else:
                logger.error(f"load data faild : File {self.filepath} is empty or does not exist.")
                raise ValueError(f"load data faild : File {self.filepath} is empty or does not exist.")
        
        # Start loading process
        logger.info(f"Loading data from ImzML file: {self.filepath}")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            # pre loading part : meta data load part
            self.preload_meta()

            # spectrum data load and create part  [find by coordinate]
            for i, c in enumerate(self.parser.coordinates):
                x, y, z = c
                c1, c2 = (self.target_locs if self.target_locs is not None else ([0, 0], [999, 999]))
                if c1[0] <= x <= c2[0] and c1[1] <= y <= c2[1]:
                    # Build (x,y,z)->index mapping and add SpectrumImzML placeholders
                    spectrum = SpectrumImzML(reader=self.reader,
                                             ibd_path=self.ibd_path,
                                             index=i,
                                             coordinates=[x, y, z],
                                             shared_mz_list=self.ms.shared_mz_list)
                    # loading part : update meta data during loading
                    self.loading_meta(x=spectrum.coordinate.x, #type: ignore
                                      y=spectrum.coordinate.y, #type: ignore
                                      index=i)

                    self.ms.add_spectrum(spectrum)
                    self.current_spectrum_num += 1

            # log warnings if any
            combined_message = "\r\n".join([f"{wm.message}" for wm in w])
            if len(combined_message) > 0:
                logger.warning(f"{combined_message}")

            # loaded part : finalize meta data load part
            self.loaded_meta()

    def preload_meta(self, *args, **kwargs):
        """Pre-load and populate metadata from the ImzML parser before loading spectra."""
        self.extract_metadata()

        # assign shared mz_list for continuous data
        if self.parser and self.ms.meta.continuous and self.ms.shared_mz_list is None:  # type: ignore
            self.ms.shared_mz_list, _ = self.parser.getspectrum(1)
            logger.info(f"Assigning shared mz_list for continuous data; len is {len(self.ms.shared_mz_list)} subsequent spectra will reuse it.")

    def loading_meta(self, *args, **kwargs):
        """Update metadata during spectrum streaming (continuous mode + min pixel tracking)."""
        # update min pixel x,y
        # need to judge meta is ImzMlMetaData
        if isinstance(self.ms.meta, ImzMlMetaData):
            self.ms.meta.min_pixel_x = min(self.ms.meta.min_pixel_x, kwargs["x"])
            self.ms.meta.min_pixel_y = min(self.ms.meta.min_pixel_y, kwargs["y"])

    def extract_metadata(self):
        """Iterate `meta_index` and populate matching attributes from the parser."""

        logger.info("Extracting metadata...")

        if self.parser is None:
            logger.error("Parser is not initialized. Please set parser or filepath first.")
            raise ValueError("Parser is not initialized. Please set parser or filepath first.")

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
            self.parser.metadata.file_description,  # File description (data type, creation time, etc.) #type: ignore
            self.parser.metadata.scan_settings,  # Scan settings (scan mode, m/z range, etc.)#type: ignore
            self.parser.metadata.instrument_configurations,  # Instrument configuration (model, ion source, etc.)#type: ignore
            self.parser.metadata.samples,  # Sample information (sample name, preparation, etc.) #type: ignore
            self.parser.metadata.softwares,  # Software information (parser, version, etc.) #type: ignore
            self.parser.metadata.data_processings,  # Data processing (peak picking, normalization, etc.) #type: ignore
            self.parser.metadata.referenceable_param_groups,  # Referenceable parameter groups (shared metadata) #type: ignore
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

    def close(self):
        """
        Safely close underlying resources associated with the ImzML parser.
        """
        if self.parser is not None:
            # Try closing memory-mapped handle if present
            m = getattr(self.parser, "m", None)
            if m is not None and hasattr(m, "close"):
                m.close()

        super().close()

    def __enter__(self):
        return super().__enter__()  # father class
        return super().__enter__()  # father class

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
