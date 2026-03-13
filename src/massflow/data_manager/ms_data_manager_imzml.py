"""
MSI Data Manager for imzML format.

Handles reading .imzML files and loading data into the MSI domain model.
Supports coordinate-based sub-region loading and a shared m/z list (shared m/z axis) for
continuous data.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
import numpy as np
from pyimzml.metadata import ParamGroup
from pyimzml.ImzMLParser import ImzMLParser
from massflow.tools.logger import get_logger
from massflow.data_manager.ms_data_manager import MSDataManager
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.ms_meta_data_imzml import MetaDataImzMl
from massflow.tools.funs import is_valid_file

logger = get_logger("massflow.datamanager")


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
        self._parser = None
        self._reader = None
        self.ms.meta  = MetaDataImzMl()

    @property
    def reader(self):
        """Get the ImzML spectrum reader."""
        if self._reader is None:
            if self.parser is not None:
                self._reader = self.parser.portable_spectrum_reader()
            else:
                logger.error("ImzML parser is not initialized. Cannot create reader.")
                raise ValueError("ImzML parser is not initialized. Cannot create reader.")
        return self._reader

    @property
    def parser(self):

        if self._parser is not None:
            return self._parser

        if is_valid_file(self.imzml_filepath) and is_valid_file(self.ibd_filepath):
            self._parser = ImzMLParser(self.imzml_filepath)
            logger.info(f"Initialized ImzML parser for file: {self.imzml_filepath}")
            return self._parser
        else:
            logger.error(f"ImzML file {self.imzml_filepath} or {self.ibd_filepath} does not exist or is empty.")
            raise FileNotFoundError(f"ImzML file {self.imzml_filepath} or {self.ibd_filepath} does not exist or is empty.")

    def load_head_data(self):
        """
        Stream spectra from .imzML and add lazy SpectrumImzML placeholders into `ms`.

        Behavior:
                - Iterates over `parser.coordinates` and creates a `SpectrumImzML` placeholder per pixel.
                - If `target_locs` is provided, only pixels within the inclusive window are added.
                - For continuous data, `ms.shared_mz_list` is populated once (on first seen spectrum)
                    and each placeholder gets `spectrum.mz_list = ms.shared_mz_list`.
        """
        # Start loading process
        logger.info(f"Loading data from ImzML file: {self.imzml_filepath}")


        # pre loading part : meta data load part
        self.preload_hook()

        # spectrum data load and create part  [find by coordinate]
        for i, c in enumerate(self.parser.coordinates):
            x, y, z = c
            if self.target_locs[0][0] <= x <= self.target_locs[1][0] and self.target_locs[0][1] <= y <= self.target_locs[1][1]:

                # Build (x,y,z)->index mapping and add SpectrumImzML placeholders
                spectrum = SpectrumImzML(reader=self.reader,
                                            ibd_path=self.ibd_filepath,
                                            index=i,
                                            coordinates=[x, y, z],
                                            shared_mz_list=self.ms.shared_mz_list)

                self.ms.add_spectrum(spectrum)
                self.current_spectrum_num += 1

        # loaded part : finalize meta data load part
        self.loaded_hook()

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

    def copy_meta(self, source_dm: MSDataManager):
        # 检查 source_dm.ms.meta 是否也是 MetaDataImzMl 类型（或其子类）
        if isinstance(source_dm.ms.meta, type(self.ms.meta)):
            self.ms.meta.copy(source_dm.ms.meta)
        else:
            logger.warning("Source and target metadata types do not match.")

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

####### hook part #########
    def preload_hook(self, *args, **kwargs):
        """Pre-load and populate metadata from the ImzML parser before loading spectra."""
        #load meatdata from parser and assign to ms.meta
        self.extract_metadata()

        # assign shared mz_list for continuous data
        if self.parser and self.ms.meta.continuous and self.ms.shared_mz_list is None:  # type: ignore
            self.ms.shared_mz_list, _ = self.parser.getspectrum(0)
            logger.info(f"Assigning shared mz_list for continuous data; len is {len(self.ms.shared_mz_list)} subsequent spectra will reuse it.")

    def loading_hook(self, *args, **kwargs):
        """Hook invoked while spectra are being loaded."""

    def loaded_hook(self):
        """Finalize metadata after loading spectra."""
        logger.info("creating ms mask.")
        self.create_ms_meta_mask()
        self.detect_sorted()
        self.inspect_data()
####### hook end  #########
