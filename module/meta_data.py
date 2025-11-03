# module/msi_meta_data.py
"""
Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""
import os
from pyimzml.metadata import ParamGroup
from pyimzml.ImzMLParser import ImzMLParser
from logger import get_logger

logger = get_logger("meta_data")

meta_index = {
    "IMS:1000042": "max_count_of_pixels_x", # Image width (pixel count)
    "IMS:1000043": "max_count_of_pixels_y",# Image height (pixel count)
    "IMS:1000044": "max_dimension_x",# Image width (physical size, µm)
    "IMS:1000045": "max_dimension_y",# Image height (physical size, µm)
    "IMS:1000046": "pixel_size_x",# Pixel width (µm)
    "IMS:1000047": "pixel_size_y",# Pixel height (µm)
    "IMS:1000053": "absolute_position_offset_x",# X-axis position offset
    "IMS:1000054": "absolute_position_offset_y",# Y-axis position offset
    "IMS:1000031": "processed",# Whether the data is processed
    "MS:1000031": "instrument_model",# Instrument model
    "MS:1000127": "centroid_spectrum",# Mass spectrum in centroid mode
    "MS:1000128": "profile_spectrum",# Mass spectrum in profile mode
    "MS:1000579": "ms1_spectrum",# MS1 spectrum
    "MS:1000580": "msn_spectrum"# MSn spectrum
}

class MetaDataBase:
    """
    Abstract base class for MSI data models.

    Manages all common metadata fields, properties, and setters
    that are shared between different data representations.

    All metadata fields starting with _meta_ are automatically synchronized
    to the metadata dictionary via property setters.
    """

    def __init__(self, name, version=None, mz_num=None, storage_mode='split'):

        self._meta = {}
        # Initialize all metadata fields via properties to trigger auto-sync
        self._name = None
        self._version = None
        self._mz_num = None
        self._storage_mode = None
        self._meta_index = None

        #set actual values through properties
        self.name = name
        self.version = version
        self.mz_num = mz_num
        self.storage_mode = storage_mode
        self.meta_index = meta_index

    def _set(self, key, value):
        self._meta[key] = value

    #Metadata properties with auto-sync
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if name is not None:
            self._name = name
            self._set('name', name)

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, version):
        if version is not None:
            assert version > 0, "Version must be positive"
            self._version = version
            self._set('version', version)

    @property
    def mz_num(self):
        return self._mz_num

    @mz_num.setter
    def mz_num(self, mz_num):
        if mz_num is not None:
            self._mz_num = int(mz_num)
            self._set('mz_num', self._mz_num)

    @property
    def storage_mode(self):
        return self._storage_mode

    @storage_mode.setter
    def storage_mode(self, mode):
        if mode is not None:
            self._storage_mode = mode
            self._set('storage_mode', self._storage_mode)

    def __getitem__(self, key):
        return self._meta[key]

    def __iter__(self):
        return iter(self._meta)

    def __len__(self):
        return len(self._meta)

    def keys(self):
        return self._meta.keys()

    def items(self):
        return self._meta.items()

    def values(self):
        return self._meta.values()

    def get(self, key, default=None):
        return self._meta.get(key, default)

    def to_dict(self):
        return dict(self._meta)

    @property
    def meta_index(self):
        return self._meta_index

    @meta_index.setter
    def meta_index(self, meta_index):
        if meta_index is not None:
            if not isinstance(meta_index, dict):
                raise TypeError("meta_index must be a dict")
            self._meta_index = meta_index


class MetaDataImzMl(MetaDataBase):
    """ImzML metadata wrapper that loads and caches frequently used fields."""
    def __init__(self,
                 name="MSI",
                 version=1.0,
                 mz_num=None,
                 storage_mode='split',
                 parser: ImzMLParser = None,
                 filepath: str = None
                 ):
        """Initialize the metadata object with either a parser or a file path."""
        super().__init__(name, version, mz_num, storage_mode)

        self._filepath = None
        self._parser = None
        self._spectrum_count_num = None
        self._max_count_of_pixels_x = None
        self._max_count_of_pixels_y = None
        self._max_dimension_x = None
        self._max_dimension_y = None
        self._pixel_size_x = None
        self._pixel_size_y = None
        self._absolute_position_offset_x = None
        self._absolute_position_offset_y = None
        self._processed = None
        self._instrument_model = None
        self._centroid_spectrum = None
        self._profile_spectrum = None
        self._ms1_spectrum = None
        self._msn_spectrum = None

        # Set actual value through property
        if parser is not None:
            self.parser = parser
        elif filepath is not None:
            self.filepath = filepath
        else:
            raise ValueError("Either parser or filepath must be provided")

        if self.parser is not None:
            self.spectrum_count_num = len(self.parser.coordinates)
            self.extract_metadata()  # Use pyimzml.ImzMLParser to retrieve metadata

    @property
    def filepath(self):
        """Return the associated imzML file path."""
        return self._filepath

    @filepath.setter
    def filepath(self, filepath: str):
        """Set the imzML file path and initialize the parser if needed."""
        if not filepath or not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        self._filepath = filepath
        if self.parser is None:
            self.parser = ImzMLParser(filepath)

    @property
    def parser(self):
        """Return the bound pyimzML parser instance."""
        return self._parser

    @parser.setter
    def parser(self, parser: ImzMLParser):
        """Set the pyimzML parser instance."""
        if parser is not None and not isinstance(parser, ImzMLParser):
            raise TypeError("parser must be an instance of pyimzML.ImzMLParser or None")
        self._parser = parser

    @property
    def spectrum_count_num(self):
        """Return the number of spectra."""
        return self._spectrum_count_num

    @spectrum_count_num.setter
    def spectrum_count_num(self, spectrum_count_num: int):
        """Persist the number of spectra and sync it to the base storage."""
        self._spectrum_count_num = spectrum_count_num
        self._set('spectrum_count_num', spectrum_count_num)

    @property
    def max_count_of_pixels_x(self):
        """Return the pixel count along the X axis."""
        return self._max_count_of_pixels_x

    @max_count_of_pixels_x.setter
    def max_count_of_pixels_x(self, max_count_of_pixels_x):
        """Set the pixel count along the X axis."""
        if max_count_of_pixels_x is not None:
            self._max_count_of_pixels_x = max_count_of_pixels_x
            self._set('max_count_of_pixels_x', max_count_of_pixels_x)

    @property
    def max_count_of_pixels_y(self):
        """Return the pixel count along the Y axis."""
        return self._max_count_of_pixels_y

    @max_count_of_pixels_y.setter
    def max_count_of_pixels_y(self, max_count_of_pixels_y):
        """Set the pixel count along the Y axis."""
        if max_count_of_pixels_y is not None:
            self._max_count_of_pixels_y = max_count_of_pixels_y
            self._set('max_count_of_pixels_y', max_count_of_pixels_y)

    @property
    def max_dimension_x(self):
        """Return the physical dimension along the X axis."""
        return self._max_dimension_x

    @max_dimension_x.setter
    def max_dimension_x(self, max_dimension_x):
        """Set the physical dimension along the X axis."""
        self._max_dimension_x = max_dimension_x
        self._set('max_dimension_x', max_dimension_x)

    @property
    def max_dimension_y(self):
        """Return the physical dimension along the Y axis."""
        return self._max_dimension_y

    @max_dimension_y.setter
    def max_dimension_y(self, max_dimension_y):
        """Set the physical dimension along the Y axis."""
        self._max_dimension_y = max_dimension_y
        self._set('max_dimension_y', max_dimension_y)

    @property
    def pixel_size_x(self):
        """Return the pixel size along the X axis."""
        return self._pixel_size_x

    @pixel_size_x.setter
    def pixel_size_x(self, pixel_size_x):
        """Set the pixel size along the X axis."""
        self._pixel_size_x = pixel_size_x
        self._set('pixel_size_x', pixel_size_x)

    @property
    def pixel_size_y(self):
        """Return the pixel size along the Y axis."""
        return self._pixel_size_y

    @pixel_size_y.setter
    def pixel_size_y(self, pixel_size_y):
        """Set the pixel size along the Y axis."""
        self._pixel_size_y = pixel_size_y
        self._set('pixel_size_y', pixel_size_y)

    @property
    def absolute_position_offset_x(self):
        """Return the absolute position offset on the X axis."""
        return self._absolute_position_offset_x

    @absolute_position_offset_x.setter
    def absolute_position_offset_x(self, absolute_position_offset_x):
        """Set the absolute position offset on the X axis."""
        self._absolute_position_offset_x = absolute_position_offset_x
        self._set('absolute_position_offset_x', absolute_position_offset_x)

    @property
    def absolute_position_offset_y(self):
        """Return the absolute position offset on the Y axis."""
        return self._absolute_position_offset_y

    @absolute_position_offset_y.setter
    def absolute_position_offset_y(self, absolute_position_offset_y):
        """Set the absolute position offset on the Y axis."""
        self._absolute_position_offset_y = absolute_position_offset_y
        self._set('absolute_position_offset_y', absolute_position_offset_y)

    @property
    def processed(self):
        """Return whether the data has been processed."""
        return self._processed
    @processed.setter
    def processed(self, processed):
        """Set whether the data has been processed."""
        self._processed = processed
        self._set('processed', processed)

    @property
    def instrument_model(self):
        """Return the mass spectrometer model."""
        return self._instrument_model

    @instrument_model.setter
    def instrument_model(self, instrument_model):
        """Set the mass spectrometer model."""
        self._instrument_model = instrument_model
        self._set('instrument_model', instrument_model)

    @property
    def centroid_spectrum(self):
        """Return whether centroid spectra are present."""
        return self._centroid_spectrum

    @centroid_spectrum.setter
    def centroid_spectrum(self, centroid_spectrum):
        """Set whether centroid spectra are present."""
        self._centroid_spectrum = centroid_spectrum
        self._set('centroid_spectrum', centroid_spectrum)

    @property
    def profile_spectrum(self):
        """Return whether profile spectra are present."""
        return self._profile_spectrum

    @profile_spectrum.setter
    def profile_spectrum(self, profile_spectrum):
        """Set whether profile spectra are present."""
        self._profile_spectrum = profile_spectrum
        self._set('profile_spectrum', profile_spectrum)

    @property
    def ms1_spectrum(self):
        """Return whether MS1 spectra are present."""
        return self._ms1_spectrum

    @ms1_spectrum.setter
    def ms1_spectrum(self, ms1_spectrum):
        """Set whether MS1 spectra are present."""
        self._ms1_spectrum = ms1_spectrum
        self._set('ms1_spectrum', ms1_spectrum)

    @property
    def msn_spectrum(self):
        """Return whether MSn spectra are present."""
        return self._msn_spectrum

    @msn_spectrum.setter
    def msn_spectrum(self, msn_spectrum):
        """Set whether MSn spectra are present."""
        self._msn_spectrum = msn_spectrum
        self._set('msn_spectrum', msn_spectrum)

    def extract_metadata(self):
        """Iterate _meta_index and populate matching attributes from the parser."""

        logger.info("Extracting metadata...")

        if self.parser is None:
            logger.error("Parser is not initialized. Please set parser or filepath first.")

        for accession_id, prop_name in self._meta_index.items():
            param_value = self.find_param_by_accession_id(accession_id)
            if param_value is not None:
                setattr(self, prop_name, param_value)


    def find_param_by_accession_id(self, accession_id: str): # Use pyimzML.ImzMLParser to fetch metadata
        """Search the predefined metadata areas for the given accession identifier."""

        search_areas = [
            self.parser.metadata.file_description,  # File description (data type, creation time, etc.)
            self.parser.metadata.scan_settings,  # Scan settings (scan mode, m/z range, etc.)
            self.parser.metadata.instrument_configurations,  # Instrument configuration (model, ion source, etc.)
            self.parser.metadata.samples,  # Sample information (sample name, preparation, etc.)
            self.parser.metadata.softwares,  # Software information (parser, version, etc.)
            self.parser.metadata.data_processings,  # Data processing (peak picking, normalization, etc.)
            self.parser.metadata.referenceable_param_groups,  # Referenceable parameter groups (shared metadata)
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