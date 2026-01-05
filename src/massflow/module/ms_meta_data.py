# module/ms_meta_data.py
"""
Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""
import os
from typing import Optional
import numpy as np
from massflow.logger import get_logger

logger = get_logger("meta_data")

meta_index = {
    # --- Image Size and Coordinates ---
    "IMS:1000042": "max_count_of_pixels_x",  # Image width (pixel count)
    "IMS:1000043": "max_count_of_pixels_y",  # Image height (pixel count)
    "IMS:1000046": "pixel_size_x",  # Pixel width (µm)
    "IMS:1000047": "pixel_size_y",  # Pixel height (µm)
    "IMS:1000053": "absolute_position_offset_x",  # X-axis position offset
    "IMS:1000054": "absolute_position_offset_y",  # Y-axis position offset
    # --- Files and Data Formats ---
    "IMS:1000030": "continuous",  # Continuous data
    "IMS:1000031": "processed",  # Whether the data is processed
    # --- Acquisition Parameters (Imaging MS Ontology) ---
    "IMS:1000041": "scan_type",
    "IMS:1000048": "scan_direction",
    "IMS:1000049": "line_scan_direction",
    "IMS:1000050": "scan_pattern",
    # --- Mass Spectrometers and Spectral Information (PSI-MS Ontology) ---
    "MS:1000031": "instrument_model",  # Instrument model
    "MS:1000127": "centroid_spectrum",  # Mass spectrum in centroid mode
    "MS:1000128": "profile_spectrum",  # Mass spectrum in profile mode
    "MS:1000579": "ms1_spectrum",  # MS1 spectrum
    "MS:1000580": "msn_spectrum",  # MSn spectrum
}


class MetaDataBase:
    """
    Abstract base class for MSI data models.

    Manages all common metadata fields, properties, and setters
    that are shared between different data representations.

    All metadata fields starting with _meta_ are automatically synchronized
    to the metadata dictionary via property setters.
    """

    def __init__(
        self,
        name: str = "default",
        version: float = 1.0,
        storage_mode: str = "split",
        max_count_of_pixels_x: Optional[int] = None,
        max_count_of_pixels_y: Optional[int] = None,
        pixel_size_x: Optional[float] = None,
        pixel_size_y: Optional[float] = None,
        mask: Optional[np.ndarray] = None,
    ):

        self._meta = {}
        # Initialize all metadata fields via properties to trigger auto-sync
        self._name = None
        self._version = None
        self._storage_mode = None
        self._meta_index = None
        self._max_count_of_pixels_x = 0
        self._max_count_of_pixels_y = 0
        self._pixel_size_x = None
        self._pixel_size_y = None
        self._processed = None
        self._continuous = None
        self._centroid_spectrum = None
        self._profile_spectrum = None
        self._peakpick = None
        self._mask = None
        self._scan_direction = None
        self._line_scan_direction = None
        self._scan_pattern = None
        self._scan_type = None

        # set actual values through properties
        self.name = name
        self.version = version
        self.storage_mode = storage_mode
        self.meta_index = meta_index
        self.max_count_of_pixels_x = max_count_of_pixels_x
        self.max_count_of_pixels_y = max_count_of_pixels_y
        self.pixel_size_x = pixel_size_x
        self.pixel_size_y = pixel_size_y
        self.mask = mask

    def _set(self, key, value):
        self._meta[key] = value

    @property
    def mask(self):
        """Return the mask of the data model."""
        return self._mask

    @mask.setter
    def mask(self, mask: Optional[np.ndarray]):
        """Set the mask of the data model."""
        if mask is None:
            return
        if not isinstance(mask, np.ndarray):
            logger.error("mask must be a numpy.ndarray")
            raise TypeError("mask must be a numpy.ndarray")
        if mask.ndim != 2:
            logger.error(f"mask must be 2D, got shape {mask.shape}")
            raise TypeError("mask must be a 2D array")
        self._mask = mask
        self._set("mask", mask)

    # Metadata properties with auto-sync
    @property
    def name(self):
        """Return the name of the data model."""
        return self._name

    @name.setter
    def name(self, name):
        if name is not None:
            self._name = name
            self._set("name", name)

    @property
    def version(self):
        """Return the version of the data model."""
        return self._version

    @version.setter
    def version(self, version):
        if version is not None:
            if version <= 0:
                logger.error("Version must be positive")
                raise ValueError("Version must be positive")

            assert version > 0, "Version must be positive"
            self._version = version
            self._set("version", version)

    @property
    def storage_mode(self):
        """Return the storage mode of the data model."""
        return self._storage_mode

    @storage_mode.setter
    def storage_mode(self, mode):
        if mode is not None:
            self._storage_mode = mode
            self._set("storage_mode", self._storage_mode)

    @property
    def max_count_of_pixels_x(self):
        """Return the pixel count along the X axis."""
        return self._max_count_of_pixels_x

    @max_count_of_pixels_x.setter
    def max_count_of_pixels_x(self, max_count_of_pixels_x):
        """Set the pixel count along the X axis."""
        if max_count_of_pixels_x is not None:
            self._max_count_of_pixels_x = max_count_of_pixels_x
            self._set("max_count_of_pixels_x", max_count_of_pixels_x)

    @property
    def max_count_of_pixels_y(self):
        """Return the pixel count along the Y axis."""
        return self._max_count_of_pixels_y

    @max_count_of_pixels_y.setter
    def max_count_of_pixels_y(self, max_count_of_pixels_y):
        """Set the pixel count along the Y axis."""
        if max_count_of_pixels_y is not None:
            self._max_count_of_pixels_y = max_count_of_pixels_y
            self._set("max_count_of_pixels_y", max_count_of_pixels_y)

    @property
    def pixel_size_x(self):
        """Return the pixel size along the X axis."""
        return self._pixel_size_x

    @pixel_size_x.setter
    def pixel_size_x(self, pixel_size_x):
        """Set the pixel size along the X axis."""
        if pixel_size_x is not None:
            self._pixel_size_x = pixel_size_x
            self._set("pixel_size_x", pixel_size_x)

    @property
    def pixel_size_y(self):
        """Return the pixel size along the Y axis."""
        return self._pixel_size_y

    @pixel_size_y.setter
    def pixel_size_y(self, pixel_size_y):
        """Set the pixel size along the Y axis."""
        if pixel_size_y is not None:
            self._pixel_size_y = pixel_size_y
            self._set("pixel_size_y", pixel_size_y)

    @property
    def continuous(self):
        """Return whether the data is continuous."""
        return self._continuous

    @continuous.setter
    def continuous(self, continuous):
        """Set whether the data is continuous."""
        self._continuous = continuous
        self._set("continuous", continuous)

    @property
    def processed(self):
        """Return whether the data has been processed."""
        return self._processed

    @processed.setter
    def processed(self, processed):
        """Set whether the data has been processed."""
        self._processed = processed
        self._set("processed", processed)

    @property
    def peakpick(self):
        """Return whether peak picking has been performed."""
        return self._peakpick

    @peakpick.setter
    def peakpick(self, peakpick):
        """Set whether peak picking has been performed."""
        self._peakpick = peakpick
        self._set("peakpick", peakpick)

    @property
    def centroid_spectrum(self):
        """Return whether centroid spectra are present."""
        return self._centroid_spectrum

    @centroid_spectrum.setter
    def centroid_spectrum(self, centroid_spectrum):
        """Set whether centroid spectra are present."""
        self._centroid_spectrum = centroid_spectrum
        self._set("centroid_spectrum", centroid_spectrum)

    @property
    def profile_spectrum(self):
        """Return whether profile spectra are present."""
        return self._profile_spectrum

    @profile_spectrum.setter
    def profile_spectrum(self, profile_spectrum):
        """Set whether profile spectra are present."""
        self._profile_spectrum = profile_spectrum
        self._set("profile_spectrum", profile_spectrum)

    @property
    def scan_direction(self):
        """Return the scan direction."""
        return self._scan_direction

    @scan_direction.setter
    def scan_direction(self, scan_direction):
        """Set the scan direction."""
        self._scan_direction = scan_direction
        self._set("scan_direction", scan_direction)

    @property
    def line_scan_direction(self):
        """Return the line scan direction."""
        return self._line_scan_direction

    @line_scan_direction.setter
    def line_scan_direction(self, line_scan_direction):
        """Set the line scan direction."""
        self._line_scan_direction = line_scan_direction
        self._set("line_scan_direction", line_scan_direction)

    @property
    def scan_pattern(self):
        """Return the scan pattern."""
        return self._scan_pattern

    @scan_pattern.setter
    def scan_pattern(self, scan_pattern):
        """Set the scan pattern."""
        self._scan_pattern = scan_pattern
        self._set("scan_pattern", scan_pattern)

    @property
    def scan_type(self):
        """Return the scan type."""
        return self._scan_type

    @scan_type.setter
    def scan_type(self, scan_type):
        """Set the scan type."""
        self._scan_type = scan_type
        self._set("scan_type", scan_type)

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
        """Return the meta index dictionary."""
        return self._meta_index

    @meta_index.setter
    def meta_index(self, meta_index_data):
        if meta_index_data is not None:
            if not isinstance(meta_index_data, dict):
                logger.error("meta_index must be a dict")
                raise TypeError("meta_index must be a dict")
            self._meta_index = meta_index_data

    def update_meta(self):
        """Update the meta index with new entries."""
        for meta_key in self._meta:
            self._meta[meta_key] = getattr(self, meta_key, None)


class MSIMetaData(MetaDataBase):
    """
    Metadata class for MSI image-matrix data.

    Extends MetaDataBase with two commonly used fields:
    - `mask`: the spatial mask of the sample (2D array or similar)
    - `need_base_mask`: whether to generate a base mask from intensity

    All fields are synchronized into the internal `_meta` dict via property setters.
    """

    def __init__(
        self,
        mask=None,
        need_base_mask: bool = False,
        name: str = "default",
        version: float = 1.0,
        storage_mode: str = "split",
        max_count_of_pixels_x: Optional[int] = None,
        max_count_of_pixels_y: Optional[int] = None,
        pixel_size_x: Optional[float] = None,
        pixel_size_y: Optional[float] = None,
        mz_num: Optional[int] = None,
    ):
        """
        Initialize MSI metadata instance.

        Parameters:
            name (str): Dataset name.
            version (float): Metadata version, must be > 0.
            storage_mode (str): Storage mode, e.g., 'split' or 'merge'.
            max_count_of_pixels_x (int): Pixel count along the X axis.
            max_count_of_pixels_y (int): Pixel count along the Y axis.
            pixel_size_x (float): Pixel size along the X axis.
            pixel_size_y (float): Pixel size along the Y axis.
            mask (Any): Spatial 2D mask used to define image dimensions.
            need_base_mask (bool): Whether base mask generation is required.

        Returns:
            None

        Raises:
            AssertionError: If `version` is not positive (checked in base class).
        """
        super().__init__(
            name=name,
            version=version,
            storage_mode=storage_mode,
            max_count_of_pixels_x=max_count_of_pixels_x,
            max_count_of_pixels_y=max_count_of_pixels_y,
            pixel_size_x=pixel_size_x,
            pixel_size_y=pixel_size_y,
            mask=mask,
        )

        self._need_base_mask = None
        self._mz_num = None

        # Set actual values through properties to sync into `_meta`
        self.need_base_mask = need_base_mask
        self.mz_num = mz_num

    @property
    def need_base_mask(self):
        """
        Get the flag indicating whether to generate a base mask.

        Returns:
            bool|None: The current flag value.
        """
        return self._need_base_mask

    @need_base_mask.setter
    def need_base_mask(self, need_base_mask: bool):
        """
        Set the flag indicating whether to generate a base mask.

        Parameters:
            need_base_mask (bool): True to enable base mask generation.

        Returns:
            None
        """
        if need_base_mask is not None:
            self._need_base_mask = bool(need_base_mask)
            self._set("need_base_mask", self._need_base_mask)

    @property
    def mz_num(self):
        """
        Get the number of m/z values.

        Returns:
            int|None: The current number of m/z values.
        """
        return self._mz_num

    @mz_num.setter
    def mz_num(self, mz_num: Optional[int]):
        """
        Set the number of m/z values.

        Parameters:
            mz_num (int): The number of m/z values.

        Returns:
            None
        """
        if mz_num is not None:
            self._mz_num = int(mz_num)
            self._set("mz_num", self._mz_num)


class ImzMlMetaData(MetaDataBase):
    """ImzML metadata wrapper that loads and caches frequently used fields."""

    def __init__(
        self,
        name="ImzML",
        version=1.0,
        storage_mode="split",
        filepath: Optional[str] = None,
        absolute_position_offset_x=None,
        absolute_position_offset_y=None,
        centroid_spectrum=None,
        profile_spectrum=None,
        ms1_spectrum=None,
        msn_spectrum=None,
        instrument_model=None,
        spectrum_count_num=None,
        min_pixel_x=None,
        min_pixel_y=None,
        mask=None,
        pixel_size_x=None,
        pixel_size_y=None,
        max_count_of_pixels_x=None,
        max_count_of_pixels_y=None,
    ):
        """Initialize the metadata object with file path."""
        super().__init__(
            name=name,
            version=version,
            storage_mode=storage_mode,
            mask=mask,
            pixel_size_x=pixel_size_x,
            pixel_size_y=pixel_size_y,
            max_count_of_pixels_x=max_count_of_pixels_x,
            max_count_of_pixels_y=max_count_of_pixels_y,
        )

        self._filepath = None
        self._spectrum_count_num = None
        self._absolute_position_offset_x = None
        self._absolute_position_offset_y = None
        self._instrument_model = None
        self._ms1_spectrum = None
        self._msn_spectrum = None
        self._min_pixel_x = 99999
        self._min_pixel_y = 99999

        # Set actual value through property
        self.filepath = filepath if filepath is not None else None

        self.absolute_position_offset_x = absolute_position_offset_x
        self.absolute_position_offset_y = absolute_position_offset_y
        self.centroid_spectrum = centroid_spectrum
        self.profile_spectrum = profile_spectrum
        self.ms1_spectrum = ms1_spectrum
        self.msn_spectrum = msn_spectrum
        self.instrument_model = instrument_model
        self.spectrum_count_num = spectrum_count_num
        self.min_pixel_x = min_pixel_x
        self.min_pixel_y = min_pixel_y

    @property
    def filepath(self):
        """Return the associated imzML file path."""
        return self._filepath

    @filepath.setter
    def filepath(self, filepath: Optional[str]):
        """Set the associated imzML file path."""
        if filepath is not None:
            if not os.path.exists(filepath):
                logger.error(f"ImzML file not found at {filepath}")
                raise FileNotFoundError(f"ImzML file not found at {filepath}")
            self._filepath = filepath

    @property
    def spectrum_count_num(self):
        """Return the number of spectra."""
        return self._spectrum_count_num

    @spectrum_count_num.setter
    def spectrum_count_num(self, spectrum_count_num: Optional[int]):
        """Persist the number of spectra and sync it to the base storage."""
        if spectrum_count_num is not None:
            self._spectrum_count_num = spectrum_count_num
            self._set("spectrum_count_num", spectrum_count_num)

    @property
    def absolute_position_offset_x(self):
        """Return the absolute position offset on the X axis."""
        return self._absolute_position_offset_x

    @absolute_position_offset_x.setter
    def absolute_position_offset_x(self, absolute_position_offset_x):
        """Set the absolute position offset on the X axis."""
        if absolute_position_offset_x is not None:
            self._absolute_position_offset_x = absolute_position_offset_x
            self._set("absolute_position_offset_x", absolute_position_offset_x)

    @property
    def absolute_position_offset_y(self):
        """Return the absolute position offset on the Y axis."""
        return self._absolute_position_offset_y

    @absolute_position_offset_y.setter
    def absolute_position_offset_y(self, absolute_position_offset_y):
        """Set the absolute position offset on the Y axis."""
        if absolute_position_offset_y is not None:
            self._absolute_position_offset_y = absolute_position_offset_y
            self._set("absolute_position_offset_y", absolute_position_offset_y)

    @property
    def instrument_model(self):
        """Return the mass spectrometer model."""
        return self._instrument_model

    @instrument_model.setter
    def instrument_model(self, instrument_model):
        """Set the mass spectrometer model."""
        if instrument_model is not None:
            self._instrument_model = instrument_model
            self._set("instrument_model", instrument_model)

    @property
    def ms1_spectrum(self):
        """Return whether MS1 spectra are present."""
        return self._ms1_spectrum

    @ms1_spectrum.setter
    def ms1_spectrum(self, ms1_spectrum):
        """Set whether MS1 spectra are present."""
        if ms1_spectrum is not None:
            self._ms1_spectrum = ms1_spectrum
            self._set("ms1_spectrum", ms1_spectrum)

    @property
    def msn_spectrum(self):
        """Return whether MSn spectra are present."""
        return self._msn_spectrum

    @msn_spectrum.setter
    def msn_spectrum(self, msn_spectrum):
        """Set whether MSn spectra are present."""
        if msn_spectrum is not None:
            self._msn_spectrum = msn_spectrum
            self._set("msn_spectrum", msn_spectrum)

    @property
    def min_pixel_x(self):
        """Return the minimum pixel X coordinate."""
        return self._min_pixel_x

    @min_pixel_x.setter
    def min_pixel_x(self, min_pixel_x):
        """Set the minimum pixel X coordinate."""
        if (
            min_pixel_x is not None
            and min_pixel_x >= 0
            and min_pixel_x <= self.max_count_of_pixels_x
        ):
            self._min_pixel_x = min_pixel_x
            self._set("min_pixel_x", min_pixel_x)

    @property
    def min_pixel_y(self):
        """Return the minimum pixel Y coordinate."""
        return self._min_pixel_y

    @min_pixel_y.setter
    def min_pixel_y(self, min_pixel_y):
        """Set the minimum pixel Y coordinate."""
        if (
            min_pixel_y is not None
            and min_pixel_y >= 0
            and min_pixel_y <= self.max_count_of_pixels_y
        ):
            self._min_pixel_y = min_pixel_y
            self._set("min_pixel_y", min_pixel_y)
