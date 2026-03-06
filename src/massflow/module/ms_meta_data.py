# module/ms_meta_data.py
"""
Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""
import os
from typing import Optional
import numpy as np
from massflow.tools.logger import get_logger

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


class MetaField:
    """Descriptor that auto-syncs attribute values to the _meta dict.

    Parameters:
        default: Default value returned when the field is not set in _meta.
        allow_none: If True, None values are stored; otherwise None is ignored.
        converter: Optional callable applied to non-None values before storing.
    """

    def __init__(self, default=None, allow_none=False, converter=None):
        self.default = default
        self.allow_none = allow_none
        self.converter = converter

    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return obj._meta.get(self.name, self.default)

    def __set__(self, obj, value):
        if value is None and not self.allow_none:
            return
        if self.converter and value is not None:
            value = self.converter(value)
        obj._meta[self.name] = value


class MetaDataBase:
    """
    Abstract base class for MSI data models.

    Manages all common metadata fields, properties, and setters
    that are shared between different data representations.

    All metadata fields are automatically synchronized to the metadata
    dictionary via MetaField descriptors or property setters.
    """

    # --- Simple synced fields (skip None) ---
    name = MetaField()
    storage_mode = MetaField()
    max_count_of_pixels_x = MetaField(default=0)
    max_count_of_pixels_y = MetaField(default=0)
    pixel_size_x = MetaField()
    pixel_size_y = MetaField()

    # --- Fields that accept None assignment ---
    continuous = MetaField(allow_none=True)
    processed = MetaField(allow_none=True)
    peakpick = MetaField(allow_none=True)
    centroid_spectrum = MetaField(allow_none=True)
    profile_spectrum = MetaField(allow_none=True)
    scan_direction = MetaField(allow_none=True)
    line_scan_direction = MetaField(allow_none=True)
    scan_pattern = MetaField(allow_none=True)
    scan_type = MetaField(allow_none=True)

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
        self._version = None
        self._mask = None
        self._meta_index = None

        # Set values through descriptors / properties
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

    # --- Properties with validation (cannot be simple descriptors) ---

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
            self._version = version
            self._set("version", version)

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

    # --- Dict-like interface delegating to _meta ---

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

    need_base_mask = MetaField(converter=bool)
    mz_num = MetaField(converter=int)

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
            ValueError: If `version` is not positive (checked in base class).
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

        self.need_base_mask = need_base_mask
        self.mz_num = mz_num


class ImzMlMetaData(MetaDataBase):
    """ImzML metadata wrapper that loads and caches frequently used fields."""

    spectrum_count_num = MetaField()
    absolute_position_offset_x = MetaField()
    absolute_position_offset_y = MetaField()
    instrument_model = MetaField()
    ms1_spectrum = MetaField()
    msn_spectrum = MetaField()

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
        self._min_pixel_x = 99999
        self._min_pixel_y = 99999

        self.filepath = filepath
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

    # --- Properties with validation (cannot be simple descriptors) ---

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
