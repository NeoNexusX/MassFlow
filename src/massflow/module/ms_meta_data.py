# module/ms_meta_data.py
"""
Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""

from typing import Optional
import numpy as np

from massflow.tools.logger import get_logger

logger = get_logger("massflow.module.meta_data")

META_INDEX = {
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
        self.name = ""  # 显式声明，消除警告

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

    def __int__(self):
        """Allows int() to be called on the descriptor itself when accessed on the class."""
        return int(self.default) if self.default is not None else 0


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
        max_count_of_pixels_x: Optional[int] = None,
        max_count_of_pixels_y: Optional[int] = None,
        pixel_size_x: Optional[float] = None,
        pixel_size_y: Optional[float] = None,
        mask: Optional[np.ndarray] = None,
    ):
        self._meta = {}
        self._version = None
        self._mask = None

        # Set values through descriptors / properties
        self.name = name
        self.version = version
        self.max_count_of_pixels_x = max_count_of_pixels_x
        self.max_count_of_pixels_y = max_count_of_pixels_y
        self.pixel_size_x = pixel_size_x
        self.pixel_size_y = pixel_size_y
        self.mask = mask
        self.meta_index = META_INDEX

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

    def __str__(self):
        """Format metadata information as a string, skipping None or empty values."""
        lines = ["MSI MetaData:"]
        for attr, value in self.items():
            # Skip None or effectively empty values
            if value is None:
                continue

            shape = getattr(value, "shape", None)
            if shape is not None:
                if len(shape) > 0:
                    val_str = str(shape)
                else:
                    # Handle 0-d arrays or empty collections if they have shape
                    continue
            else:
                val_str = str(value)
            lines.append(f"  meta_{attr}: {val_str}")

        return "\r\n".join(lines)

    def copy(self, other: "MetaDataBase"):
        pass
