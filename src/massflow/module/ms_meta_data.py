# module/ms_meta_data.py
"""
Author: MassFlow Development Team Bionet/dre,NeoNeuxs
License: See LICENSE file in project root
"""

from typing import Any, Iterator, Optional
import numpy as np
from pyimzml.metadata import ParamGroup

from massflow.tools.logger import get_logger

logger = get_logger("massflow.module.meta_data")

# accession -> field name (leaf CV terms, value copied as-is)
# accession -> (field name, canonical short key) (category CV terms, whose
#   presence -- not value -- selects one of several mutually exclusive child
#   terms; the short key matches the vocabulary used by stream_imzml_writer's
#   obo_codes so a round-tripped value can be written back without translation)
META_INDEX = {
    "fileDescription": {
        "IMS:1000030": "continuous",
        "IMS:1000031": "processed",
    },
    "referenceableParamGroupList": {
        "MS:1000127": "centroid_spectrum",
        "MS:1000128": "profile_spectrum",
        "MS:1000579": "ms1_spectrum",
    },
    "scanSettingsList": {
        "IMS:1000042": "max_count_of_pixels_x",
        "IMS:1000043": "max_count_of_pixels_y",
        "IMS:1000046": "pixel_size_x",
        "IMS:1000047": "pixel_size_y",
        "IMS:1000053": "absolute_position_offset_x",
        "IMS:1000054": "absolute_position_offset_y",
        # scan pattern children
        "IMS:1000410": ("scan_pattern", "meandering"),
        "IMS:1000411": ("scan_pattern", "one_way"),
        "IMS:1000412": ("scan_pattern", "random_access"),
        "IMS:1000413": ("scan_pattern", "flyback"),
        # scan type children
        "IMS:1000480": ("scan_type", "horizontal_line"),
        "IMS:1000481": ("scan_type", "vertical_line"),
        # scan direction children
        "IMS:1000400": ("scan_direction", "bottom_up"),
        "IMS:1000401": ("scan_direction", "top_down"),
        "IMS:1000402": ("scan_direction", "left_right"),
        "IMS:1000403": ("scan_direction", "right_left"),
        "IMS:1000404": ("scan_direction", "no_direction"),
        # line scan direction children
        "IMS:1000490": ("line_scan_direction", "line_right_left"),
        "IMS:1000491": ("line_scan_direction", "line_left_right"),
        "IMS:1000492": ("line_scan_direction", "line_bottom_up"),
        "IMS:1000493": ("line_scan_direction", "line_top_down"),
    },
    
}

# META_INDEX top-level keys use imzML element names; translate to the
# pyimzml.metadata.Metadata attribute that holds the parsed equivalent.
AREA_TO_PYIMZML = {
    "fileDescription": "file_description",
    "referenceableParamGroupList": "referenceable_param_groups",
    "scanSettingsList": "scan_settings",
}

# Virtual section grouping every populated field that no META_INDEX
# section covers (e.g. name, version, contact_name, instrument_model).
OTHER_SECTION = "other"


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
    contact_name = MetaField()
    max_count_of_pixels_x = MetaField(default=0)
    max_count_of_pixels_y = MetaField(default=0)
    pixel_size_x = MetaField()
    pixel_size_y = MetaField()
    absolute_position_offset_x = MetaField()
    absolute_position_offset_y = MetaField()
    spectrum_count_num = MetaField()
    instrument_model = MetaField()  # not auto-extracted; see populate_from_parser
    ms1_spectrum = MetaField()
    msn_spectrum = MetaField()  # not auto-extracted; see populate_from_parser

    # --- Fields that accept None assignment ---
    continuous = MetaField(allow_none=True)
    processed = MetaField(allow_none=True)
    centroid_spectrum = MetaField(allow_none=True)
    profile_spectrum = MetaField(allow_none=True)
    scan_direction = MetaField(allow_none=True)
    line_scan_direction = MetaField(allow_none=True)
    scan_pattern = MetaField(allow_none=True)
    scan_type = MetaField(allow_none=True)

    def __init__(
        self,
        name: str = "default",
        contact_name: Optional[str] = None,
        version: float = 1.0,
        max_count_of_pixels_x: Optional[int] = None,
        max_count_of_pixels_y: Optional[int] = None,
        pixel_size_x: Optional[float] = None,
        pixel_size_y: Optional[float] = None,
        absolute_position_offset_x: Optional[float] = None,
        absolute_position_offset_y: Optional[float] = None,
        instrument_model: Optional[str] = None,
        spectrum_count_num: Optional[int] = None,
        ms1_spectrum: Optional[bool] = None,
        msn_spectrum: Optional[bool] = None,
        mask: Optional[np.ndarray] = None,
    ):
        self._meta = {}
        self._version = None
        self._mask = None

        # Set values through descriptors / properties
        self.name = name
        self.contact_name = contact_name or "MassFlow"
        self.version = version
        self.max_count_of_pixels_x = max_count_of_pixels_x
        self.max_count_of_pixels_y = max_count_of_pixels_y
        self.pixel_size_x = pixel_size_x
        self.pixel_size_y = pixel_size_y
        self.absolute_position_offset_x = absolute_position_offset_x
        self.absolute_position_offset_y = absolute_position_offset_y
        self.instrument_model = instrument_model
        self.spectrum_count_num = spectrum_count_num
        self.ms1_spectrum = ms1_spectrum
        self.msn_spectrum = msn_spectrum
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

    # --- Section-based access following the META_INDEX segmentation ---

    def sections(self) -> list[str]:
        """Return the available section names, including OTHER_SECTION."""
        return [*self.meta_index, OTHER_SECTION]

    def get_section(self, section: str) -> dict[str, Any]:
        """Return the populated metadata fields belonging to one section.

        Parameters:
            section: A top-level META_INDEX key (e.g. "scanSettingsList"),
                or OTHER_SECTION for populated fields no section covers.

        Returns:
            Mapping of field name -> stored value. Fields that were never
            populated (or hold None) are omitted.

        Raises:
            KeyError: If `section` is not one of `sections()`.
        """
        if section in self.meta_index:
            result = {}
            for target in self.meta_index[section].values():
                field = target[0] if isinstance(target, tuple) else target
                value = self._meta.get(field)
                if field not in result and value is not None:
                    result[field] = value
            return result
        if section == OTHER_SECTION:
            indexed = {
                target[0] if isinstance(target, tuple) else target
                for accession_map in self.meta_index.values()
                for target in accession_map.values()
            }
            return {
                key: value
                for key, value in self._meta.items()
                if key not in indexed and value is not None
            }
        valid = ", ".join(self.sections())
        logger.error(f"Unknown metadata section {section!r}; expected one of: {valid}")
        raise KeyError(
            f"Unknown metadata section {section!r}; expected one of: {valid}"
        )

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

    def copy(self, other: "MetaDataBase", exclude: bool = True) -> None:
        """Copy compatible metadata fields from another metadata object."""
        if other is None:
            logger.warning("Source MS has no metadata to copy or meta is none.")
            return

        exclude_keys = {
            "mask",
            "filepath",
            "spectrum_count_num",
            "continuous",
            "processed",
            "centroid_spectrum",
            "profile_spectrum",
        } if exclude else set()

        # Iterate through all attributes and copy if they exist in self and are not excluded
        for key, value in other.items():
            if key in exclude_keys:
                continue

            if hasattr(self, key):
                setattr(self, key, value)

        logger.info("Metadata copy completed.")

    def populate_from_parser(self, parser: Any) -> None:
        """Batch-extract metadata from an ImzMLParser and fill fields defined by `meta_index`."""
        if parser is None:
            logger.error("Parser is not initialized. Please set parser or filepath first.")
            raise ValueError("Parser is not initialized. Please set parser or filepath first.")

        md = parser.metadata
        for area_name, accession_map in self.meta_index.items():
            area = getattr(md, AREA_TO_PYIMZML[area_name], None)
            for param_group in self._iter_param_groups(area):
                for _, accession, value, *_ in param_group.cv_params:
                    target = accession_map.get(accession)
                    if target is None:
                        continue
                    if isinstance(target, tuple):
                        field, canonical = target
                        setattr(self, field, canonical)
                    else:
                        setattr(self, target, value)

    @staticmethod
    def _iter_param_groups(area: Any) -> Iterator[ParamGroup]:
        """Yield each ParamGroup contained in one metadata area (single group or dict of groups)."""
        if area is None:
            return
        if isinstance(area, ParamGroup):
            yield area
            return
        if isinstance(area, dict):
            yield from (pg for pg in area.values() if isinstance(pg, ParamGroup))
