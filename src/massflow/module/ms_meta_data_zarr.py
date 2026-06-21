"""Zarr-specific metadata adapter."""

from __future__ import annotations

from typing import Any

import numpy as np

from massflow.module.ms_meta_data import MetaDataBase, MetaField


class MetaDataZarr(MetaDataBase):
    """Metadata model persisted inside MassFlow spectra Zarr stores."""

    spectrum_count_num = MetaField()
    absolute_position_offset_x = MetaField()
    absolute_position_offset_y = MetaField()
    instrument_model = MetaField()
    ms1_spectrum = MetaField()
    msn_spectrum = MetaField()

    def __init__(
        self,
        name="zarr_data",
        version=1.0,
        absolute_position_offset_x=None,
        absolute_position_offset_y=None,
        centroid_spectrum=None,
        profile_spectrum=None,
        ms1_spectrum=None,
        msn_spectrum=None,
        instrument_model=None,
        spectrum_count_num=None,
        mask=None,
        pixel_size_x=None,
        pixel_size_y=None,
        max_count_of_pixels_x=None,
        max_count_of_pixels_y=None,
    ):
        super().__init__(
            name=name,
            version=version,
            mask=mask,
            pixel_size_x=pixel_size_x,
            pixel_size_y=pixel_size_y,
            max_count_of_pixels_x=max_count_of_pixels_x,
            max_count_of_pixels_y=max_count_of_pixels_y,
        )
        self.absolute_position_offset_x = absolute_position_offset_x
        self.absolute_position_offset_y = absolute_position_offset_y
        self.centroid_spectrum = centroid_spectrum
        self.profile_spectrum = profile_spectrum
        self.ms1_spectrum = ms1_spectrum
        self.msn_spectrum = msn_spectrum
        self.instrument_model = instrument_model
        self.spectrum_count_num = spectrum_count_num

    @classmethod
    def from_metadata(cls, source: MetaDataBase | None) -> "MetaDataZarr":
        """Create Zarr metadata from any MassFlow metadata object."""
        meta = cls()
        if source is None:
            return meta
        meta.copy(source)
        return meta

    @classmethod
    def from_zarr(cls, root: Any) -> "MetaDataZarr":
        """Load metadata from a MassFlow spectra Zarr root group."""
        group = root["metadata"] if "metadata" in root else root
        meta = cls()
        for key, value in dict(group.attrs).items():
            if hasattr(type(meta), key):
                setattr(meta, key, value)
            else:
                meta._set(key, value)
        if "mask" in group:
            meta.mask = np.asarray(group["mask"][:])
        return meta

    def to_zarr(self, root: Any) -> None:
        """Persist metadata into a MassFlow spectra Zarr root group."""
        group = root.create_group("metadata", overwrite=True)
        for key, value in self.items():
            if key == "mask" or value is None:
                continue
            json_value = _json_scalar(value)
            if json_value is not None:
                group.attrs[key] = json_value
        if self.mask is not None:
            group.create_array("mask", data=np.asarray(self.mask), overwrite=True)

    def copy(self, other: MetaDataBase):
        """Copy compatible metadata fields from another metadata object."""
        if other is None:
            return
        for key, value in other.items():
            if key == "mask":
                if value is not None:
                    self.mask = np.asarray(value).copy()
                continue
            if hasattr(type(self), key):
                setattr(self, key, value)
            else:
                self._set(key, value)


def _json_scalar(value: Any):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)):
        return value
    return None
