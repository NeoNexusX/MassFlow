"""Zarr-specific metadata adapter."""

from __future__ import annotations

from typing import Any

import numpy as np

from massflow.module.ms_meta_data import MetaDataBase
from massflow.tools.logger import get_logger
logger = get_logger("massflow.module.meta_data_zarr")


class MetaDataZarr(MetaDataBase):
    """Metadata model persisted inside MassFlow spectra Zarr stores.

    Field schema lives on `MetaDataBase`; this subclass only supplies the
    zarr-specific default name and zarr (de)serialization.
    """

    def __init__(self, name="zarr_data", **kwargs):
        super().__init__(name=name, **kwargs)

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

    def copy(self, other: MetaDataBase, exclude: bool = True):
        """Copy compatible metadata fields from another metadata object."""
        if other is None:
            logger.warning("Source MS has no metadata to copy or meta is none.")
            return

        if exclude:
            exclude_keys = {
                "mask",
                "filepath",
                "spectrum_count_num",
                "continuous",
                "processed",
                "centroid_spectrum",
                "profile_spectrum",
            }
        else:
            exclude_keys = set()

        # Iterate through all attributes and copy if they exist in self and are not excluded
        for key, value in other.items():
            if key in exclude_keys:
                continue

            if hasattr(self, key):
                setattr(self, key, value)

        logger.info("Metadata copy completed.")


def _json_scalar(value: Any):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)):
        return value
    return None
