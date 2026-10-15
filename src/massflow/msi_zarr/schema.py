"""Schema definitions for the MassFlow MSI Zarr format."""

import os
from typing import Any, Literal

import numpy as np


ZARR_FORMAT_NAME = "massflow.msi_zarr"
ZARR_FORMAT_VERSION = "1.1"
ZARR_MODE_CONTINUOUS: Literal["continuous"] = "continuous"
ZARR_MODE_PROCESSED: Literal["processed"] = "processed"
ZARR_MODE_AUTO: Literal["auto"] = "auto"
ZarrStorageMode = Literal["continuous", "processed"]
ZarrWriteMode = Literal["auto", "continuous", "processed"]
ZarrRowAxis = Literal["pixel", "ion"]
ZARR_GROUP_ION_IMAGE = "ion_image"
ZARR_GROUP_SPECTRA = "spectra"


def normalize_zarr_path(path: str | os.PathLike[str]) -> str:
    """Return ``path`` as a string ending in ``.zarr``."""
    value = str(path)
    return value if value.endswith(".zarr") else f"{value}.zarr"


def read_array_int(array: Any, index: int) -> int:
    """Read one integer scalar without exposing Zarr's broad selection type."""
    return int(np.asarray(array[index]).item())
