"""Schema definitions for the MassFlow MSI Zarr format."""

from typing import Literal


ZARR_FORMAT_NAME = "massflow.msi_zarr"
ZARR_FORMAT_VERSION = "1.0"
ZARR_MODE_CONTINUOUS: Literal["continuous"] = "continuous"
ZARR_MODE_PROCESSED: Literal["processed"] = "processed"
ZARR_MODE_AUTO: Literal["auto"] = "auto"
ZarrStorageMode = Literal["continuous", "processed"]
ZarrWriteMode = Literal["auto", "continuous", "processed"]
ZarrRowAxis = Literal["pixel", "ion"]
