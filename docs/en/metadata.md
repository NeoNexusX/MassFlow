# MassFlow Metadata Module

This document describes the metadata subsystem in MassFlow, focusing on the three classes defined in `module/ms_meta_data.py`: `MetaDataBase`, `MSIMetaData`, and `ImzMlMetaData`. It covers their fields, properties, typical usage, and key considerations.

## Overview

- Design
  - Metadata (dataset info, instrument, coordinates, pixel sizes) is modeled as an object independent of the spectrum collection `MassSpectrumSet`. Data managers (e.g., `MSDataManagerImzML`) attach or update it while reading.
  - All attributes exposed via properties automatically synchronize to the internal dictionary `_meta`, enabling serialization and dict-like access.
- Capabilities
  - Record image size (number of pixels) and physical pixel size (µm).
  - Store occupancy mask (`mask`).
  - Cache common imzML metadata (spectrum count, instrument model, centroid/profile mode, etc.).
- Core classes
  - `MetaDataFileBase`: abstract base providing auto-sync and dict interfaces.
  - `MSIMetaData`: concrete subclass for MSI matrix-style metadata.
  - `ImzMlMetaData`: concrete subclass for imzML metadata, managing an `ImzMLParser`.

```mermaid
classDiagram
    MetaDataBase <|-- MSIMetaData
    MetaDataBase <|-- ImzMlMetaData

    class MetaDataBase {
        <<abstract base>>
    }

    class MSIMetaData {
        <<concrete class>>
    }

    class ImzMlMetaData {
        <<concrete class>>
    }
```

## Core Types

### MetaDataBase (abstract)

```python
class massflow.module.ms_meta_data.MetaDataBase(
    name: str = "default",
    version: float = 1.0,
    storage_mode: str = "split",
    max_count_of_pixels_x: int | None = None,
    max_count_of_pixels_y: int | None = None,
    pixel_size_x: float | None = None,
    pixel_size_y: float | None = None,
    mask: np.ndarray | None = None,
)
```

- Purpose
  - Provide common metadata fields and property wrappers; setters synchronize values into `_meta` via `self._set(key, value)`.
- Key fields/properties (partial)
  - `name`, `version`, `storage_mode`
  - `max_count_of_pixels_x`, `max_count_of_pixels_y` (pixel count)
  - `pixel_size_x`, `pixel_size_y` (µm)
  - `processed`, `peakpick`
  - `centroid_spectrum`, `profile_spectrum`
  - `mask` (2D array; shape must be `(max_count_of_pixels_y, max_count_of_pixels_x)`)
  - `meta_index` (CV mapping; used for extraction)
- Access
  - Dict-like: `__getitem__`, `keys()`, `items()`, `values()`, `get()`, `to_dict()`.

### MSIMetaData (matrix MSI)

```python
class massflow.module.ms_meta_data.MSIMetaData(
    mask=None,
    need_base_mask: bool = False,
    name: str = "default",
    version: float = 1.0,
    storage_mode: str = "split",
    max_count_of_pixels_x: int | None = None,
    max_count_of_pixels_y: int | None = None,
    pixel_size_x: float | None = None,
    pixel_size_y: float | None = None,
    mz_num: int | None = None,
)
```

- Extra fields
  - `need_base_mask`: whether to compute a base mask from intensities.
  - `mz_num`: global m/z count (optional, for statistics/control).
- Usage
  - Non-imzML sources or MSI converted to dense matrix format.
- Behavior
  - Inherits auto-sync: setting attributes writes into `_meta`.

### ImzMlMetaData (imzML)

```python
class massflow.module.ms_meta_data.ImzMlMetaData(
    name: str = "ImzML",
    version: float = 1.0,
    storage_mode: str = "split",
    filepath: str | None = None,
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
)
```

- Initialization
  - Provide `filepath: str` pointing to an existing `.imzML` file; the path is validated when set.
- Key fields
  - `filepath`, `spectrum_count_num`
  - `absolute_position_offset_x`, `absolute_position_offset_y`
  - `instrument_model`, `ms1_spectrum`, `msn_spectrum`
  - `min_pixel_x`, `min_pixel_y`
- Validation/sync
  - `filepath` must exist; otherwise a `FileNotFoundError` is raised.
  - `min_pixel_x/min_pixel_y` must satisfy `0 ≤ value ≤ max_count_of_pixels_*`.

## Usage Examples

### Scenario 1: MSI matrix

```python
>>> from massflow.module.ms_meta_data import MSIMetaData
>>> import numpy as np
>>>
>>> # 1. Create metadata object
>>> meta = MSIMetaData(name="test_dataset", mz_num=100, need_base_mask=True)
>>>
>>> # 2. Update metadata (e.g., set mask)
>>> my_mask = np.zeros((10, 10))
>>> meta.mask = my_mask
>>>
>>> # 3. Access via attribute and dict
>>> print(meta.name)
test_dataset
>>> print(meta["name"])
test_dataset
```

### Scenario 2: imzML

```python
>>> from massflow.module.mass_spectrum_set import MassSpectrumSet
>>> from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
>>>
>>> ms = MassSpectrumSet()
>>> ms_dm = MSDataManagerImzML(ms=ms, filepath="data/example.imzML")

>>> # Load spectra and extract metadata
>>> ms_dm.load_full_data_from_file()
>>>
>>> # Access metadata populated via CV index mapping
>>> print(f"X-axis pixels: {ms.meta.max_count_of_pixels_x}")
X-axis pixels: [X value from example.imzML]
>>> print(f"Y-axis pixels: {ms.meta.max_count_of_pixels_y}")
Y-axis pixels: [Y value from example.imzML]
>>> print(f"Min pixel x: {ms.meta.min_pixel_x}")
Min pixel x: 0
```

## Relation to Data Managers
- `MSDataManagerImzML` constructs and maintains `ImzMlMetaData` while reading, and attaches it to the underlying `MassSpectrumSet` (`ms.meta`) for subsequent analysis/plotting.
- For coordinate filtering or visualization, read `mask` and pixel size from metadata.

## Notes & Recommendations
- `version` must be positive.
- `filepath` must exist; for `ImzMlMetaData` this is validated when setting the `filepath` property.
- `mask` must be a 2D NumPy array; conventionally its shape should match `(max_count_of_pixels_y, max_count_of_pixels_x)`.
- Initialize `max_count_of_pixels_*` before setting `min_pixel_*` to pass bounds checks.
