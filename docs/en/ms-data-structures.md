# MassFlow MS Data Structures

This document describes the core data structures used in MassFlow for mass spectrometry (MS) workflows, with focus on loading and organizing spectra from ImzML. It covers the following modules and their relationships:

- `module/spectrum.py`
- `module/mass_spectrum_set.py`
- `module/ms_data_manager.py`
- `module/ms_data_manager_imzml.py`

The intent is to present implementation-focused details about classes, attributes, indexing patterns, and the lazy-loading data flow.

## Overview

MassFlow separates domain models (data structures) from data managers (I/O and orchestration):

- Domain model
  - `Spectrum`: Represents a single spectrum with spatial coordinates.
  - `SpectrumImzML`: Specialized spectrum class for lazy-loading data from ImzML.
  - `MassSpectrumSet`: Collection of spectra with efficient coordinate-based indexing.
- Data manager
  - `MSDataManager`: Abstract base class defining common options and interface.
  - `MSDataManagerImzML`: Concrete manager that reads `.imzML` files and populates a `MassSpectrumSet`.

Metadata is handled separately (see `module/ms_meta_data.py`) and attached to managers/models as needed.


### Class Diagram (Inheritance)

The following diagram shows the inheritance relationships among the core classes described above. `SpectrumImzML` extends `Spectrum`, and `MSDataManagerImzML` extends `MSDataManager`. `MassSpectrumSet` is a standalone container class with no inheritance.

The `MSDataManager` and its concrete implementations (like `MSDataManagerImzML`) are responsible for populating the `MassSpectrumSet` container with spectrum data. The manager holds a reference to a `MassSpectrumSet` instance and loads spectra into it through the `load_full_data_from_file()` method. The `MassSpectrumSet` container holds multiple `Spectrum` instances (or its subclasses like `SpectrumImzML`), organizing them by spatial coordinates for efficient access.

```mermaid
classDiagram
    Spectrum <|-- SpectrumImzML
    MSDataManager <|-- MSDataManagerImzML
    MSDataManager --> MassSpectrumSet : populates
    MassSpectrumSet o-- Spectrum : contains
    
    class MassSpectrumSet {
        <<container>>
    }
    
    class Spectrum {
        <<base class>>
    }

    class SpectrumImzML {
        <<concrete class>>
    }

    class MSDataManager {
        <<abstract>>
    }

    class MSDataManagerImzML {
        <<concrete class>>
    }
```



## Core Types

### *Spectrum* Class

```python
class massflow.module.spectrum.Spectrum(mz_list, intensity, coordinate, sort_by_mz=True, shared_mz_list=None)
```

Represents a single mass spectrum bound to spatial coordinates. Key characteristics:

- Parameters
  - `mz_list` (*Optional[np.ndarray]*) — Array of m/z values. Can be `None` for lazy loading.
  - `intensity` (*Optional[np.ndarray]*) — Array of intensity values. Can be `None` for lazy loading.
  - `coordinate` (*Sequence[int]*) — Three integers `[x, y, z]` representing spatial coordinates.
  - `sort_by_mz` (*bool, optional*) — Whether the data is sorted by m/z. Defaults to `True`.
  - `shared_mz_list` (*Optional[np.ndarray]*) — Shared m/z axis for continuous data (managed by `MassSpectrumSet`).
- Properties
  - `mz_list: np.ndarray` — lazily resolved getter/setter; may be `None` until loaded or set.
  - `intensity: np.ndarray` — lazily resolved getter/setter; may be `None` until loaded or set.
- Utilities
  - `crop_range(x_range, sort_by_mz=True, mode="new")` — crop spectrum to a given m/z range.
  - `is_sorted()` — check whether the spectrum is sorted by m/z.
- Invariants and notes
  - `coordinate` must contain exactly 3 integers; indexing assumes `[x, y, z]` order.
  - `mz_list` and `intensity` must be same length when both are present.
  - Lazy loading is supported via `None` initialization, deferring data acquisition or swap-in from disk.

- Example

```python
>>> from massflow.tools.plot import plot_spectrum
>>> spectrum = ms[0] or ms.get_spectrum(0, 0, 0)
>>> mz_list = spectrum.mz_list
>>> print(mz_list)
output:
[16441.998  938.1308  2318.6423 ...  1174.1575  1333.138   1488.291 ]
>>> # Plot the spectrum using the standalone plotting utility
>>> plot_spectrum(base=spectrum)
```

- Example2

```python
if __name__ == "__main__":
    from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
    from massflow.module.mass_spectrum_set import MassSpectrumSet
    from massflow.tools.plot import plot_spectrum
    FILE_PATH = "data/example.imzML"
    ms = MassSpectrumSet()
    # Create MS collection and manager
    # Auto-created from parser
    with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:

        # Load data with lazy-loading placeholders
        manager.load_full_data_from_file()
        spectrum = ms[0]
        plot_spectrum(spectrum)

output:
```

![image-20251106110148417](https://s2.loli.net/2025/11/06/9u8I7wkNsvlT4jS.png)

### *SpectrumImzML* Class

```python
class massflow.module.spectrum_imzml.SpectrumImzML(coordinates, index=None, reader=None, ibd_path=None, mz_list=None, intensity=None, shared_mz_list=None, sort_by_mz=True)
```

Specialized mass spectrum class for ImzML format with lazy loading capabilities. This class extends `Spectrum` to provide efficient handling of ImzML (Imaging Mass Spectrometry Markup Language) format data by implementing lazy loading to minimize memory usage.

Note: use [MSDataManagerImzML](#msdatamanagerimzml-class) to automatically create and manage [SpectrumImzML](#spectrumimzml-class).

- Parameters
  - `coordinates` (*Sequence[int]*) — Three integers `[x, y, z]` representing spatial coordinates of the spectrum.
  - `index` (*int*) — Index of the spectrum within the ImzML file.
  - `reader` (*PortableSpectrumReader | None*) — Low-level reader used to stream spectra from the `.ibd` file.
  - `ibd_path` (*str | None*) — Path to the corresponding `.ibd` binary file.
  - `mz_list` (*Optional[np.ndarray]*) — Optional m/z axis; if `None`, loaded lazily on first access or shared via `shared_mz_list`.
  - `intensity` (*Optional[np.ndarray]*) — Optional intensity array; if `None`, loaded lazily on first access.
  - `shared_mz_list` (*Optional[np.ndarray]*) — Shared m/z axis injected by `MassSpectrumSet` for continuous data.
  - `sort_by_mz` (*bool, optional*) — Whether to keep data sorted by m/z. Defaults to `True`.
- Inherited Attributes
  - `coordinate` (*PixelCoordinates*) — 3D coordinates `[x, y, z]` of the spectrum.
  - `x`, `y`, `z` (*int*) — Individual coordinate components.
  - `sort_by_mz` (*bool*) — Flag indicating if data is sorted by m/z values.
- Properties
  - `mz_list: np.ndarray` — Lazily loaded array of m/z values; triggers data loading on first access.
  - `intensity: np.ndarray` — Lazily loaded array of intensity values; ensures data is loaded when accessed.
- Lazy loading behavior
  - On first access to `mz_list` or `intensity`, opens the `.ibd` file and uses `reader.read_spectrum_from_file(ibd_file, index)` to load both arrays at once.
  - `mz_list` is cached after the first load; subsequent accesses reuse cached data.
  - Data loading is deferred until first property access to minimize memory footprint.
- Invariants and notes
  - The actual m/z and intensity data are not loaded during initialization when `mz_list`/`intensity` are `None`.
  - Both `mz_list` and `intensity` are loaded together for efficiency when either is first accessed.
  - Inherits all visualization and manipulation methods from `Spectrum`.
- Example

```python
>>> # Create SpectrumImzML instance (data not yet loaded)
>>> spectrum = SpectrumImzML(parser, index=0, coordinates=[0, 0, 0])
>>> # First access triggers lazy loading
>>> mz_values = spectrum.mz_list
>>> print(mz_values[:5])
output:
[100.05    150.12    200.34    250.67    300.89]
```

### *MassSpectrumSet* Class

```python
class massflow.module.mass_spectrum_set.MassSpectrumSet()
```

Collection class for managing multiple mass spectra with coordinate-based indexing. This class serves as a container and manager for multiple `Spectrum` (or `SpectrumImzML`) instances, providing efficient storage, retrieval, and manipulation of mass spectrometry data organized by 3D spatial coordinates.

- Parameters
  - No parameters required for initialization.
- Attributes
  - `meta` (*Optional*) — Metadata object associated with the spectrum collection.
- Methods
  - `add_spectrum(spectrum)` — Add a mass spectrum to the collection with automatic coordinate indexing.
  - `get_spectrum(x, y, z=0)` — Retrieve a mass spectrum by its 3D coordinates.
  - `plot_ms_mask(save_path, figsize, dpi, origin, cmap)` — Plot the occupancy mask stored in metadata.
- Indexing patterns
  - `ms[index]` — Sequential access by integer index into `_queue`.
  - `ms[x, y]` — 2D coordinate access (z defaults to 0).
  - `ms[x, y, z]` — Full 3D coordinate access.
  - `ms[x, y, z] = spectrum` — Direct assignment with automatic indexing.
- Special methods
  - `__len__()` — Returns the total number of spectra in the collection.
  - `__iter__()` — Returns an iterator over all spectra in insertion order.
  - `__getitem__(key)` — Supports flexible indexing by integer or coordinate tuple.
  - `__setitem__(key, spectrum)` — Assigns spectrum to coordinates with automatic indexing updates.
- Invariants and notes
  - Maintains two internal data structures for efficient sequential and coordinate-based access.
  - Supports both 2D (x, y) and 3D (x, y, z) coordinate systems.
  - Coordinates are automatically extracted and indexed when adding spectra.
  - Assignment operations automatically update spectrum coordinates to match the key.
- Example

```python
>>> from massflow.module.mass_spectrum_set import MassSpectrumSet
>>> # Create MS collection
>>> ms = MassSpectrumSet()
>>> # Add spectra (user should not load by yourself, use data managers first)
>>> ms.add_spectrum(spectrum1)
>>> # Access by index
>>> spec = ms[0]
>>> # Access by coordinates
>>> spec = ms[10, 20, 0]
>>> # Direct assignment
>>> ms[5, 5, 0] = new_spectrum
>>> # Iteration
>>> for spectrum in ms:
...     print(spectrum.coordinates)
output:
[0, 0, 0]
[1, 0, 0]
[2, 0, 0]
```

- Example2

Since the coordinates may be non-contiguous, it is recommended to check the available coordinate range using a mask.

```python
>>> ms_md.load_full_data_from_file()
>>> ms_md.inspect_data()
>>> ms.plot_ms_mask()
```

<img src="https://s2.loli.net/2025/11/06/nOG2oBmLd715qup.png" alt="image-20251106120024072" style="zoom: 50%;" />

## Data Managers

### *MSDataManager* Class (Abstract)

```python
class massflow.module.ms_data_manager.MSDataManager(
    ms=None,
    target_mz_range=None,
    target_locs=None,
    filepath=None,
    temp_dir=None,
    max_threads: int = 8,
    mz_dtype=np.float64,
    intensity_dtype=np.float32,
)
```

Abstract base class that defines common configuration and the contract for loading mass spectrometry data into a `MassSpectrumSet` collection. This class provides the foundation for concrete data manager implementations that handle different file formats.

- Parameters
  - `ms` (*MassSpectrumSet | None*) — The target `MassSpectrumSet` instance to populate with spectra. If `None`, a new one is created.
  - `target_mz_range` (*Optional[Tuple[float, float]]*) — Inclusive range `(min_mz, max_mz)` to filter peaks. Can be `None` to load all m/z values.
  - `target_locs` (*Optional[List[Tuple[int, int] | Tuple[int, int, int]]]*) — A bounding region defined by two coordinates `[(x1, y1), (x2, y2)]` or 3D equivalents; used to limit loaded spectra. Can be `None` to load all locations.
  - `filepath` (*Optional[str]*) — Path to the input data file.
  - `temp_dir` (*Optional[str]*) — Directory for temporary swap `.imzML` files when `filepath` is not provided.
  - `max_threads` (*int*) — Maximum number of worker threads for parallel loading utilities.
  - `mz_dtype` — NumPy dtype for m/z values.
  - `intensity_dtype` — NumPy dtype for intensity values.
- Methods
  - `load_full_data_from_file()` — Abstract method; must be implemented by concrete managers to load data.
  - `inspect_data(inpect_num=10)` — Logs dataset information including count, sample spectra lengths and ranges.
- Validation
  - Ensures `target_locs` contains at least two coordinate points when provided.
  - Validates that `x1 < x2` and `y1 < y2` for bounding box coordinates.
  - Raises appropriate errors for invalid configuration.
- Invariants and notes
  - This is an abstract class and cannot be instantiated directly.
  - Concrete subclasses must implement `load_full_data_from_file()` method.
  - Provides common infrastructure for spatial and m/z range filtering.
  - Maintains a counter for tracking loading progress.
- Example

```python
# Cannot instantiate abstract class directly
# Use concrete implementations like MSDataManagerImzML
if __name__ == "__main__":
    from module.ms_data_manager_imzml import MSDataManagerImzML
    from massflow.module.mass_spectrum_set import MassSpectrumSet
    FILE_PATH = "data/example.imzML"
    ms = MassSpectrumSet()
    # Create MS collection and manager
    # Auto-created from parser
    with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:

        # Load data with lazy-loading placeholders
        manager.load_full_data_from_file()

        print(manager.current_spectrum_num)

output:
10000
```

### *MSDataManagerImzML* Class

```python
class massflow.module.ms_data_manager_imzml.MSDataManagerImzML(
    ms=None,
    target_locs=None,
    filepath=None,
    max_threads: int = 8,
    temp_dir=None,
    mz_dtype=np.float64,
    intensity_dtype=np.float32,
)
```

Concrete data manager for `.imzML` files. This class extends `MSDataManager` to handle ImzML format mass spectrometry imaging data, managing metadata initialization and lazy population of the `MassSpectrumSet` collection with `SpectrumImzML` placeholders.

- Parameters
  - `ms` (*MassSpectrumSet | None*) — The target `MassSpectrumSet` instance to populate with spectra.
  - `target_locs` (*Optional[List[Tuple[int, int] | Tuple[int, int, int]]]*) — Bounding region for spatial filtering. Can be `None`.
  - `filepath` (*Optional[str]*) — Path to the `.imzML` file to load.
  - `max_threads` (*int*) — Maximum worker threads for batch loading utilities.
  - `temp_dir` (*Optional[str]*) — Directory for temporary swap files when `filepath` is not provided.
  - `mz_dtype` — NumPy dtype for m/z values.
  - `intensity_dtype` — NumPy dtype for intensity values.
- Attributes
  - `parser` (*ImzMLParser | None*) — ImzML parser instance created from `filepath` for reading spectrum metadata.
  - `reader` (*PortableSpectrumReader | None*) — Binary reader bound to the `.ibd` file.
  - `meta` (*ImzMlMetaData*) — Metadata wrapper bound to the parser, caching commonly used fields like image dimensions and instrument information.
- Inherited Attributes
  - `ms` (*MassSpectrumSet*) — The target domain model to populate.
  - `target_mz_range`, `target_locs`, `filepath`, `current_spectrum_num` — From `MSDataManager`.
- Methods
  - `load_full_data_from_file()` — Implements the abstract method to load ImzML data with lazy-loading spectra.
  - `inspect_data(inpect_num=10)` — Inherited method for dataset inspection.
- Initialization logic
  - Validates and creates an `ImzMLParser` from `filepath`.
  - Initializes `ImzMlMetaData` and attaches it to the underlying `MassSpectrumSet`.
- Invariants and notes
  - Only supports `.imzML` file format; raises error for other extensions.
  - Uses lazy loading to minimize memory usage—actual spectrum data is not loaded until accessed.
  - Metadata is tightly coupled and provides essential information for downstream processing.
- Example

```python
from massflow.module.mass_spectrum_set import MassSpectrumSet
from module.ms_data_manager_imzml import MSDataManagerImzML

# Run examples when executing this file directly
if __name__ == "__main__":

    FILE_PATH = "data/example.imzML"
    ms = MassSpectrumSet()
    # Create MS collection and manager
    # Auto-created from parser
    with MSDataManagerImzML(ms=ms,target_locs=[(1, 1), (50, 50)],filepath=FILE_PATH) as manager:

        # Load data with lazy-loading placeholders
        manager.load_full_data_from_file()

        # Inspect loaded data
        manager.inspect_data(inpect_num=5)

        # Access spectrum (triggers lazy load)
        spectrum = ms[40, 13]
        print(spectrum.mz_list[:5])

inspect_data output:
INFO:     25-11-10 19:34 202 ms_data_manager - creating ms mask.
INFO:     25-11-10 19:34 102 ms_data_manager - MS meta data:
                                                 target_mz_range: None
                                                 target_locs: [(1, 1), (50, 50)]
                                                 filepath: data/example.imzML
                                                 current_spectrum_num: 1910
                                                 meta_name: ImzML
                                                 meta_version: 1.0
                                                 meta_storage_mode: split
                                                 meta_centroid_spectrum: None
                                                 meta_profile_spectrum: True
                                                 meta_max_count_of_pixels_x: 227
                                                 meta_max_count_of_pixels_y: 93
                                                 meta_pixel_size_x: 100.0
                                                 meta_pixel_size_y: 100.0
                                                 meta_absolute_position_offset_x: 0.0
                                                 meta_absolute_position_offset_y: 0.0
                                                 meta_min_pixel_x: 0
                                                 meta_min_pixel_y: 2
                                                 meta_mask: (93, 227)
                                               
INFO:     25-11-10 19:34 115 ms_data_manager - MS  information:
                                                 MS len: 74749
                                                 MS range: 400.0 - 1000.0
                                                 MS coord: (0, 38, 0)
                                                 max and min mz_list: 1000.0 - 400.0
                                                 max intensity: 16.307722091674805
                                               
                                                 MS len: 74749
                                                 MS range: 400.0 - 1000.0
                                                 MS coord: (0, 39, 0)
                                                 max and min mz_list: 1000.0 - 400.0
                                                 max intensity: 17.022132873535156
                                               
                                                 MS len: 74749
                                                 MS range: 400.0 - 1000.0
                                                 MS coord: (0, 40, 0)
                                                 max and min mz_list: 1000.0 - 400.0
                                                 max intensity: 16.470420837402344
                                               
                                                 MS len: 74749
                                                 MS range: 400.0 - 1000.0
                                                 MS coord: (0, 41, 0)
                                                 max and min mz_list: 1000.0 - 400.0
                                                 max intensity: 19.481334686279297
                                               
                                                 MS len: 74749
                                                 MS range: 400.0 - 1000.0
                                                 MS coord: (0, 42, 0)
                                                 max and min mz_list: 1000.0 - 400.0
                                                 max intensity: 18.155176162719727
                                               
                                               
[400.         400.00802697 400.01605394 400.02408091 400.03210788]
```

## Metadata Dependencies

While not part of the three focus modules, `ImzMlMetaData` (in `module/ms_meta_data.py`) is tightly coupled to `MSDataManagerImzML`:

- It stores and exposes commonly used fields (e.g., image dimensions, pixel sizes, instrument model) and computes `spectrum_count_num` from `parser.coordinates`.
- It can be initialized with either a parser or a file path and will extract metadata via `pyimzml`.

## Data Flow Summary

1. Create an empty `MassSpectrumSet` collection.
2. Initialize `MSDataManagerImzML` with the `MassSpectrumSet` instance, optional `target_locs`/`target_mz_range`, and `filepath`.
3. Call `load_full_data_from_file()` to populate the `MassSpectrumSet` with `SpectrumImzML` placeholders for spectra within the target region.
4. When accessing `mz_list` or `intensity` on a `SpectrumImzML` spectrum, lazy loading fetches the data from the ImzML file via the parser.
5. Metadata remains available through `meta` for consumers that need image dimensions, pixel sizes, instrument model, etc.

## Extensibility Notes

- New managers can subclass `MSDataManager` to support additional formats; adhere to the `load_full_data_from_file()` contract and, if practical, use lazy-loading spectra like `SpectrumImzML`.
- Consider using `target_mz_range` during load to prefilter peaks for large datasets, or leave it for downstream processing depending on performance needs.

## Glossary

- ImzML: An XML-based standard for MSI data storage.
- Spectrum: A pair of arrays `(mz_list, intensity)` associated with specific spatial coordinates.
- Lazy loading: Deferring actual data reading until first access to reduce memory and I/O overhead.