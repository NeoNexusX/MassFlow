# MassFlow Naming Conventions

This guide uses the classes, modules, and storage layouts in the current source tree.

## Python names

- Modules, functions, parameters, and variables use `snake_case`: `ms_data_manager_imzml.py`, `load_head_data`, `write_flat_batch`, `mz_axis`.
- Classes use `PascalCase`, preserving established abbreviations: `MSDataManagerImzML`, `MSDataManagerZarr`, `IonImageDataManagerH5`, `MetaDataImzMl`.
- Constants use `ALL_CAPS`: `ZARR_MODE_CONTINUOUS`, `ZARR_FORMAT_VERSION`.
- Internal state and helpers use a leading underscore: `_reader`, `_queue`, `_require_reader`.
- Use existing public exports, for example `from massflow.data_manager import MSDataManagerImzML`.

## Classes and responsibilities

| Responsibility | Current classes |
| --- | --- |
| Individual spectrum | `Spectrum`, `SpectrumImzML`, `SpectrumZarr` |
| Spectrum collection | `MassSpectrumSet` |
| Pixel coordinates | `PixelCoordinates` |
| Spectrum managers | `MSDataManager`, `MSDataManagerImzML`, `MSDataManagerZarr` |
| Ion-image managers | `IonImageDataManager`, `IonImageDataManagerBinary`, `IonImageDataManagerH5`, `IonImageDataManagerZarr` |
| Metadata | `MetaDataBase`, `MetaDataImzMl`, `MetaDataZarr` |
| Preprocessing | `Preprocessor`, `PreprocessorAPI`, `FlatPreprocess` |

Source files are grouped by responsibility: `module/`, `data_manager/`, `preprocess/`, `msi_zarr/`, `segmentation/`, and `tools/`.

## Operations and parameters

Use verb-first operation names from the current APIs:

- `load_head_data()` loads metadata and lazy spectrum placeholders.
- `extract_metadata()` reads metadata into the manager's model.
- `flat_generator()` and `matrix_generator()` yield batches.
- `write_flat_batch()` writes flat batches; writer classes also expose `add_spectrum()`.
- `read_ion_image()` / `plot_ion_image()` read or display an image through visualization helpers.
- `get_ion_image()` is the ion-image manager operation.
- `finalize()` completes supported writers, and `close()` releases resources.

Use `mz`, `intensity`, `lengths`, `coordinates`, and `mz_axis` consistently. Verify signatures before copying examples between APIs: the low-level Zarr writer uses `add_spectrum(mz, intensity, coordinates)`.

## Metadata fields

Public fields omit the `meta_` prefix: `meta.pixel_size_x`, `meta.continuous`, `meta.profile_spectrum`. Descriptors and properties synchronize values into the `_meta` dictionary. Printed output adds `meta_` to labels; these printed labels are not Python property names.

Use `to_dict()`, `get()`, `items()`, `keys()`, and `values()` for metadata access. See [Metadata](./metadata.md) for the current classes and field semantics.

`MassSpectrumSet` exposes `add_spectrum()`, `update_spectrum_with_coord()`, `update_spectrum_with_index()`, iteration, indexing, and `len()`. Its queue is an implementation detail.

## Storage names

- Zarr group names are `axes`, `metadata`, `spectra`, `ion_image`, and `stats`. Array names include `intensity`, `offsets`, `coordinates`, and `mz`; see [Zarr schema](./zarr_schema.md).
- H5 ion-image datasets are `mz_axis`, `coordinates`, and `ion_images`, with root attributes for dimensions and counts.
- Encoding names are `continuous` and `processed`; writers may accept `auto` for inference. These names describe m/z organization.
- File access uses `mode="r"` or `mode="w"` where supported. Do not reuse encoding names as access modes.

## Examples and validation

Use meaningful variables such as `source_dm`, `processed_dm`, `mz_axis`, and `coordinates`. Keep file paths explicit in examples and obtain method defaults from the source signature.

Use exceptions at public input and I/O boundaries; assertions express internal invariants. Keep existing project naming when extending an API rather than introducing parallel aliases.
