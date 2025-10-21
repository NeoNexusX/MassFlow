# MassFlow Naming Conventions

This guide summarizes the existing naming style across the codebase (`module/msi_module.py`, `module/msi_data_manager.py`, `module/msi_data_manager_zys.py`, `example.py`, and `README*`) and standardizes it for Python code, HDF5 dataset names, and identifiers in examples/docs.

## General Rules
- Prefer clear, readable, and semantically meaningful names; avoid ambiguity and excessive abbreviations.
- Keep style consistent: use `snake_case` for files/functions/variables and `PascalCase` for classes.
- Preserve domain abbreviations in uppercase inside class names, e.g., `MSI`, `HDF5`, `ZYS`.
- Method names should start with verbs to express action semantics (e.g., `load_*`, `write_*`, `inspect_*`).

## Files and Directories
- Directory names are lowercase words; use underscores if needed: `module/`.
- Python files use `snake_case`: `msi_module.py`, `msi_data_manager.py`, `msi_data_manager_zys.py`.
- Example file: `example.py`; avoid placing business logic in examples.

## Modules and Packages
- Maintain public exports in `module/__init__.py` when appropriate.
- Prefer concise relative imports: `from .msi_module import MSI, MSIBaseModule`.

## Class Naming (PascalCase)
- Domain model classes: `MSI`.
- Manager/coordination classes: `MSIDataManager`, `MSIDataManagerZYS`.
- Primitive/slice classes: `MSIBaseModule`.
- Rules:
  - Base names derive from the domain concept or responsibility.
  - Keep domain abbreviations uppercase embedded in names, e.g., `MSIDataManager` (not `MsiDataManager`).

## Methods and Functions (snake_case)
- Verb-first names to reflect operations and targets: `load_full_data_from_file`, `load_meta_from_file`, `allocate_data_from_meta`, `add_msi_slice`, `get_msi_by_mz`, `plot_msi`, `inspect_data`, `write2local`.
- Helper/internal methods may use a leading underscore: `_inspect_hdf5_structure`, `_write_meta_data`, `_create_datasets`.
- Use clear parameter phrases for ranges/filters: `target_mz_range`, `display_threshold_percent`.

## Variable Naming (snake_case)
- Common variables: `mz`, `msroi`, `base_mask`, `coords`, `num_pixels`, `selected_channels`.
- Booleans represent conditions/states: `need_base_mask`.
- Constants use ALL_CAPS with underscores: `FILE_PATH`.
- Avoid single-letter names except in loops/math contexts (e.g., `i`).

## Attributes and Metadata
- Private metadata fields use `_meta_*` prefix: `_meta_name`, `_meta_version`, `_meta_mask`, `_meta_mz_num`, `_meta_storage_mode`, `_meta_need_base_mask`.
- Public properties use `meta_*` and are exposed via `@property`: `meta_name`, `meta_version`, `meta_mask`, `meta_mz_num`, `meta_storage_mode`, `meta_need_base_mask`.
- Data matrix attribute: public `data`, private `_data`.
- Use `update_metadata()` to automatically collect `_meta_*` fields into the metadata map.
- Queue access is unified: `add_msi_slice()`, `get_queue()`, `__len__()`, `__iter__()`, `__getitem__()`.

## Constants and Modes
- Use ALL_CAPS with underscores: `FILE_PATH`, `FILE_PATH_ZYS` (in examples).
- Mode strings are lowercase with clear semantics: `mode in {"merge", "split"}`; set `meta_storage_mode` before writing.

## Abbreviations and Proper Nouns
- `MSI` (Mass Spectrometry Imaging): keep uppercase in class/method names.
- `HDF5`: use the standard uppercase abbreviation in method/variable names (e.g., `_inspect_hdf5_structure`).
- `ZYS`: source data abbreviation can be used in derived class names `MSIDataManagerZYS` and relevant methods.
- Domain terms like `mz`, `msroi` follow dataset naming; do not translate or rename them.

## Private/Internal Conventions
- Internal methods/fields use a leading underscore: `_create_datasets`, `_write_meta_data`, `_metadata`, `_queue`.
- Access internal state via public methods/properties; avoid direct use of private fields.
- Assertions and error messages:
  - Use `assert` with explicit messages: `assert self.meta_mz_num > 0, "meta_mz_num must be greater than 0"`.
  - Parameter/file-type checks: `assert self.filepath.endswith('.mat'), "Error: filepath is not a .mat file."`.

## HDF5/Dataset Naming
- Group names: `mz_{mz_value:.4f}` (four decimal places).
- Datasets:
  - `mz`: the single m/z value.
  - `msroi`: 2D image matrix (gzip compression allowed).
  - Metadata datasets use `meta_*` names aligned with public properties.
- Metadata writing:
  - Public keys omit the leading underscore: `meta_version`, `meta_mz_num`.
  - Types: numeric for `version`/`num` (e.g., `np.float32`), strings via `h5py.string_dtype('utf-8')`.

## Parameters and Type Hints
- Add type hints where helpful, prioritizing key interfaces and return types: `def get_msi(self) -> MSI`, `h5_data_zys: Optional[h5py.File]`.
- Defaults should reflect common usage: `target_mz_range=None`, `storage_mode='split'`, `need_base_mask=False`.

## Example Consistency
- Example variables/instances follow the same style: `msi_dm`, `msi_dm_zys`, `target_mz_range=[101, 102]`.
- Entry point follows the standard idiom:
  ```python
  if __name__ == "__main__":
      ...
  ```

## Anti-Patterns (Avoid)
- Mixed casing styles: `MsiDataManager` (incorrect; use `MSIDataManager`).
- CamelCase or uppercase file names: `MsiModule.py` (incorrect; use `msi_module.py`).
- Non-semantic abbreviations: `tmp1`, `val2` (avoid; use semantic names).

## Practical Guidance
- Determine class/method responsibilities and data semantics first, then settle on names.
- Register domain abbreviations (e.g., `MSI`, `HDF5`, `ZYS`) and keep consistent uppercase usage.
- For metadata, adopt the `meta_*` prefix for properties and maintain mapping to private `_meta_*` fields.
- When writing/reading HDF5, strictly follow group/dataset naming and type conventions to prevent parsing inconsistencies.

---
To extend this guide for tests, logging, and exception categorization, add corresponding sections with code samples to this file as needed.