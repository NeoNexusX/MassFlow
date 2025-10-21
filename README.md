# MassFlow

English | [简体中文](README_CN.md)

MassFlow is a modular and high-performance preprocessing framework for Mass Spectrometry Imaging (MSI) data. It provides efficient data management, processing, and visualization capabilities for MSI research.

## Features

### ✅ Implemented Features

#### Data Import and Export
- **Supported formats**: `.h5`, `.msi`, `.mat` (MATLAB format)
- **Batch processing**: Import entire directories of MSI files
- **Flexible storage modes**: 
  - Split mode: Each m/z value saved as a separate file
  - Merge mode: All data combined into a single file
- **Data filtering**: Load specific m/z ranges to reduce memory usage
- **Metadata management**: Automatic metadata extraction and storage

#### Data Management
- **MSI Domain Model**: Clean API with metadata, slice queue, and optional data matrix
- **Memory-efficient loading**: Selective m/z range loading to minimize memory footprint
- **Data normalization**: Per-channel min-max normalization
- **Base mask generation**: Automatic generation of tissue region masks
- **Sparsity filtering**: Filter channels based on sparsity threshold

#### Visualization
- **MSI image plotting**: Display images for specified m/z ranges
- **Customizable display**: Adjustable thresholds, colormaps, and figure sizes
- **Batch visualization**: Plot multiple images within a specified range
- **Export support**: Save plots to file or display interactively

#### Specialized Processors
- **MSIDataManager**: General-purpose data manager for `.h5`/`.msi` files
- **MSIDataManagerZYS**: Specialized processor for MATLAB `.mat` format with custom preprocessing

### 📋 Planned Features (TODO)

#### Data Import and output
- .imzML format support
- .msi format support
- .mat format support only import may output

#### Baseline Correction

For example：TopHat Filter.

#### Denoising & Smoothing

#### Peak Detection

#### Peak Alignment

#### Normalization

- Total Ion Current (TIC)
- Median Normalization
- Reference Ion Normalization
- Binning（Resampling & Spectral Binning

## Installation

```bash
# Clone the repository
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow

# Install dependencies
pip install numpy h5py matplotlib pympler
```

## Quick Start

### Example 1: Load Complete MSI Data

```python
from module.msi_data_manager import MSIDataManager
from module.msi_module import MSI

# Create MSI instance
msi = MSI(name='example1', version=1.0, mask=None, need_base_mask=True)

# Initialize data manager
msi_dm = MSIDataManager(msi, filepath="./data/your_data.mat")

# Load all data
msi_dm.load_full_data_from_file()

# Inspect data
msi_dm.inspect_data()
```

### Example 2: Load Specific m/z Range

```python
# Load only m/z values between 101 and 102
msi2 = MSI(name='example2', version=1.0, mask=None, need_base_mask=True)
msi_dm2 = MSIDataManager(msi2, target_mz_range=[101, 102], filepath="./data/your_data.mat")
msi_dm2.load_full_data_from_file()

# Plot images in the loaded range
msi2.plot_msi()
```

### Example 3: Process MATLAB .mat Files

```python
from module.msi_data_manager_zys import MSIDataManagerZYS

# Create MSI instance for .mat file
msi3 = MSI(name='20250329_sample', version=1.0, mask=None, need_base_mask=True)

# Initialize ZYS data manager
msi_dm_zys = MSIDataManagerZYS(msi3, filepath="./data/sample.mat")

# Load and rebuild data
msi_dm_zys.load_data_from_zys_mat()
msi_dm_zys.rebuild_hdf5_file_from_zys()

# Plot specific m/z range
msi3.plot_msi(target_mz_range=[100, 150])
```

### Example 4: Export Data

```python
# Export as merged file
msi_dm_zys.write2local(mode='merge', output_fold='./data')

# Export as split files (one file per m/z)
msi_dm_zys.write2local(mode='split')
```

### Example 5: Query Specific m/z Values

```python
# Get MSI slices for a specific m/z value
msi_list = msi3.get_msi_by_mz(mz_value_min=88.1122, tol=1e-3)
image = msi_list[0].msroi.T
base_mask = msi_list[0].base_mask.T
```

## Project Structure

```
MassFlow/
├── module/
│   ├── __init__.py
│   ├── msi_module.py              # Core MSI domain model
│   ├── msi_data_manager.py         # General data manager
│   └── msi_data_manager_zys.py     # MATLAB .mat format processor
├── example.py                       # Usage examples
├── README.md
└── LICENSE                          # GNU GPL v3
```

## Core Components

### MSI Class
The main domain model that encapsulates:
- Metadata management
- MSI slice queue
- Data matrix storage
- Query and visualization methods

### MSIDataManager
General-purpose data manager supporting:
- `.h5` and `.msi` file formats
- Batch import from directories
- Split/merge export modes
- m/z range filtering

### MSIDataManagerZYS
Specialized manager for MATLAB `.mat` files with:
- Custom data structure handling
- Per-channel normalization
- Sparsity-based filtering
- HDF5 conversion capabilities

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## References

- [MSI Data Processing Workflow](https://pleinelune-r.github.io/2025/08/05/MSI%E6%95%B0%E6%8D%AE%E5%A4%84%E7%90%86%E6%B5%81%E7%A8%8B/)
- [Cardinal MSI](https://cardinalmsi.org/)
- [Cardinal GitHub Repository](https://github.com/kuwisdelu/Cardinal/tree/devel/R)
- [MATLAB: Preprocessing Raw Mass Spectrometry Data](https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html)
- [Mass Spectrometry Imaging Preprocessing Review](https://www.sciencedirect.com/science/article/pii/S0169743921001015)
- [PyOpenMS Documentation](https://pyopenms.readthedocs.io/en/latest/user_guide/background.html#why-use-openms)
- [Recent Advances in MSI Analysis](https://pubs.acs.org/doi/10.1021/jasms.4c00314)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

For questions and support, please open an issue on the GitHub repository.