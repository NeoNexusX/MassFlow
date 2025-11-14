# MassFlow

English | [简体中文](README_CN.md)

MassFlow is a modular, high-performance framework for Mass Spectrometry Imaging (MSI) and Mass Spectrometry (MS) data. It provides lazy imzML reading, denoising/smoothing, baseline correction, normalization, peak alignment/picking, and MSI data management and export (.msi/.h5).

## Installation

Requirements: Python >= 3.9 (recommend Python 3.12)

```bash
# Clone the repository
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow

# Install dependencies
pip install -r requirements.txt
```

## Development and Contributing

### Contribution Workflow & Standards
- Contributing guide: `docs/en/contribution.md` (English) and `docs/zh/contribution.md` (Chinese)
- Naming conventions: `docs/en/naming-conventions.md` and `docs/zh/naming-conventions.md`
- Issue templates: `.github/ISSUE_TEMPLATE/feature.md`, `.github/ISSUE_TEMPLATE/bug.md`, `feature_en.md`, `bug_en.md`
- Commit messages: Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `refactor:`, `test:`). Example: `feat(data-manager): support split/merge write modes`
- Local checks: `ruff .`, `black .`, `isort .`, `pylint module/`
- PR checklist: consistent interfaces and naming, no performance/memory issues, assertions and error handling in place, docs/examples updated, CI passes, clear change description

### Recommended VSCode Extensions

Install these extensions to improve Python development and data handling:

- `ms-python.python`
- `ms-python.vscode-pylance`
- `ms-python.pylint`
- `charliermarsh.ruff`
- `ms-python.black-formatter`
- `ms-python.isort`
- `h5web.vscode-h5web`

Command line install (optional):

```bash
code --install-extension ms-python.python \
  && code --install-extension ms-python.vscode-pylance \
  && code --install-extension ms-python.pylint \
  && code --install-extension charliermarsh.ruff \
  && code --install-extension ms-python.black-formatter \
  && code --install-extension ms-python.isort \
  && code --install-extension h5web.vscode-h5web
```

## Quick Start

Open `example.ipynb` in Jupyter (recommended), or run the snippet below to verify data loading:

```python
from module.ms_module import MS
from module.ms_data_manager_imzml import MSDataManagerImzML

FILE_PATH = "data/your_file.imzML"
ms = MS()
with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:
    manager.load_full_data_from_file()
    manager.inspect_data()
    ms.plot_ms_mask()
```

Online docs: https://neonexusx.github.io/MassFlow/

Or minimal usage (read imzML and denoise one spectrum):


## Project Structure

```
MassFlow/
├── module/
│   ├── __init__.py
│   ├── meta_data.py
│   ├── ms_data_manager.py
│   ├── ms_data_manager_imzml.py
│   ├── ms_module.py
│   ├── msi_data_manager_base.py
│   ├── msi_data_manager_msi.py
│   ├── msi_data_manager_zys.py
│   └── msi_module.py
├── preprocess/
│   ├── ms_preprocess.py
│   ├── filter_helper.py
│   ├── baseline_correction_helper.py
│   ├── normalizer_helper.py
│   ├── peak_alignment.py
│   └── peak_pick_helper.py
├── docs/
│   ├── en/...
│   └── zh/...
├── .github/
│   └── ISSUE_TEMPLATE/...
├── tools/
│   └── plot.py
├── example.ipynb
├── example.py
├── logger.py
├── requirements.txt
├── LICENSE
└── README*.md
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

- Cardinal MSI: https://cardinalmsi.org/
- MATLAB Mass Spectrometry Preprocessing: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- PyOpenMS: https://pyopenms.readthedocs.io/

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

For questions and support, please open an issue on the GitHub repository.