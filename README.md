# MassFlow

English | [简体中文](README_CN.md)

MassFlow is a modular and high-performance preprocessing framework for Mass Spectrometry Imaging (MSI) data. It provides efficient data management, processing, and visualization capabilities for MSI research.

## Installation

```bash
# Clone the repository
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow

# Install dependencies
pip install -r requirements.txt
```

## Development and Contributing

### Contribution Workflow & Standards
- Contributing guide: see `docs/CONTRIBUTING_EN.md` (English) and `docs/CONTRIBUTING.md` (Chinese)
- Issue templates: `.github/ISSUE_TEMPLATE/feature.md`, `.github/ISSUE_TEMPLATE/bug.md` and English templates `feature_en.md`, `bug_en.md` (available in GitHub “New issue” page)
- Naming conventions: see `docs/NAMING_CONVENTIONS_EN.md` and `docs/NAMING_CONVENTIONS.md`
- Commit messages: follow Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `refactor:`, `test:`). Example: `feat(data-manager): support split/merge write modes`
- Local checks (sample commands): `ruff .`, `black .`, `isort .`, `pylint module/`
- PR checklist highlights: consistent interfaces and naming, no obvious performance/memory issues, assertions and error handling in place, docs/examples/tests updated, CI passes, clear change description

### Recommended Trae/VSCode Extensions (Python & Code Quality)

To improve code quality and developer experience, install the following extensions. The workspace file `.vscode/extensions.json` also recommends them automatically.

- `ms-python.python` — Python support
- `ms-python.vscode-pylance` — IntelliSense and type analysis
- `ms-python.pylint` — Static analysis (uses repo `.pylintrc`)
- `charliermarsh.ruff` — Fast linting and style enforcement
- `ms-python.black-formatter` — Code formatting (Black)
- `ms-python.isort` — Import sorting
- `h5web.vscode-h5web` — HDF5 visualization

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

## Project Structure

```
MassFlow/
├── module/                           # Source code package (Python modules for the current project)
│   ├── __init__.py                   # Package initialization (exports main classes/interfaces)
│   ├── msi_module.py                 # Core MSI domain models (MSI class, slices, etc.)
│   ├── msi_data_manager.py           # Generic data manager (.h5 / .msi files, etc.)
│   ├── msi_data_manager_msi.py       # .msi/.h5 specific implementation (if exists)
│   ├── msi_data_manager_zys.py       # MATLAB .mat specific manager
│   └── msi_preprocess.py             # Preprocessing utilities and workflows
├── example.py                         # Usage examples / quick start script
├── data/                              # Sample data (can be added to .gitignore)
├── docs/                              # Documentation (contribution guidelines, naming conventions, etc.)
│   ├── CONTRIBUTING.md
│   ├── CONTRIBUTING_EN.md
│   ├── NAMING_CONVENTIONS.md
│   └── NAMING_CONVENTIONS_EN.md
├── tests/                             # Unit tests
│   └── test_imports.py
├── .github/                           # GitHub configuration (Issue templates, CI, etc.)
│   └── ISSUE_TEMPLATE/
├── .vscode/                           # VSCode workspace suggestions (extensions.json, etc.)
├── README.md                          # English documentation
├── README_CN.md                       # Chinese documentation (current file)
├── requirements.txt                   # Dependencies list
├── .pylintrc                          # pylint configuration
├── LICENSE                            # License (GNU GPL v3)
└── setup.cfg / pyproject.toml         # Optional packaging/build configuration
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