# MassFlow

English | [简体中文](README_CN.md)

MassFlow is a modular preprocessing and data management framework for Mass Spectrometry Imaging (MSI) and Mass Spectrometry (MS) data. Currently supports:

## Installation

Requirements: Python >= 3.12 (Version 3.12 is recommended)

```bash
# Clone the repository
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
```

## Quick Start

It is recommended to open `example.ipynb` in Jupyter, or run the following code snippet directly to verify data loading:

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

Online Documentation: https://neonexusx.github.io/MassFlow/

## Project Structure

```
MassFlow/
├── .github/                 # GitHub configuration (Issue templates, Workflows)
├── .vscode/                 # VSCode configuration
├── data/                    # Example data
├── docs/                    # Documentation (VitePress)
│   ├── en/
│   └── zh/
├── logs/                    # Runtime logs
├── src/
│   └── massflow/            # Core source code
│       ├── module/          # Data models and managers
│       │   ├── ms_module.py # MS/ImzML base types
│       │   └── ...
│       ├── preprocess/      # Preprocessing algorithms
│       │   ├── ms_preprocess.py
│       │   └── ...
│       └── tools/           # Utility functions
├── tests/                   # Test cases
├── LICENSE
├── main.py
├── package.json
├── pyproject.toml           # Project configuration and dependencies
└── README.md
```

## Development & Contribution

- Contribution Guide: `docs/en/contribution.md` and `docs/zh/contribution.md`
- Naming Conventions: `docs/en/naming-conventions.md` and `docs/zh/naming-conventions.md`
- Issue Templates: `.github/ISSUE_TEMPLATE/feature_en.md`, `bug_en.md`, `feature.md`, `bug.md`
- Local Checks: `ruff .`, `black .`, `isort .`, `pylint module/`
- Commit Convention: Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `refactor:`, `test:`)
- Recommended Extensions: Python, Pylance, Ruff, Black, isort, Pylint, Markdownlint, GitLens, H5Web

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## References

- Cardinal MSI: https://cardinalmsi.org/
- MATLAB Mass Spectrometry Preprocessing: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- PyOpenMS: https://pyopenms.readthedocs.io/

## Feedback

For support or to report issues, please open an issue on the GitHub repository.