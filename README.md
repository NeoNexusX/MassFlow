# MassFlow

English | [简体中文](README_CN.md)

MassFlow is a modular preprocessing and data management framework for Mass Spectrometry Imaging (MSI) data.

## Get the code

```bash
# Clone the repository
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
```

## Online Documentation

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
│       ├── module/               # Core data models and data managers
│       │   ├── spectrum.py       # Spectrum base type
│       │   ├── spectrum_imzml.py # ImzML spectrum with lazy loading
│       │   ├── mass_spectrum_set.py  # Collection of spectra
│       │   ├── ms_meta_data.py       # Metadata structures (ImzMlMetaData, etc.)
│       │   └── ...
│       ├── preprocess/              # Preprocessing entry points and helpers
│       │   ├── spectrum_pre_fun.py  # Spectrum-level preprocessing API (SpectrumPreprocess)
│       │   ├── dm_pre_fun.py        # Data-manager-level preprocessing API (Preprocess)
│       │   ├── batch_pre_fun.py     # Batch processing utilities (BatchPreprocess)
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

### Quick Dev Start

It is recommended to  run the following code (main.py) directly to verify data loading:

```bash
# Close conda first:
conda deactivate

# Install uv:
please follow :https://docs.astral.sh/uv/getting-started/installation/

# For Example:
# Linux && macOS:
curl -LsSf https://astral.sh/uv/install.sh | sh 
# Windows:
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# Install dependencies:
uv sync 
uv pip install -e .

#run the code
uv run python main.py
```

- please read Contribution Guide: `docs/en/contribution.md` and `docs/zh/contribution.md`  first 
- Naming Conventions: `docs/en/naming-conventions.md` and `docs/zh/naming-conventions.md`
- Commit Convention: Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `refactor:`, `test:`)
- Recommended Extensions: Python, Pylance, Pylint, H5Web

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## References

- MATLAB Mass Spectrometry Preprocessing: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- Cardinal MSI: https://cardinalmsi.org/
  - Cardinal  Guide Book;https://bioconductor.org/packages/devel/bioc/vignettes/Cardinal/inst/doc/Cardinal3-guide.html
  - Cardinal  Github：https://github.com/kuwisdelu/Cardinal
  - Matter Github：https://github.com/kuwisdelu/matter
- PyOpenMS: https://pyopenms.readthedocs.io/
  - github：https://github.com/OpenMS/OpenMS
  - docs：https://github.com/OpenMS/OpenMS-docs
  - contributing libraries ：https://github.com/OpenMS/contrib
  - Flash APP ：https://github.com/OpenMS/FLASHApp

## Feedback

For support or to report issues, please open an issue on the GitHub repository.