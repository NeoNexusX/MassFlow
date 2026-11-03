# MassFlow

English | [简体中文](README_CN.md)

MassFlow is a modular preprocessing and data management framework for Mass Spectrometry Imaging (MSI) and Mass Spectrometry (MS) data.

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
├── logs/                    # Logs created by explicit initialization
├── src/
│   └── massflow/            # Core source code
│       ├── data_manager/         # imzML / spectra Zarr and Binary / H5 / Zarr ion-image managers
│       ├── module/               # Core data models
│       │   ├── spectrum.py       # Spectrum base type
│       │   ├── spectrum_imzml.py # ImzML spectrum with lazy loading (SpectrumImzML)
│       │   ├── mass_spectrum_set.py  # Collection of spectra (MassSpectrumSet)
│       │   ├── ms_meta_data.py       # Shared metadata (MetaDataBase)
│       │   └── ...
│       ├── preprocess/              # Preprocessing entry points and helper functions
│       │   ├── api.py               # Chainable task registration API (PreprocessorAPI)
│       │   ├── preprocessor.py      # Pipeline executor (Preprocessor)
│       │   ├── flat_pre_fun.py      # Batch processing utilities (FlatPreprocess)
│       │   └── ...
│       ├── r_preprocess/       # R/Cardinal adapter
│       ├── msi_zarr/           # MSI Zarr schema, reader, and writer
│       └── tools/              # Utility functions
├── tests/                   # Test cases
├── LICENSE
├── main.py
├── package.json
├── pyproject.toml           # Project configuration and dependencies
└── README.md
```

## Development & Contribution

### Quick Dev Start

It is recommended to run the following code (`main.py`) directly to verify data loading:

```bash
# First, close conda:
conda deactivate

# Install uv:
# Please refer to: https://docs.astral.sh/uv/getting-started/installation/

# For example:
# Linux && macOS:
curl -LsSf https://astral.sh/uv/install.sh | sh
# Windows:
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# Install dependencies:
uv sync

# Install optional R backend dependencies (requires local R):
uv sync --extra r

# Run the code
uv run python main.py
```

- `uv sync` installs the project itself by default, in editable mode (equivalent to `pip install -e .`); use `--no-editable` to disable editable installation.
- The default installation (`uv sync`) does not include R/Cardinal backend dependencies.
- Development dependencies are declared in the `dependency-groups.dev` section of `pyproject.toml`; install them with `uv sync --group dev`.
- To use the Cardinal backend (`backend="cardinal"`), install the `r` optional dependency and make sure R can be found on the system (on `PATH`, or via `R_HOME`).

- Please read the Contribution Guide first: `docs/en/contribution.md` and `docs/zh/contribution.md`
- Naming Conventions: `docs/en/naming-conventions.md` and `docs/zh/naming-conventions.md`
- Commit Convention: Conventional Commits (e.g., `feat:`, `fix:`, `docs:`, `refactor:`, `test:`)
- Recommended Extensions: Python, Pylance, Pylint, H5Web

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## References

- MATLAB Mass Spectrometry Preprocessing: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- Cardinal MSI: https://cardinalmsi.org/
  - Cardinal Guide: https://bioconductor.org/packages/devel/bioc/vignettes/Cardinal/inst/doc/Cardinal3-guide.html
  - Cardinal Github: https://github.com/kuwisdelu/Cardinal
  - Matter Github: https://github.com/kuwisdelu/matter
- PyOpenMS: https://pyopenms.readthedocs.io/
  - Github: https://github.com/OpenMS/OpenMS
  - Docs: https://github.com/OpenMS/OpenMS-docs
  - Contributed libraries: https://github.com/OpenMS/contrib
  - Flash APP: https://github.com/OpenMS/FLASHApp

## Feedback

For support or to report issues, please open an issue on the GitHub repository.
