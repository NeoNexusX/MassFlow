# MassFlow

English | [简体中文](README_CN.md)

MassFlow is a modular, high-performance preprocessing and data management framework for Mass Spectrometry Imaging (MSI) data.

## Project structure

```
MassFlow/                        # uv workspace root
├── src/massflow/                # Core MSI framework (massflow package)
│   ├── data_manager/            # Data managers (MSDataManager, MSDataManagerImzML)
│   ├── module/                  # Core data models
│   │   ├── spectrum.py          # Spectrum base type
│   │   ├── spectrum_imzml.py    # ImzML spectrum with lazy loading
│   │   ├── mass_spectrum_set.py # Collection of spectra
│   │   └── ms_meta_data.py      # Metadata structures (ImzMlMetaData, …)
│   ├── preprocess/              # Preprocessing pipeline
│   │   ├── api.py               # Chainable task registration API (PreprocessorAPI)
│   │   ├── preprocessor.py      # Pipeline executor (Preprocessor)
│   │   └── flat_pre_fun.py      # Flat-array preprocessing (FlatPreprocess)
│   ├── r_preprocess/            # R/Cardinal adapter
│   └── tools/                   # Utility functions
├── mzorbit/                     # Agent layer workspace member (mzorbit package)
│   └── README.md                # mzorbit developer guide
├── tests/                       # massflow test suite
├── docs/                        # VitePress documentation (en + zh)
├── main.py                      # Quick smoke-test entry point
├── pyproject.toml               # massflow package + workspace root config
├── uv.lock                      # Locked dependency graph (commit this)
└── LICENSE
```

### Workspace layout

```
massflow  (src/massflow)  ←── core framework, no agent knowledge
    ▲
    │ depends on
mzorbit  (mzorbit/) ←── agent layer: tools, MCP server, annotation
```

Both packages share a single `uv.lock` managed by the workspace root.

---

## Online documentation

https://neonexusx.github.io/MassFlow/

---

## Quick start (development)

### 1. Install uv

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### 2. Clone and sync

```bash
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
uv sync
```

`uv sync` installs **both** `massflow` and `mzorbit` as editable packages in one step. No separate `pip install -e .` needed.

### 3. Optional extras

```bash
# R/Cardinal backend (requires a local R installation)
uv sync --extra r

# MCP server support for mzorbit
uv sync --extra mcp
```

### 4. Smoke test

```bash
uv run python main.py
```

---

## Running tests

```bash
# massflow core tests
uv run pytest tests/

# mzorbit tests
uv run --package mzorbit pytest mzorbit/tests/

# All tests in one shot
uv run pytest tests/ mzorbit/tests/
```

---

## mzorbit (agent layer)

`mzorbit` is the agent interface on top of MassFlow. It provides:

- CLI commands (`python -m mzorbit <cmd>`) for loading, QC, preprocessing, annotation, and end-to-end workflow
- An MCP server (`mzorbit-mcp-server`) that exposes these tools to MCP clients such as Claude Code and Cursor

See [mzorbit/README.md](mzorbit/README.md) for the full developer guide, CLI reference, and MCP setup.

---

## Development guidelines

- **Commit convention**: Conventional Commits — `feat:`, `fix:`, `docs:`, `refactor:`, `test:`
- **Naming conventions**: see `docs/en/naming-conventions.md` / `docs/zh/naming-conventions.md`
- **Contribution guide**: see `docs/en/contribution.md` / `docs/zh/contribution.md`
- **Recommended VS Code extensions**: Python, Pylance, Pylint, H5Web
- **R backend**: only needed when using `backend="cardinal"`. Set `R_HOME` if R is not on `PATH`.

---

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

---

## References

- [MATLAB Mass Spectrometry Preprocessing](https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html)
- [Cardinal MSI](https://cardinalmsi.org/) · [Guide](https://bioconductor.org/packages/devel/bioc/vignettes/Cardinal/inst/doc/Cardinal3-guide.html) · [GitHub](https://github.com/kuwisdelu/Cardinal)
- [PyOpenMS](https://pyopenms.readthedocs.io/) · [GitHub](https://github.com/OpenMS/OpenMS)

## Feedback

Open an issue on the [GitHub repository](https://github.com/NeoNexusX/MassFlow/issues).
