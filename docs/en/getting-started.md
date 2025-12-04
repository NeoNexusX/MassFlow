---
title: Getting Started
---

# Getting Started

MassFlow is a modular computational framework for Mass Spectrometry Imaging (MSI) and Mass Spectrometry (MS) data, providing efficient data reading, preprocessing, and management capabilities.
Note: This project is still under development. Please refer to the code comments in the dev branch. Technical documentation is for reference only.

## Prerequisites
- uv (Python project management tool)

## Get the Source Code
It is recommended to obtain the repository via:
- Clone (Recommended):
  ```bash
  git clone https://github.com/NeoNexusX/MassFlow.git
  cd MassFlow
  ```
- Fork then Clone (For contributing):
  ```bash
  git fork https://github.com/NeoNexusX/MassFlow.git
  # Then clone your Fork
  ```
- Download ZIP: Click GitHub "Code" → "Download ZIP" and extract locally.

## Set up Python Environment

```bash
close conda first:
  conda deactivate

uv install part:
  https://docs.astral.sh/uv/getting-started/installation/

For Example :
  Linux && Macos :
  curl -LsSf https://astral.sh/uv/install.sh | sh 
  Windows:
  powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

Install dependencies:
uv sync 

uv pip install -e .
```

## Run Examples and Tests

Tests use the pytest framework. See: https://docs.pytest.org/en/stable/

Run all tests:

```bash
uv run pytest
```

Run a specific test:

```bash
uv run pytest tests/test_read.py
```

If the virtual environment managed by uv is already activated, you can omit `uv run`:

```bash
pytest tests/test_read.py
```

Example Code Snippet:

```python
# Quick snippet: Read .imzML and plot occupancy mask
from massflow.module.ms_module import MS
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML

FILE_PATH = "data/your_file.imzML"
ms = MS()
with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:
    manager.load_full_data_from_file()
    manager.inspect_data()
    ms.plot_ms_mask()
```

The code above demonstrates MSI/MS data reading and basic visualization. Logs are output to `logs/`.

## Common Commands
Documentation build commands:
  - Build static site: `npm run docs:build`
  - Preview built site: `npm run docs:preview`

## Next Steps
- Contribution Guide: `/en/contribution`
- Naming Conventions: `/en/naming-conventions`
- Data Structures: `/en/ms-data-structures`
- Noise Reduction: `/en/noise_reduction`
- Baseline Correction: `/en/baseline_correction`
- Collaboration Guide: `/en/collaboration_guide`
- More documents...