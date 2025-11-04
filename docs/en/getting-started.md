---
title: Getting Started
---

# Getting Started

MassFlow is a modular framework for Mass Spectrometry Imaging (MSI) and Mass Spectrometry (MS) data. It provides efficient data loading, preprocessing (denoising/smoothing), and basic data management for MSI workflows.

## Prerequisites
- Python `>= 3.9` (tested on `3.12`)

## Get the Source Code
You can obtain the repository in one of the following ways:
- Clone (recommended):
  ```bash
  git clone https://github.com/NeoNexusX/MassFlow.git
  cd MassFlow
  ```
- Fork then clone your fork (for contributions):
  ```bash
  git fork https://github.com/NeoNexusX/MassFlow.git
  # then clone your fork
  ```
- Download ZIP: click "Code" → "Download ZIP" on GitHub and extract locally.

## Set Up Python Environment
It is recommended to use a virtual environment on macOS:

```bash
# In the repository root
python3 -m venv .venv
source .venv/bin/activate  # macOS/Linux

# Install dependencies
pip install -r requirements.txt
```

If you use `conda`:
```bash
conda create -n massflow python=3.12 -y
conda activate massflow
pip install -r requirements.txt
```

## Quick Start (Run the Example)
Run the bundled example to verify your environment:

```bash
python example.py
```

This script demonstrates reading MSI/MS data and basic preprocessing. Logs are written to `logs/`.


Other useful commands:
- Build static site: `npm run docs:build`
- Preview built site: `npm run docs:preview`

## Troubleshooting
- Ensure the correct Python version (`python3 --version`) and that your virtual environment is activated (`which python`).
- If `pip` fails due to SSL or networking, try `pip install -r requirements.txt --no-cache-dir` or set a trusted mirror.
- On Apple Silicon (M1/M2/M3), prefer Python 3.11+ and native wheels; update `pip` with `python -m pip install --upgrade pip`.

## Next Steps
- Contribution guidelines: `docs/CONTRIBUTING_EN.md` / `docs/CONTRIBUTING.md`
- Naming conventions: `docs/NAMING_CONVENTIONS_EN.md` / `docs/NAMING_CONVENTIONS.md`
- Collaboration guide: `docs/Collaboration_Guide.md` / `docs/协作指北.md`