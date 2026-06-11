---
title: Claude Code Integration
---

# Claude Code Integration

MassFlow exposes its analysis tools as an MCP (Model Context Protocol) server, allowing Claude Code to call MassFlow functions directly during a coding session.

## Prerequisites

- Python >= 3.12
- MassFlow installed in editable mode with MCP support:
  ```bash
  pip install -e ".[mcp]"
  ```
- [Claude Code](https://docs.anthropic.com/en/docs/claude-code) installed and available on PATH

## Quick Start

### 1. Generate the MCP config

From the MassFlow project root, run:

```bash
python -m massflow.mcp_server.server --print-claude-config
```

This outputs a JSON snippet like:

```json
{
  "mcpServers": {
    "massflow": {
      "command": "/path/to/python",
      "args": ["-m", "massflow.mcp_server.server"]
    }
  }
}
```

### 2. Add the config to Claude Code

**Option A** -- Project-level (recommended): save the output as `.mcp.json` in the project root:

```bash
python -m massflow.mcp_server.server --print-claude-config > .mcp.json
```

Claude Code picks up `.mcp.json` automatically when started from this directory.

**Option B** -- Global: paste the `mcpServers` block into your Claude Code settings.

### 3. Start Claude Code

```bash
cd /path/to/MassFlow
claude
```

Claude Code will start the MassFlow MCP server automatically via stdio. You can then ask Claude to inspect datasets, run preprocessing, etc.

## Available MCP Tools

Once connected, Claude Code has access to these MassFlow tools:

| Tool | Description |
|------|-------------|
| `load_imzml` | Load imzML header and dataset index |
| `dataset_qc_summary` | Build a sampled-spectrum QC summary |
| `recommend_preprocessing_plan` | Recommend preprocessing operations |
| `run_preprocessing` | Run preprocessing pipeline |
| `build_spatial_preview` | Generate TIC and base-peak spatial preview images |
| `apply_spatial_png_mask` | Apply a PNG mask to spatial data |
| `analyze_imzml_end_to_end` | Full end-to-end workflow (load, preprocess, spatial, annotate, enrich) |
| `validate_integrated_conclusion` | Validate the integrated conclusion |

## Example Usage in Claude Code

Once connected, you can interact naturally:

```
> Load the imzML file at data/sample.imzML and show me the dataset summary

> Run preprocessing on data/sample.imzML with noise reduction and normalization

> Run the full end-to-end workflow on data/sample.imzML
```

## Configuration

MassFlow reads optional runtime config from `massflow.toml` (or `.massflow.toml`) in the project root. See `massflow.toml.example` for available options:

```toml
[agent_tools]
default_output_dir = "outputs"
max_threads = 8

[mcp_server]
transport = "stdio"
host = "127.0.0.1"
port = 8000
```

When `output_dir` is not specified, outputs are written to `massflow_output/` next to the input imzML file.

## Manual Debug Mode (HTTP Transport)

For debugging, you can start the MCP server manually over HTTP:

```bash
python -m massflow.mcp_server.server --transport streamable-http --host 127.0.0.1 --port 8000
```

Then connect your MCP client to `http://127.0.0.1:8000`.

## Troubleshooting

**`ModuleNotFoundError: No module named 'mcp'`**

Install the MCP dependency:

```bash
pip install -e ".[mcp]"
```

**Server starts but tools don't appear in Claude Code**

- Verify `.mcp.json` is in the project root and contains the correct `mcpServers` block
- Check that the `command` path points to the Python interpreter where MassFlow is installed
- Restart Claude Code after adding the config

**Check your setup**

```bash
python -m massflow.agent_tools doctor --require-mcp
```

This runs a self-check and reports any missing dependencies.
