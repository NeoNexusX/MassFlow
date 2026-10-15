# Installing the MassFlow Agent in opencode (Skills + MCP)

The MassFlow agent has two parts, and you need both:

| Component | What it does | Location |
| --- | --- | --- |
| **mzorbit MCP server** | 27 tools: read imzML, QC, preprocessing, clustering, compound annotation | `mzorbit/` (git submodule) |
| **3 skills** | Workflows and decision rules: which tool to call when, and when to ask you first | `.claude/skills/` |

With only the MCP server, the model can call the tools but does not know the workflow. With only the skills, the tools the skills tell the model to call do not exist.

This guide uses opencode as the default client. Every step was verified on opencode 1.17.20. Replace `/path/to/MassFlow` with the absolute path of your local clone.

## 1. Set up the environment

Install [git](https://git-scm.com/), [uv](https://docs.astral.sh/uv/) and [opencode](https://opencode.ai/) first.

```bash
git clone --recurse-submodules https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
uv sync --all-packages --extra mcp
```

If you cloned without `--recurse-submodules` (the `mzorbit/` folder is empty), run `git submodule update --init`, then run the `uv sync` command above.

::: warning Use only this sync command
A plain `uv sync` (without `--all-packages`) installs only `massflow` and **uninstalls** `mzorbit` and `mcp` from `.venv`, which breaks the MCP server. Always sync with `uv sync --all-packages --extra mcp`.
:::

Check the installation:

```bash
.venv/bin/python -c "import mzorbit, mcp; print('ok')"
```

## 2. Register the MCP server

### Generate the config

From the repository root, run:

```bash
uv run mzorbit-mcp-server --print-config opencode
```

The output looks like this:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "mzorbit": {
      "type": "local",
      "command": ["/path/to/MassFlow/.venv/bin/python", "-m", "mzorbit.interfaces.mcp_server"],
      "enabled": true,
      "timeout": 1800000
    }
  }
}
```

- `command`: the absolute path of this environment's Python. opencode starts it directly, so you do not need to activate the virtual environment. Do not start the server yourself.
- `timeout`: 1,800,000 ms, which is 30 minutes. In opencode 1.17 this value is also the timeout for **tool calls**. Without it, every call is cut off after 60 seconds, while running `finalize_preprocessing`, `run_umap_clustering` or `export_dataset` on a full dataset usually takes several minutes to tens of minutes.

### Write the config file

| Where you start opencode | File to write |
| --- | --- |
| Only inside the MassFlow repository | `opencode.json` at the repository root (already in `.gitignore`, so it is never committed) |
| Anywhere, for example a data folder | Global config `~/.config/opencode/opencode.json` |

If the target file does not exist yet, you can write it directly, for example:

```bash
uv run mzorbit-mcp-server --print-config opencode > opencode.json
```

If the target file already exists, **do not** overwrite it with `>`. Merge only the `"mzorbit"` entry under `"mcp"` into it.

### Verify

```bash
opencode mcp list
```

You should see `✓ mzorbit connected`.

## 3. Install the skills

| Skill | When it triggers |
| --- | --- |
| `msi-dataset-evaluation` | New data arrives, quality checks, comparing data before and after processing |
| `msi-preprocessing` | Preprocessing: peak picking, alignment, normalization, denoising, background removal, export |
| `msi-spatial-analysis` | Spatial segmentation, region marker ions, compound and literature searches |

### Inside the repository: nothing to install

opencode walks up from the folder it was started in to the git repository root and loads the skills in `.claude/skills/` automatically.

### Outside the repository: add them to the global config

Add the following to `~/.config/opencode/opencode.json`, merged with what is already there:

```json
{
  "skills": {
    "paths": ["/path/to/MassFlow/.claude/skills"]
  }
}
```

This points at the repository itself rather than a copy, so the skills update with every `git pull`. When you start opencode inside the repository, this setting and the automatic discovery coexist without loading anything twice.

If you prefer symlinks, this works too (use one method or the other):

```bash
mkdir -p ~/.config/opencode/skills
ln -sfn /path/to/MassFlow/.claude/skills/msi-* ~/.config/opencode/skills/
```

### Verify

```bash
opencode debug skill > /tmp/opencode-skills.json
grep -o '"location": "[^"]*msi-[^"]*"' /tmp/opencode-skills.json
```

This should list three `SKILL.md` files. The command's output is large, so write it to a file before searching: in opencode 1.17, piping it straight into `grep` passes only the first 64 KB, which can cut off the `msi-*` entries.

After changing the config or the skills, **restart opencode**. The config is not reloaded while it runs.

## 4. Start using it

Describe the task in plain language in opencode. The model loads the matching skill and then calls the mzorbit tools. For example:

> Check the data quality of /data/msi/kidney.imzML. Is it good enough for lipid imaging?

> Preprocess /data/msi/kidney.imzML with peak picking and alignment. I want to choose the normalization myself.

> Find the spatially distinct regions in /data/msi/kidney_aligned.imzML and the ions that mark each region.

Keep in mind:

- Give the imzML as an absolute path. The `.ibd` file with the same name must be in the same folder.
- Guided mode is the default: you decide on normalization, baseline correction, peak filters and tissue masks. For fully automatic processing, say "process it fully automatically, don't ask me". The model then decides on its own and reports each choice with its evidence.
- Outputs never overwrite inputs. The `<file name>.metadata.json` next to each output records the ion polarity, the biological metadata and the full processing history.

## 5. Optional: change the tool defaults

Create `massflow.toml` at the repository root to change the defaults of every tool:

```toml
[agent_tools]
max_threads = 8                  # default thread count for every tool
batch_size = 256                 # default batch size for finalize / tissue mask / export
default_output_dir = "./output"  # default output folder, relative to this file
```

Create this file first, then generate the config as in step 2: the generated command then includes `--config /path/to/MassFlow/massflow.toml` automatically.

## 6. Update

```bash
cd /path/to/MassFlow
git pull --recurse-submodules
uv sync --all-packages --extra mcp
```

Then restart opencode.

## 7. Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| `opencode mcp list` shows mzorbit as failed | Most often `mzorbit` or `mcp` is missing from `.venv`: run `uv sync --all-packages --extra mcp` again. If you moved the repository or rebuilt `.venv`, the Python path in `command` is stale: regenerate the config as in step 2. Running the `command` directly in a terminal shows the exact error |
| A tool times out after about 60 seconds | The config has no `"timeout"`: add `1800000` as in step 2 |
| The model does not follow the workflow or load a skill | opencode was started outside the repository without the global skills config (see step 3), or opencode was not restarted after a config change |
| "polarity unknown" | Processed files no longer record the ion polarity. Tell the model whether the data are positive or negative mode; it records this with `set_acquisition` |
| Preprocessing stalls or runs out of memory | `finalize_preprocessing` already retries once with a quarter of the batch size. If it still fails, ask the model to rerun with a smaller `batch_size` (e.g. 64) or fewer `max_threads` |

## Other clients

The same MCP server works with other clients; generate their config with a different argument:

```bash
uv run mzorbit-mcp-server --print-config claude   # Claude Code: write to .mcp.json at the repository root
uv run mzorbit-mcp-server --print-config codex    # Codex: write to .codex/config.toml
```

Claude Code reads `.claude/skills/` natively, so its skills need no separate installation.
