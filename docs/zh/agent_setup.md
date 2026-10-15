# 在 opencode 中安装 MassFlow 智能体（Skills + MCP）

MassFlow 的智能体由两部分组成，两者都需要安装：

| 组件 | 作用 | 位置 |
| --- | --- | --- |
| **mzorbit MCP 服务** | 提供 27 个工具：读取 imzML、质控、预处理、聚类、化合物注释 | `mzorbit/`（git 子模块） |
| **3 个 Skills** | 提供工作流和决策规则：什么时候调用哪个工具，什么时候先问你 | `.claude/skills/` |

只装 MCP 时，模型能调用工具，但不知道正确的流程；只装 Skills 时，Skills 让模型调用的工具并不存在。

本文以 opencode 为默认客户端，所有步骤已在 opencode 1.17.20 上验证。下文的 `/path/to/MassFlow` 请替换为你本机仓库的绝对路径。

## 1. 准备环境

需要事先安装 [git](https://git-scm.com/)、[uv](https://docs.astral.sh/uv/) 和 [opencode](https://opencode.ai/)。

```bash
git clone --recurse-submodules https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
uv sync --all-packages --extra mcp
```

如果 clone 时漏了 `--recurse-submodules`（`mzorbit/` 目录为空），补一条 `git submodule update --init`，再运行上面的 `uv sync`。

::: warning 只用这一条 sync 命令
单独运行 `uv sync`（不带 `--all-packages`）只会安装 `massflow`，并把 `mzorbit` 和 `mcp` 从 `.venv` 中**卸载**，MCP 服务随之失效。以后每次同步依赖都用 `uv sync --all-packages --extra mcp`。
:::

检查安装结果：

```bash
.venv/bin/python -c "import mzorbit, mcp; print('ok')"
```

## 2. 注册 MCP 服务

### 生成配置

在仓库根目录运行：

```bash
uv run mzorbit-mcp-server --print-config opencode
```

输出类似：

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

- `command`：本环境 Python 的绝对路径。opencode 直接启动它，不需要先激活虚拟环境，也不要手动启动服务。
- `timeout`：1800000 毫秒，即 30 分钟。在 opencode 1.17 中，这个值同时是**工具调用**的超时。不设置时，每次调用 60 秒就会被中断，而对整套数据运行 `finalize_preprocessing`、`run_umap_clustering`、`export_dataset` 通常需要几分钟到几十分钟。

### 写入配置文件

| 你在哪里启动 opencode | 写入的文件 |
| --- | --- |
| 只在 MassFlow 仓库目录内 | 仓库根目录的 `opencode.json`（已在 `.gitignore` 中，不会被提交） |
| 任意目录，例如数据目录 | 全局配置 `~/.config/opencode/opencode.json` |

目标文件还不存在时，可以直接写入，例如：

```bash
uv run mzorbit-mcp-server --print-config opencode > opencode.json
```

目标文件已经存在时，**不要**用 `>` 覆盖它。只需把输出中 `"mcp"` 下面的 `"mzorbit"` 这一段合并进去。

### 验证

```bash
opencode mcp list
```

应看到 `✓ mzorbit connected`。

## 3. 安装 Skills

| Skill | 何时触发 |
| --- | --- |
| `msi-dataset-evaluation` | 拿到新数据、做质控、比较处理前后的数据 |
| `msi-preprocessing` | 预处理：峰提取、对齐、归一化、去噪、去背景、导出 |
| `msi-spatial-analysis` | 空间分割、区域标志离子、化合物与文献检索 |

### 在仓库内使用：无需安装

opencode 会从启动目录向上查找，直到 git 仓库根目录，并自动加载其中 `.claude/skills/` 下的 Skills。

### 在仓库外使用：加入全局配置

在 `~/.config/opencode/opencode.json` 中加入下面这段，并与已有内容合并：

```json
{
  "skills": {
    "paths": ["/path/to/MassFlow/.claude/skills"]
  }
}
```

这里指向的是仓库本身，而不是一份拷贝，所以 `git pull` 之后 Skills 会自动更新。在仓库内启动时，这个配置和自动发现可以共存，不会重复加载。

如果更习惯使用软链接，也可以这样做（与上面的方式二选一）：

```bash
mkdir -p ~/.config/opencode/skills
ln -sfn /path/to/MassFlow/.claude/skills/msi-* ~/.config/opencode/skills/
```

### 验证

```bash
opencode debug skill > /tmp/opencode-skills.json
grep -o '"location": "[^"]*msi-[^"]*"' /tmp/opencode-skills.json
```

应列出 3 个 `SKILL.md`。这条命令的输出很大，请务必先写入文件再查找：在 opencode 1.17 中，直接用管道接 `grep` 时只会传出前 64 KB，`msi-*` 可能恰好被截掉。

修改配置或 Skills 后，需要**重启 opencode**，配置不会热加载。

## 4. 开始使用

在 opencode 中直接用自然语言描述任务即可。模型会根据描述加载对应的 Skill，再调用 mzorbit 工具。例如：

> 帮我看看 /data/msi/kidney.imzML 的数据质量，能不能用来做脂质成像

> 把 /data/msi/kidney.imzML 做一下预处理，峰提取加对齐，归一化我想自己选

> 找出 /data/msi/kidney_aligned.imzML 里空间上不同的区域，以及标志每个区域的离子

使用时注意：

- imzML 请写绝对路径，同名的 `.ibd` 文件必须放在同一目录下。
- 默认是引导模式：归一化、基线校正、峰过滤、组织掩膜由你决定。如果想全自动处理，直接说"全自动处理，不要问我"，模型会自行决定，并逐项汇报每个选择及其依据。
- 输出永远不会覆盖输入。每个输出文件旁边的 `<文件名>.metadata.json` 记录了离子极性、生物学元数据和完整的处理历史。

## 5. 可选：调整工具默认值

在仓库根目录创建 `massflow.toml`，可以修改所有工具的默认参数：

```toml
[agent_tools]
max_threads = 8                  # 各工具默认线程数
batch_size = 256                 # finalize / 组织掩膜 / 导出 的默认批大小
default_output_dir = "./output"  # 默认输出目录，相对于本文件
```

请先创建这个文件，再按第 2 步生成配置：生成的命令会自动带上 `--config /path/to/MassFlow/massflow.toml`。

## 6. 更新

```bash
cd /path/to/MassFlow
git pull --recurse-submodules
uv sync --all-packages --extra mcp
```

然后重启 opencode。

## 7. 常见问题

| 现象 | 原因与处理 |
| --- | --- |
| `opencode mcp list` 显示 mzorbit 失败 | 最常见的原因是 `.venv` 中缺少 `mzorbit` 或 `mcp`，重新运行 `uv sync --all-packages --extra mcp` 即可。如果移动过仓库或重建过 `.venv`，`command` 中的 Python 路径会失效，需要按第 2 步重新生成配置。也可以在终端直接运行 `command` 中的命令，查看具体报错 |
| 工具运行约 60 秒后报超时 | 配置中缺少 `"timeout"`，按第 2 步补上 `1800000` |
| 模型没有按流程工作，没有加载 Skill | opencode 在仓库外启动，但没有配置全局 Skills（见第 3 步）；或者修改配置后没有重启 opencode |
| 提示 polarity unknown | 处理后的文件不再记录离子极性。告诉模型是正离子还是负离子模式，它会用 `set_acquisition` 记录下来 |
| 预处理中途卡住或内存不足 | `finalize_preprocessing` 已经会自动用 1/4 的批大小重试一次。仍然失败时，让模型用更小的 `batch_size`（如 64）或更少的 `max_threads` 重新运行 |

## 其他客户端

同一个 MCP 服务也可以用于其他客户端，只需换一个参数生成配置：

```bash
uv run mzorbit-mcp-server --print-config claude   # Claude Code：写入仓库根目录的 .mcp.json
uv run mzorbit-mcp-server --print-config codex    # Codex：写入 .codex/config.toml
```

Claude Code 会原生读取 `.claude/skills/`，不需要另外安装 Skills。
