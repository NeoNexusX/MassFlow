---
title: Claude Code 集成
---

# Claude Code 集成

MassFlow 将分析工具暴露为 MCP（Model Context Protocol）服务器，使 Claude Code 能在编码会话中直接调用 MassFlow 函数。

## 前置条件

- Python >= 3.12
- 以 MCP 支持模式安装 MassFlow：
  ```bash
  pip install -e ".[mcp]"
  ```
- 已安装 [Claude Code](https://docs.anthropic.com/en/docs/claude-code) 并在 PATH 中可用

## 快速开始

### 1. 生成 MCP 配置

在 MassFlow 项目根目录下运行：

```bash
python -m massflow.mcp_server.server --print-claude-config
```

输出如下 JSON 片段：

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

### 2. 将配置添加到 Claude Code

**方式 A** — 项目级配置（推荐）：将输出保存为项目根目录下的 `.mcp.json`：

```bash
python -m massflow.mcp_server.server --print-claude-config > .mcp.json
```

Claude Code 在此目录启动时会自动读取 `.mcp.json`。

**方式 B** — 全局配置：将 `mcpServers` 块粘贴到 Claude Code 的全局设置中。

### 3. 启动 Claude Code

```bash
cd /path/to/MassFlow
claude
```

Claude Code 会通过 stdio 自动启动 MassFlow MCP 服务器。之后即可让 Claude 检查数据集、运行预处理等。

## 可用的 MCP 工具

连接后，Claude Code 可使用以下 MassFlow 工具：

| 工具 | 说明 |
|------|------|
| `load_imzml` | 加载 imzML 头部和数据集索引 |
| `dataset_qc_summary` | 构建抽样谱 QC 摘要 |
| `recommend_preprocessing_plan` | 推荐预处理操作 |
| `run_preprocessing` | 运行预处理流水线 |
| `build_spatial_preview` | 生成 TIC 和基峰空间预览图 |
| `apply_spatial_png_mask` | 对空间数据应用 PNG 掩模 |
| `analyze_imzml_end_to_end` | 完整端到端流程（加载、预处理、空间分析、注释、富集） |
| `validate_integrated_conclusion` | 验证综合结论 |

## 在 Claude Code 中的使用示例

连接后可直接自然语言交互：

```
> 加载 data/sample.imzML 并展示数据集摘要

> 对 data/sample.imzML 运行去噪和归一化预处理

> 对 data/sample.imzML 运行完整端到端流程
```

## 配置

MassFlow 从项目根目录的 `massflow.toml`（或 `.massflow.toml`）读取可选运行时配置。可用选项参见 `massflow.toml.example`：

```toml
[agent_tools]
default_output_dir = "outputs"
max_threads = 8

[mcp_server]
transport = "stdio"
host = "127.0.0.1"
port = 8000
```

未指定 `output_dir` 时，输出写入 imzML 文件同级的 `massflow_output/` 目录。

## 手动调试模式（HTTP 传输）

调试时可手动启动 MCP 服务器：

```bash
python -m massflow.mcp_server.server --transport streamable-http --host 127.0.0.1 --port 8000
```

然后将 MCP 客户端连接到 `http://127.0.0.1:8000`。

## 常见问题

**`ModuleNotFoundError: No module named 'mcp'`**

安装 MCP 依赖：

```bash
pip install -e ".[mcp]"
```

**服务器启动但 Claude Code 中看不到工具**

- 确认 `.mcp.json` 在项目根目录且包含正确的 `mcpServers` 块
- 检查 `command` 路径指向安装了 MassFlow 的 Python 解释器
- 添加配置后重启 Claude Code

**检查安装状态**

```bash
python -m massflow.agent_tools doctor --require-mcp
```

运行自检并报告缺失的依赖。
