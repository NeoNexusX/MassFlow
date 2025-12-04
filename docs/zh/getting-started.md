# 快速开始

MassFlow 是一个面向质谱成像（MSI）与质谱（MS）数据的模块化计算框架，提供高效的数据读取、预处理，与管理功能。
注意：本项目还在开发中，请以dev分支的代码注释为准，技术文档仅作参考。

## 前置条件
- uv  python项目管理软件

## 获取源码
推荐通过以下方式获取仓库：
- 克隆（推荐）：
  ```bash
  git clone https://github.com/NeoNexusX/MassFlow.git
  cd MassFlow
  ```
- 先 Fork 再克隆（用于贡献代码）：
  ```bash
  git fork https://github.com/NeoNexusX/MassFlow.git
  # 然后克隆你的 Fork
  ```
- 下载 ZIP：点击 GitHub “Code” → “Download ZIP” 并在本地解压。

## 设置 Python 环境

```bash
close conda first:
  conda deactivate

uv install part:
  https://uv.doczh.com/getting-started/installation/

For Example :
  Linux && Macos :
  curl -LsSf https://astral.sh/uv/install.sh | sh 
  Windows:
  powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

Install dependencies:
uv sync 

uv pip install -e .
```

## 运行示例 和测试：

测试使用了pytest框架，详见：https://docs.pytest.org/en/stable/ 或者 中文：https://docs.pytest.cn/en/stable/

运行所有测试：

```bash
uv run pytest
```

运行某一个测试：

```python
uv run pytest tests/test_read.py
```

如果已经启动了由uv管理的虚拟环境则可以省略uv run：

```python
pytest tests/test_read.py
```

测试样例代码：

```python
# 快速片段：读取 .imzML 并绘制占用掩码
from massflow.module.ms_module import MS
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML

FILE_PATH = "data/your_file.imzML"
ms = MS()
with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:
    manager.load_full_data_from_file()
    manager.inspect_data()
    ms.plot_ms_mask()
```

上述代码演示 MSI/MS 数据读取与基础可视化。日志输出到 `logs/`。

## 常用命令
文档构建命令：
  - 构建静态站点：`npm run docs:build`
  - 预览已构建站点：`npm run docs:preview`

## 下一步
- 贡献说明：`/zh/contribution`
- 命名规范：`/zh/naming-conventions`
- 数据结构：`/zh/ms-data-structures`
- 噪声抑制：`/zh/noise_reduction`
- 基线校正：`/zh/baseline_correction`
- 协作指南：`/zh/collaboration_guide`
- 等等文档