# MassFlow

[English](README.md) | 简体中文

MassFlow 是一个用于质谱成像 (MSI) 和质谱 (MS) 数据的模块化预处理和数据管理框架。目前支持：

## 安装

要求：Python >= 3.12（推荐使用 3.12 版本）

```bash
# 克隆仓库
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
```

## 快速开发入门

建议直接运行以下代码 (main.py) 来验证数据加载：

```bash
# 首先关闭 conda：
conda deactivate

# 安装 uv：
# https://docs.astral.sh/uv/getting-started/installation/

# 例如：
# Linux && macOS:
curl -LsSf https://astral.sh/uv/install.sh | sh 
# Windows:
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# 安装依赖：
uv sync && uv pip install -e .
```

在线文档：https://neonexusx.github.io/MassFlow/

## 项目结构

```
MassFlow/
├── .github/                 # GitHub 配置 (Issue 模板, Workflows)
├── .vscode/                 # VSCode 配置
├── data/                    # 示例数据
├── docs/                    # 文档 (VitePress)
│   ├── en/
│   └── zh/
├── logs/                    # 运行时日志
├── src/
│   └── massflow/            # 核心源代码
│       ├── module/          # 数据模型和管理器
│       │   ├── ms_module.py # MS/ImzML 基础类型
│       │   └── ...
│       ├── preprocess/      # 预处理算法
│       │   ├── ms_preprocess.py
│       │   └── ...
│       └── tools/           # 工具函数
├── tests/                   # 测试用例
├── LICENSE
├── main.py
├── package.json
├── pyproject.toml           # 项目配置和依赖
└── README.md
```

## 开发与贡献

- 贡献指南：`docs/en/contribution.md` 和 `docs/zh/contribution.md`
- 命名规范：`docs/en/naming-conventions.md` 和 `docs/zh/naming-conventions.md`
- Issue 模板：`.github/ISSUE_TEMPLATE/feature_en.md`, `bug_en.md`, `feature.md`, `bug.md`
- 本地检查：`ruff .`, `black .`, `isort .`, `pylint module/`
- 提交规范：Conventional Commits (例如：`feat:`, `fix:`, `docs:`, `refactor:`, `test:`)
- 推荐扩展：Python, Pylance, Ruff, Black, isort, Pylint, Markdownlint, GitLens, H5Web

## 许可证

本项目采用 GNU General Public License v3.0 许可证 - 详情请参阅 [LICENSE](LICENSE) 文件。

## 参考资料

- Cardinal MSI: https://cardinalmsi.org/
- MATLAB Mass Spectrometry Preprocessing: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- PyOpenMS: https://pyopenms.readthedocs.io/

## 反馈

如需支持或报告问题，请提交 GitHub Issue。
