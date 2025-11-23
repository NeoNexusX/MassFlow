# MassFlow

[English](README.md) | 简体中文

MassFlow 是一个面向质谱成像（MSI）与质谱（MS）数据的模块化预处理与数据管理框架。目前支持：

## 安装

要求：Python >= 3.12（建议使用 3.12 版本）

```bash
# 克隆仓库
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
```

## 快速开始

推荐在 Jupyter 打开 `example.ipynb`，或直接运行以下代码片段验证数据读取：

```python
from module.ms_module import MS
from module.ms_data_manager_imzml import MSDataManagerImzML

FILE_PATH = "data/your_file.imzML"
ms = MS()
with MSDataManagerImzML(ms=ms, target_locs=[(1, 1), (50, 50)], filepath=FILE_PATH) as manager:
    manager.load_full_data_from_file()
    manager.inspect_data()
    ms.plot_ms_mask()
```

在线文档: https://neonexusx.github.io/MassFlow/

## 项目结构

```
MassFlow/
├── .github/                 # GitHub 配置（Issue 模板、Workflows）
├── .vscode/                 # VSCode 配置
├── data/                    # 示例数据
├── docs/                    # 文档（VitePress）
│   ├── en/
│   └── zh/
├── logs/                    # 运行日志
├── src/
│   └── massflow/            # 核心源码
│       ├── module/          # 数据模型与管理器
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
├── pyproject.toml           # 项目配置与依赖
└── README_CN.md
```

## 开发与贡献

- 贡献指南：`docs/zh/contribution.md` 与 `docs/en/contribution.md`
- 命名规范：`docs/zh/naming-conventions.md` 与 `docs/en/naming-conventions.md`
- Issue 模板：`.github/ISSUE_TEMPLATE/feature.md`、`bug.md`、`feature_en.md`、`bug_en.md`
- 本地检查：`ruff .`、`black .`、`isort .`、`pylint module/`
- 提交规范：Conventional Commits（如 `feat:`、`fix:`、`docs:`、`refactor:`、`test:`）
- 推荐扩展：Python、Pylance、Ruff、Black、isort、Pylint、Markdownlint、GitLens、H5Web

## 许可证

本项目采用 GNU 通用公共许可证 v3.0 - 详见 [LICENSE](LICENSE)。

## 参考资料

- Cardinal MSI: https://cardinalmsi.org/
- MATLAB 质谱预处理: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- PyOpenMS: https://pyopenms.readthedocs.io/

## 反馈

如需支持或发现问题，请提交 GitHub Issue。
