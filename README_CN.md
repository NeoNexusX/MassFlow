# MassFlow

[English](README.md) | 简体中文

MassFlow 是一个面向质谱成像（MSI）与质谱（MS）数据的模块化预处理与数据管理框架。目前支持：
- 读取 imzML（MS 成像光谱，按像素惰性加载）
- 光谱去噪/平滑：移动平均（MA）、高斯（Gaussian）、双边（Bilateral），以及基于邻域搜索（kNN）的变体
- MSI 数据的基本管理与写出（.msi/.h5，按 m/z 分组写入）


## 安装

要求：Python >= 3.9（已在 3.12 上测试）

```bash
# 克隆仓库
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow

# 安装依赖
pip install -r requirements.txt
```

## 快速开始

推荐直接运行示例：

```bash
python example.py
```

## 项目结构

```
MassFlow/
├── example.py
├── logger.py
├── module/
│   ├── __init__.py
│   ├── ms_data_manager.py           # MS 基类数据管理器
│   ├── ms_module.py                 # MS/ImzML 基础类型（光谱、惰性加载）
│   ├── msi_data_manager.py          # MSI 抽象/通用管理器
│   ├── msi_data_manager_msi.py      # MSI 写出（.msi/.h5）
│   ├── msi_data_manager_zys.py      # ZYS (.mat) 数据管理
│   ├── msi_module.py                # MSI 领域模型与可视化
│   └── __pycache__/...
├── preprocess/
│   ├── filter.py                    # 去噪/平滑函数（MA/高斯/双边 + NS 变体）
│   └── ms_preprocess.py             # 预处理入口（噪声抑制 API 等）
├── data/
│   ├── example.imzML
│   └── example.ibd
├── docs/
│   ├── CONTRIBUTING.md / EN.md
│   ├── NAMING_CONVENTIONS.md / EN.md
│   ├── Collaboration_Guide.md
│   └── 协作指北.md
├── logs/
│   └── *.log
├── requirements.txt
├── LICENSE
└── README*.md
```

与此前文档不同点：
- 预处理模块文件名为 `preprocess/ms_preprocess.py`（非 `msi_preprocess.py`）
- MS（imzML）读取在 `module/ms_data_manager_imzml.py` 中（示例已导入）
- 当前仓库未包含 `.github/`、`.vscode/`、`tests/` 等目录

## 开发与贡献

- 贡献指南：见 `docs/CONTRIBUTING.md` 与 `docs/CONTRIBUTING_EN.md`
- 命名规范：见 `docs/NAMING_CONVENTIONS*.md`
- 推荐扩展：Python、Pylance、Ruff、Black、isort、Pylint、Markdownlint、GitLens、H5Web

## 许可证

本项目采用 GNU 通用公共许可证 v3.0 - 详见 [LICENSE](LICENSE)。

## 参考资料

- Cardinal MSI: https://cardinalmsi.org/
- MATLAB 质谱预处理: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- PyOpenMS: https://pyopenms.readthedocs.io/

## 反馈

如需支持或发现问题，请提交 GitHub Issue。
