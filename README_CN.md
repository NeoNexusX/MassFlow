# MassFlow

[English](README.md) | 简体中文

MassFlow 是一个面向质谱成像（MSI）与质谱（MS）数据的模块化预处理与数据管理框架。目前支持：
- 读取 imzML（MS 成像质谱，按像素惰性加载）
- 质谱去噪/平滑：移动平均（MA）、高斯（Gaussian）、双边（Bilateral），以及基于邻域搜索（kNN）的变体
- 基线校正、归一化、峰提取/峰对齐
- MSI 数据管理与写出（.msi/.h5，按 m/z 分组写入）


## 安装

要求：Python >= 3.9（建议使用 3.12 版本）

```bash
# 克隆仓库
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow

# 安装依赖
pip install -r requirements.txt
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
├── module/
│   ├── __init__.py
│   ├── meta_data.py
│   ├── ms_data_manager.py
│   ├── ms_data_manager_imzml.py
│   ├── ms_module.py                 # MS/ImzML 基础类型（质谱、惰性加载）
│   ├── msi_data_manager_base.py
│   ├── msi_data_manager_msi.py
│   ├── msi_data_manager_zys.py
│   └── msi_module.py
├── preprocess/
│   ├── ms_preprocess.py
│   ├── filter_helper.py
│   ├── baseline_correction_helper.py
│   ├── normalizer_helper.py
│   ├── peak_alignment.py
│   └── peak_pick_helper.py
├── docs/
│   ├── en/...
│   └── zh/...
├── .github/
│   └── ISSUE_TEMPLATE/...
├── tools/
│   └── plot.py
├── example.ipynb
├── example.py
├── logger.py
├── requirements.txt
├── LICENSE
└── README*.md
```

与此前文档不同点：
- 预处理模块入口为 `preprocess/ms_preprocess.py`（非 `msi_preprocess.py`）
- imzML 读取在 `module/ms_data_manager_imzml.py` 中
- 仓库包含 `.github/`，不包含 `.vscode/` 与 `tests/` 目录

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
