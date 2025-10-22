# MassFlow

[English](README.md) | 简体中文

MassFlow 是一个模块化的高性能质谱成像（MSI）数据预处理框架。它为 MSI 研究提供高效的数据管理、处理和可视化功能。


## 安装

```bash
# 克隆仓库
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow

# 安装依赖
pip install -r requirements.txt
```

## 快速开始

### 示例 1：加载完整 MSI 数据

```python
from module.msi_data_manager import MSIDataManager
from module.msi_module import MSI

# Create MSI instance
msi = MSI(name='example1', version=1.0, mask=None, need_base_mask=True)

# Initialize data manager
msi_dm = MSIDataManager(msi, filepath="./data/your_data.mat")

# Load all data
msi_dm.load_full_data_from_file()

# Inspect data
msi_dm.inspect_data()
```

### 示例 2：加载特定 m/z 范围

```python
# Load only m/z values between 101 and 102
msi2 = MSI(name='example2', version=1.0, mask=None, need_base_mask=True)
msi_dm2 = MSIDataManager(msi2, target_mz_range=[101, 102], filepath="./data/your_data.mat")
msi_dm2.load_full_data_from_file()

# Plot images in the loaded range
msi2.plot_msi()
```

## 项目结构

```
MassFlow/
├── module/
│   ├── __init__.py
│   ├── msi_module.py              # 核心 MSI 领域模型
│   ├── msi_data_manager.py         # 通用数据管理器
│   └── msi_data_manager_zys.py     # MATLAB .mat 格式处理器
├── example.py                       # 使用示例
├── README.md                        # 英文文档
├── README_CN.md                     # 中文文档
└── LICENSE                          # GNU GPL v3 许可证
```

## 核心组件

### MSI 类
主要领域模型，封装了：
- 元数据管理
- MSI 切片队列
- 数据矩阵存储
- 查询和可视化方法

### MSIDataManager
通用数据管理器，支持：
- `.h5` 和 `.msi` 文件格式
- 从目录批量导入
- Split/Merge 导出模式
- m/z 范围过滤

### MSIDataManagerZYS
MATLAB `.mat` 文件的专用管理器，具有：
- 自定义数据结构处理
- 逐通道归一化
- 基于稀疏性的过滤
- HDF5 转换功能

## 开发与贡献指南

### 贡献流程与规范
- 贡献指南：参见 `docs/CONTRIBUTING.md`（中文）与 `docs/CONTRIBUTING_EN.md`（英文）
- Issue 模板：`.github/ISSUE_TEMPLATE/feature.md`、`.github/ISSUE_TEMPLATE/bug.md` 及英文模板 `feature_en.md`、`bug_en.md`（在 GitHub 的 `New issue` 页面自动可选）
- 命名规范：参见 `docs/NAMING_CONVENTIONS.md` 与 `docs/NAMING_CONVENTIONS_EN.md`
- 提交信息：推荐使用 Conventional Commits（如 `feat:`、`fix:`、`docs:`、`refactor:`、`test:`），示例：`feat(data-manager): support split/merge write modes`
- PR 检查清单要点：接口与命名一致、性能/内存无明显问题、断言与错误处理到位、文档/示例/测试同步更新、CI 通过、变更说明清晰

### 推荐的 Trae/VSCode 扩展（Python 与代码质量控制）

为提升代码质量和开发体验，建议安装以下扩展。工作区文件 `.vscode/extensions.json` 会自动推荐这些扩展。

- `ms-python.python` — Python 支持
- `ms-python.vscode-pylance` — 智能提示与类型分析
- `ms-python.pylint` — 静态检查（使用仓库 `.pylintrc`）
- `charliermarsh.ruff` — 快速代码规范与风格检查
- `ms-python.black-formatter` — 代码格式化（Black）
- `ms-python.isort` — 导入排序
- `h5web.vscode-h5web` — HDF5 可视化

命令行安装（可选）：

```bash
code --install-extension ms-python.python \
  && code --install-extension ms-python.vscode-pylance \
  && code --install-extension ms-python.pylint \
  && code --install-extension charliermarsh.ruff \
  && code --install-extension ms-python.black-formatter \
  && code --install-extension ms-python.isort \
  && code --install-extension h5web.vscode-h5web
```

## 许可证

本项目采用 GNU 通用公共许可证 v3.0 - 详见 [LICENSE](LICENSE) 文件。

## 参考资料

- [MSI 数据处理工作流程](https://pleinelune-r.github.io/2025/08/05/MSI%E6%95%B0%E6%8D%AE%E5%A4%84%E7%90%86%E6%B5%81%E7%A8%8B/)
- [Cardinal MSI](https://cardinalmsi.org/)
- [Cardinal GitHub 仓库](https://github.com/kuwisdelu/Cardinal/tree/devel/R)
- [MATLAB: 预处理原始质谱数据](https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html)
- [质谱成像预处理综述](https://www.sciencedirect.com/science/article/pii/S0169743921001015)
- [PyOpenMS 文档](https://pyopenms.readthedocs.io/en/latest/user_guide/background.html#why-use-openms)
- [MSI 分析最新进展](https://pubs.acs.org/doi/10.1021/jasms.4c00314)

## 贡献

欢迎贡献！请随时提交 Pull Request。

## 联系方式

如有问题和需要支持，请在 GitHub 仓库中提交 issue。
