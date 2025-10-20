# MassFlow

[English](README.md) | 简体中文

MassFlow 是一个模块化的高性能质谱成像（MSI）数据预处理框架。它为 MSI 研究提供高效的数据管理、处理和可视化功能。

## 功能特性

### ✅ 已实现功能

#### 数据导入和导出
- **支持格式**：`.h5`、`.msi`、`.mat`（MATLAB 格式）
- **批量处理**：导入整个目录的 MSI 文件
- **灵活的存储模式**：
  - Split 模式：每个 m/z 值保存为单独文件
  - Merge 模式：所有数据合并到单个文件
- **数据过滤**：加载特定 m/z 范围以减少内存使用
- **元数据管理**：自动提取和存储元数据

#### 数据管理
- **MSI 领域模型**：提供清晰的 API，包含元数据、切片队列和可选数据矩阵
- **内存高效加载**：选择性 m/z 范围加载以最小化内存占用
- **数据归一化**：逐通道最小-最大归一化
- **基础掩码生成**：自动生成组织区域掩码
- **稀疏性过滤**：基于稀疏性阈值过滤通道

#### 可视化
- **MSI 图像绘制**：显示指定 m/z 范围的图像
- **可自定义显示**：可调节阈值、色图和图形尺寸
- **批量可视化**：在指定范围内绘制多个图像
- **导出支持**：保存图像到文件或交互式显示

#### 专用处理器
- **MSIDataManager**：用于 `.h5`/`.msi` 文件的通用数据管理器
- **MSIDataManagerZYS**：用于 MATLAB `.mat` 格式的专用处理器，具有自定义预处理功能

### 📋 计划功能（TODO）

#### 数据导入
- imzML 格式支持
- mzML 格式支持
- Bruker `.d` 格式支持
- Agilent `.bd` 格式支持

#### 基线校正
例如：TopHat 滤波器

#### 去噪与平滑

#### 峰检测

#### 峰对齐

#### 归一化
- 总离子流（TIC）
- 中位数归一化
- 参考离子归一化
- 分箱（重采样和光谱分箱）

## 安装

```bash
# 克隆仓库
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow

# 安装依赖
pip install numpy h5py matplotlib pympler
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

### 示例 3：处理 MATLAB .mat 文件

```python
from module.msi_data_manager_zys import MSIDataManagerZYS

# Create MSI instance for .mat file
msi3 = MSI(name='20250329_sample', version=1.0, mask=None, need_base_mask=True)

# Initialize ZYS data manager
msi_dm_zys = MSIDataManagerZYS(msi3, filepath="./data/sample.mat")

# Load and rebuild data
msi_dm_zys.load_data_from_zys_mat()
msi_dm_zys.rebuild_hdf5_file_from_zys()

# Plot specific m/z range
msi3.plot_msi(target_mz_range=[100, 150])
```

### 示例 4：导出数据

```python
# Export as merged file
msi_dm_zys.write2local(mode='merge', output_fold='./data')

# Export as split files (one file per m/z)
msi_dm_zys.write2local(mode='split')
```

### 示例 5：查询特定 m/z 值

```python
# Get MSI slices for a specific m/z value
msi_list = msi3.get_msi_by_mz(mz_value_min=88.1122, tol=1e-3)
image = msi_list[0].msroi.T
base_mask = msi_list[0].base_mask.T
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
