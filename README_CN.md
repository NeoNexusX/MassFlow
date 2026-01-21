# MassFlow

[English](README.md) | 简体中文

MassFlow 是一个用于质谱成像 (MSI) 和质谱 (MS) 数据的模块化预处理和数据管理框架。

## 获取代码

```bash
# 克隆仓库
git clone https://github.com/NeoNexusX/MassFlow.git
cd MassFlow
```

## 在线文档

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
│       ├── module/               # 核心数据模型和数据管理器
│       │   ├── spectrum.py       # 光谱基础类型 (Spectrum)
│       │   ├── spectrum_imzml.py # 带惰性加载的 ImzML 光谱 (SpectrumImzML)
│       │   ├── mass_spectrum_set.py  # 光谱集合容器 (MassSpectrumSet)
│       │   ├── ms_meta_data.py       # 元数据结构 (ImzMlMetaData 等)
│       │   └── ...
│       ├── preprocess/              # 预处理入口与辅助函数
│       │   ├── spectrum_pre_fun.py  # 光谱级预处理 API (SpectrumPreprocess)
│       │   ├── dm_pre_fun.py        # 数据管理器级预处理 API (Preprocess)
│       │   ├── batch_pre_fun.py     # 批处理工具 (BatchPreprocess)
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

### 快速开发入门

建议直接运行以下代码 (main.py) 来验证数据加载：

```bash
# 首先关闭 conda：
conda deactivate

# 安装 uv：
# 请参考：https://docs.astral.sh/uv/getting-started/installation/

# 例如：
# Linux && macOS:
curl -LsSf https://astral.sh/uv/install.sh | sh 
# Windows:
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# 安装依赖：
uv sync 
uv pip install -e .

# 运行代码
uv run python main.py
```

- 请先阅读贡献指南：`docs/en/contribution.md` 和 `docs/zh/contribution.md`
- 命名规范：`docs/en/naming-conventions.md` 和 `docs/zh/naming-conventions.md`
- 提交规范：Conventional Commits (例如：`feat:`, `fix:`, `docs:`, `refactor:`, `test:`)
- 推荐扩展：Python, Pylance, Pylint, H5Web

## 许可证

本项目采用 GNU General Public License v3.0 许可证 - 详情请参阅 [LICENSE](LICENSE) 文件。

## 参考资料

- MATLAB Mass Spectrometry Preprocessing: https://www.mathworks.com/help/bioinfo/ug/preprocessing-raw-mass-spectrometry-data.html
- Cardinal MSI: https://cardinalmsi.org/
  - Cardinal 指南：https://bioconductor.org/packages/devel/bioc/vignettes/Cardinal/inst/doc/Cardinal3-guide.html
  - Cardinal Github：https://github.com/kuwisdelu/Cardinal
  - Matter Github：https://github.com/kuwisdelu/matter
- PyOpenMS: https://pyopenms.readthedocs.io/
  - github：https://github.com/OpenMS/OpenMS
  - docs：https://github.com/OpenMS/OpenMS-docs
  - 贡献库：https://github.com/OpenMS/contrib
  - Flash APP：https://github.com/OpenMS/FLASHApp

## 反馈

如需支持或报告问题，请提交 GitHub Issue。
