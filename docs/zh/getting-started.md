---
title: 快速开始
---

# 快速开始

MassFlow 是一个面向质谱成像（MSI）与质谱（MS）数据的模块化框架，提供高效的数据读取、预处理（去噪/平滑）与基础数据管理能力。

## 前置条件
- Python `>= 3.9`（开发测试版本为 `3.12`）

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

## 设置 Python 环境（macOS）
建议使用虚拟环境：

```bash
# 在仓库根目录
python3 -m venv .venv
source .venv/bin/activate  # macOS/Linux

# 安装依赖
pip install -r requirements.txt
```

如果使用 `conda`：
```bash
conda create -n massflow python=3.12 -y
conda activate massflow
pip install -r requirements.txt
```

## 运行示例
执行示例脚本验证环境：

```bash
python example.py
```

该脚本演示 MSI/MS 数据读取与基础预处理。日志保存在 `logs/`。

## 常用命令
- 构建静态站点：`npm run docs:build`
- 预览已构建站点：`npm run docs:preview`

## 故障排查
- 确认 Python 版本（`python3 --version`）以及虚拟环境激活状态（`which python`）。
- 若 `pip` 受网络或 SSL 影响，可尝试 `pip install -r requirements.txt --no-cache-dir` 或设置可信镜像。
- Apple Silicon（M1/M2/M3）建议使用 Python 3.11+ 与原生轮子；先升级 pip：`python -m pip install --upgrade pip`。

## 下一步
- 贡献指南：`docs/CONTRIBUTING.md` / `docs/CONTRIBUTING_EN.md`
- 命名规范：`docs/NAMING_CONVENTIONS.md` / `docs/NAMING_CONVENTIONS_EN.md`
- 协作指引：`docs/协作指北.md` / `docs/Collaboration_Guide.md`