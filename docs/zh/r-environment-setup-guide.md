# R 环境设置指南

MassFlow 的默认预处理功能由自身 Python 实现完成。同时MassFlow也提供 Cardinal 的接口。

只有在明确选择 `backend="cardinal"`，需要兼容已有 Cardinal 工作流或进行方法比较时，才需要安装并配置 R环境。

本文档介绍了如何配置 R 环境，包括设置 R 的安装路径。

[TOC]

## 1. R 环境检测机制

在初始化 R 环境时，MassFlow 会按照以下优先级顺序查找 R 的安装路径（`R_HOME`）:

1. **函数参数**：调用初始化函数时直接传递的路径。
2. **环境变量 `R_HOME`**：在操作系统环境变量中设置的 `R_HOME`。
3. **默认配置**：通过 `set_default_r_home` 设置的路径。
4. **自动检测**：尝试通过系统命令（`R RHOME`）自动获取 R 安装路径。

如果通过上述任何一种方法都无法找到有效的 R 路径，程序将抛出错误。

## 2. 如何配置 R 路径

### 方法 1：使用环境变量（推荐）

在运行 Python 脚本之前，请将系统环境变量 `R_HOME` 设置为指向您的 R 安装目录。

**Windows PowerShell 示例：**

```powershell
$env:R_HOME = "C:\Program Files\R\R-4.3.1"
python main.py
```

**Linux/macOS 示例:**

```bash
export R_HOME=/usr/lib/R
python main.py
```

### 方法 2：在代码中设置默认路径

可以在导入其他 MassFlow 模块之前，使用 `set_default_r_home` 函数来设置默认的 R 路径。

```python
from massflow.r_preprocess import set_default_r_home

# Set R installation path
set_default_r_home("C:\\Program Files\\R\\R-4.3.1")

# Proceed with other operations
from massflow.r_preprocess import init_r_environment
env = init_r_environment()
```

### 方法 3：在初始化时指定路径

在初始化 R 环境时，直接传入 `r_home` 参数。

```python
from massflow.r_preprocess import init_r_environment

# Initialize and specify path
env = init_r_environment(r_home="C:\\Program Files\\R\\R-4.3.1")
```

## 3. 常见问题

### `R_HOME 未找到` 错误

如果您遇到此错误，说明 MassFlow 无法自动找到 R。请检查 R 是否已安装，或尝试按照上文“如何配置 R 路径”中的方法之一手动指定路径。

### 与 `UnicodeDecodeError` 相关的错误

MassFlow 会自动将 R 的区域设置设为英语（`LANGUAGE=en`，`LC_ALL=C`），以避免字符编码问题。通常情况下，无需用户手动干预。





作者：董乐强，董京伟