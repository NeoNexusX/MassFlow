# R Environment Setup Guide

MassFlow's default preprocessing functions are implemented in Python by MassFlow itself. MassFlow also provides an interface to Cardinal.

You only need to install and configure an R environment when you explicitly select `backend="cardinal"`, need compatibility with an existing Cardinal workflow, or want to compare methods.

This page explains how to configure the R environment, including how to set the R installation path.

## 1. R Environment Detection

When the R environment is initialized, MassFlow searches for the R installation path (`R_HOME`) in the following order of priority:

1. **Function parameter**: a path passed directly when the initialization function is called.
2. **`R_HOME` environment variable**: the `R_HOME` value configured in the operating-system environment variables.
3. **Default configuration**: the path set through `set_default_r_home`.
4. **Automatic detection**: MassFlow attempts to obtain the R installation path with a system command (`R RHOME`).

If no valid R path can be found by any of these methods, the program raises an error.

## 2. How to Configure the R Path

### Method 1: Use an Environment Variable (Recommended)

Before running your Python script, set the system environment variable `R_HOME` to the R installation directory.

**Windows PowerShell example:**

```powershell
$env:R_HOME = "C:\Program Files\R\R-4.3.1"
python main.py
```

**Linux/macOS example:**

```bash
export R_HOME=/usr/lib/R
python main.py
```

### Method 2: Set the Default Path in Code

Before importing other MassFlow modules, use the `set_default_r_home` function to set the default R path.

```python
from massflow.r_preprocess import set_default_r_home

# Set R installation path
set_default_r_home("C:\\Program Files\\R\\R-4.3.1")

# Proceed with other operations
from massflow.r_preprocess import init_r_environment
env = init_r_environment()
```

### Method 3: Specify the Path During Initialization

Pass the `r_home` argument directly when initializing the R environment.

```python
from massflow.r_preprocess import init_r_environment

# Initialize and specify path
env = init_r_environment(r_home="C:\\Program Files\\R\\R-4.3.1")
```

## 3. Common Problems

### `R_HOME Not Found` Error

If you encounter this error, MassFlow was unable to locate R automatically. Check whether R is installed, or manually specify the path using one of the methods described above under “How to Configure the R Path.”

### Errors Related to `UnicodeDecodeError`

MassFlow automatically sets the R locale to English (`LANGUAGE=en`, `LC_ALL=C`) to avoid character-encoding problems. In most cases, no manual user intervention is required.





Authors: Leqiang Dong, Jingwei Dong
