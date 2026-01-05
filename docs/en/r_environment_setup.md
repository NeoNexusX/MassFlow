# R Environment Setup Guide

Some preprocessing features of MassFlow (such as Cardinal-based peak alignment) rely on the R language environment. To ensure these features work correctly, MassFlow provides a flexible R environment configuration mechanism.

This document describes how to configure the R environment, including setting the R installation path.

## 1. R Environment Detection Mechanism

When initializing the R environment, MassFlow looks for the R installation path (`R_HOME`) in the following order of priority:

1. **Function Argument**: The path passed directly when calling the initialization function.
2. **Environment Variable `R_HOME`**: The `R_HOME` set in the operating system environment variables.
3. **Default Configuration**: The path set via `set_default_r_home`.
4. **Automatic Detection**: Attempts to automatically obtain the R installation path via the system command (`R RHOME`).

If a valid R path cannot be found by any of the above methods, the program will raise an error.

## 2. How to Configure the R Path

### Method 1: Using Environment Variables (Recommended)

Before running the Python script, set the system environment variable `R_HOME` to point to your R installation directory.

**Windows PowerShell Example:**
```powershell
$env:R_HOME = "C:\Program Files\R\R-4.3.1"
python main.py
```

**Linux/macOS Example:**

```bash
export R_HOME=/usr/lib/R
python main.py
```

### Method 2: Setting Default Path in Code

You can use the `set_default_r_home` function to set the default R path before importing other MassFlow modules.

```python
from massflow.r_preprocess import set_default_r_home

# Set R installation path
set_default_r_home("C:\\Program Files\\R\\R-4.3.1")

# Proceed with other operations
from massflow.r_preprocess import init_r_environment
env = init_r_environment()
```

### Method 3: Specifying Path at Initialization

Pass the `r_home` argument directly when initializing the R environment.

```python
from massflow.r_preprocess import init_r_environment

# Initialize and specify path
env = init_r_environment(r_home="C:\\Program Files\\R\\R-4.3.1")
```

## 3. Common Issues

### `R_HOME not found` Error
If you encounter this error, it means MassFlow cannot automatically find R. Please check if R is installed, or try manually specifying the path using one of the methods in "How to Configure the R Path" above.

### `UnicodeDecodeError` Related Errors
MassFlow automatically sets the R locale to English (`LANGUAGE=en`, `LC_ALL=C`) to avoid character encoding issues. Usually, no manual user intervention is required.