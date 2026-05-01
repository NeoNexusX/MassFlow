# MassFlow

本文介绍 MassFlow 中的归一化（normalization）模块，重点包括：
- 流水线级入口 `Preprocessor.normalization`
- `preprocess/helper/normalizer_helper.py` 中的核心函数

## 概述

- 输入与输出
  - 输入为扁平数组形式（`mz_data`、`intensity`、`lengths`），表示一批谱图的 m/z 和强度数据。
  - 输出为归一化后的扁平数组，m/z 与坐标保持不变，强度替换为归一化后的值。
- 方法
  - TIC（Total Ion Current）归一化（`tic_numba`）：将每条谱的强度缩放使总和（sum）等于指定值。
  - RMS（Root Mean Square）归一化（`rms_numba`）：将每条谱的强度缩放使均方根（RMS）等于指定值。
  - 参考峰归一化（`ref_numba`）：将每条谱的强度缩放使指定 m/z 附近的信号归一化到目标值。
- 所有方法均在 Numba JIT 内核中按谱执行，利用 `prange` 并行循环。

### 函数关系示意图

```mermaid
graph LR
  A["Preprocessor.normalization"] --> B["FlatPreprocess.normalization_flat"]
  B --> C["normalizer"]
  C --> D["normalizer_numba<br/>(Numba JIT)"]
  D --> E["归一化后数据"]
```

## 核心 API

### Preprocessor.normalization（流水线级）

```python
massflow.preprocess.api.PreprocessorAPI.normalization(
    *,
    method: str = "tic_numba",
    scale: float | None = None,
    mz_flat: np.ndarray | None = None,
    ref: float | None = None,
    ref_tolerance: float = 0.1,
) -> Self
```

- 描述：流水线级归一化入口。通过 `Preprocessor` 链式调用，按批次对所有光谱执行归一化。
- 输入：通过 `Preprocessor(dm)` 绑定数据管理器。
- 输出：调用 `.start()` 后返回新的 `MSDataManagerImzML`，其光谱的强度为归一化后的结果。
- 支持方法：`"tic_numba"`、`"rms_numba"`、`"ref_numba"`。
- 参数：
  - `method`：归一化方法，默认 `"tic_numba"`。
  - `scale`：归一化后的幅度缩放因子（可选）。若提供，强度会被乘以此值。
  - `mz_flat`：扁平 m/z 数组（`ref_numba` 方法必需）。
  - `ref`：参考 m/z 值（`ref_numba` 方法必需）。
  - `ref_tolerance`：参考峰匹配容差（`ref_numba` 方法使用），默认 0.1。
- 异常：
  - `ValueError`：不支持的 method；scale 非法；ref_numba 缺少必需参数。

示例：

```python
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor
from massflow.tools.plot import plot_spectrum

FILE_PATH = "data/example.imzML"
dm = MSDataManagerImzML(filepath=FILE_PATH)
dm.load_head_data()

normalized = (
    Preprocessor(dm)
    .noise_reduction(method="gaussian_numba", window=7)
    .normalization(method="tic_numba")
    .start()
)

plot_spectrum(processed.ms[0])
normalized.close()
dm.close()
```

![image-20260430140609289](https://neonexus-picture.oss-ap-southeast-1.aliyuncs.com/test/image-20260430140609289.png)

### normalizer（辅助函数级）

```python
preprocess.helper.normalizer_helper.normalizer(
    intensity: np.ndarray,
    method: str = "tic_numba",
    scale: float | None = None,
    mz_flat: np.ndarray | None = None,
    ref: float | None = None,
    ref_tolerance: float = 0.1,
    lengths: np.ndarray | None = None,
) -> np.ndarray
```

- 描述：扁平数组归一化分发器。根据 `method` 调用对应的 Numba 后端。
- 参数：
  - `intensity`：一维扁平强度数组。
  - `lengths`：光谱长度数组，用于在 Numba 内部按谱拆分处理。
  - 其余参数与流水线 API 一致。
- 返回值：归一化后的强度数组。
- 支持的方法仅限 `"tic_numba"`、`"rms_numba"` 和 `"ref_numba"`。

### FlatPreprocess.normalization_flat（扁平批处理级）

```python
massflow.preprocess.flat_pre_fun.FlatPreprocess.normalization_flat(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    method: str = "tic_numba",
    scale: float | None = None,
    mz_flat: np.ndarray | None = None,
    ref: float | None = None,
    ref_tolerance: float = 0.1,
) -> FlatBatchResult
```

- 描述：`Preprocessor.normalization` 的底层实现。返回包含 mz_data、归一化后 intensity 和 lengths 的 `FlatBatchResult`。

## 方法详解

### TIC 归一化（tic_numba）

每条谱的强度除以该谱的强度总和，使 sum = 1。若提供了 `scale`，结果再乘以 scale。

适用场景：需要像素间总强度可比，适合可视化和相对定量分析。

### RMS 归一化（rms_numba）

每条谱的强度除以该谱的均方根（`sqrt(mean(x^2))`），使 RMS = 1。若提供了 `scale`，结果再乘以 scale。

适用场景：对强峰不敏感，在含高动态范围信号的数据中更稳定。

### 参考峰归一化（ref_numba）

在每条谱中查找容差窗口内距离 `ref` 最近的峰，用该峰的强度值作为归一化因子。若提供了 `scale`，结果再乘以 scale。

适用场景：已知某内标 m/z，希望相对于该内标归一化。需要同时提供 `mz_flat` 和 `ref` 参数。

## 参数说明与调参建议

- `method` 选择
  - 需要像素间总强度可比时，优先选 `"tic_numba"`。
  - 希望对极端值更稳健时，选 `"rms_numba"`。
  - 已知内标峰时，选 `"ref_numba"`。
- `scale`（幅度缩放因子）
  - 归一化后的额外缩放。常用场景：TIC 归一化后设为 1000 或 10000 便于数值展示。
  - 必须为有限非负数。
- `ref` + `mz_flat`（参考峰归一化）
  - `ref`：目标内标峰 m/z 值。
  - `mz_flat`：与 intensity 同长度的 m/z 扁平数组。
  - `ref_tolerance`：内标峰搜索容差（绝对差值），默认 0.1 Da。

## 参考

- `preprocess/helper/normalizer_helper.py`：归一化分发器（normalizer）
- `preprocess/numba/normalization_numba.py`：Numba 加速实现（normalizer_numba）
- `preprocess/flat_pre_fun.py`：扁平预处理入口（FlatPreprocess.normalization_flat）
- `preprocess/api.py`：链式 API（PreprocessorAPI.normalization）

## 错误处理与日志

- 输入校验：intensity 必须为非空一维 numpy 数组，否则抛出 `TypeError`。
- scale 校验：必须为有限非负数。
- ref_tolerance 校验：必须为有限非负数。
- method 校验：仅支持 `"tic_numba"`、`"rms_numba"`、`"ref_numba"`。
- ref_numba 缺参校验：缺少 `mz_flat` 或 `ref` 时抛出 `ValueError`。
