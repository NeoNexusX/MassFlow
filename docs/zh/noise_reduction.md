# MassFlow 噪声抑制与滤波（filter.py）

本文档介绍 MassFlow 中的噪声抑制与滤波模块，重点说明 `preprocess/filter.py` 的函数及其与 `MSIPreprocessor.noise_reduction` 的协同使用。内容包含 API 说明、示例代码、参数及调参建议、适用场景与常见问题。
- 相关模块
  - `preprocess/filter.py`
  - `preprocess/ms_preprocess.py`（提供 `MSIPreprocessor.noise_reduction` 的统一入口）

## 概览

- 输入与输出
  - 输入：`module.ms_module.SpectrumBaseModule`（一维 `intensity` 必须存在；`mz_list` 可选）
  - 输出：一维平滑后的 `intensity` 或新的 `SpectrumBaseModule`（通过 `MSIPreprocessor.noise_reduction` 返回）
- 算法类别
  - 时域卷积平滑：`smooth_signal_ma`（移动平均/自定义核）、`smooth_signal_gaussian`（离散高斯）
  - 多项式拟合平滑：`smooth_signal_savgol`（Savitzky-Golay）
  - 小波去噪：`smooth_signal_wavelet`（基于 PyWavelets 的阈值化）
  - 基于 m/z 邻域搜索的平滑：`smooth_ns_signal_ma`、`smooth_ns_signal_gaussian`、`smooth_ns_signal_bi`（双边）
- 预处理
  - `smooth_preprocess`：将负数强度置零，清理数据引用，避免后续算法异常
  - 邻域搜索预处理（内部使用）：NaN 强度填充、`mz_list` 缺失时回退为等距索引

### 函数关系图

```mermaid
graph LR
  A[MSIPreprocessor.noise_reduction(data, method,...)] --> P[smooth_preprocess]
  A --> MA[smooth_signal_ma]
  A --> G[smooth_signal_gaussian]
  A --> SG[smooth_signal_savgol]
  A --> WT[smooth_signal_wavelet]
  A --> MA_NS[smooth_ns_signal_ma]
  A --> G_NS[smooth_ns_signal_gaussian]
  A --> BI_NS[smooth_ns_signal_bi]
```

## 核心 API

### MSIPreprocessor.noise_reduction

```python
preprocess.ms_preprocess.MSIPreprocessor.noise_reduction(
  data: SpectrumBaseModule | SpectrumImzML,
  method: str = "ma",
  window: int = 2,
  sd: float = 2,
  sd_intensity: float | None = None,
  p: int = 2,
  coef: np.ndarray | None = None,
  polyorder: int = 2,
  wavelet: str = "db4",
  threshold_mode: str = "soft"
) -> SpectrumBaseModule | SpectrumImzML
```

- 说明：统一的降噪入口。根据 `method` 分派到具体滤波实现，返回与输入同坐标的光谱对象，`intensity` 为平滑后的结果。
- 支持的 `method`：
  - `"ma"`、`"gaussian"`、`"savgol"`、`"wavelet"`
  - `"ma_ns"`、`"gaussian_ns"`、`"bi_ns"`
- 返回：新的 `SpectrumBaseModule`（或 `SpectrumImzML`）实例，`mz_list` 与坐标保留，`intensity` 替换为平滑结果
- 异常：`ValueError`（不支持的方法）、`TypeError`（输入类型不合法）

### smooth_signal_ma

```python
preprocess.filter.smooth_signal_ma(
  x: SpectrumBaseModule,
  coef: np.ndarray | None = None,
  window: int = 5
) -> np.ndarray
```

- 说明：移动平均（或自定义卷积核）平滑。一维边界使用 `edge` 填充，自动归一化权重。
- 参数：
  - `coef`：卷积核；为 `None` 时使用长度为 `window` 的均匀核
  - `window`：窗口长度（正整数；若为偶数自动调至奇数）
- 返回：与输入等长的一维 `intensity`
- 异常：`ValueError`（无 `coef` 时 `window <= 0`）、`TypeError`（`intensity` 非 1D）

### smooth_signal_gaussian

```python
preprocess.filter.smooth_signal_gaussian(
  x: SpectrumBaseModule,
  sd: float | None = None,
  window: int = 5
) -> np.ndarray
```

- 说明：离散高斯核平滑；默认 `sd = window / 4`
- 依赖：`scipy.stats.norm`（若不可用则退化为 `numpy` 实现）
- 参数：`sd`（高斯标准差）、`window`（奇数长度）
- 异常：`ValueError`（窗口非正）、`TypeError`（`intensity` 非 1D）

### smooth_signal_savgol

```python
preprocess.filter.smooth_signal_savgol(
  x: SpectrumBaseModule,
  window: int = 5,
  polyorder: int = 2
) -> np.ndarray
```

- 说明：Savitzky-Golay 多项式拟合平滑；自动保证 `window` 为奇数且 `window >= 3`
- 依赖：`scipy.signal.savgol_filter`
- 参数：`polyorder < window`（若不满足自动调低）
- 返回：与输入等长的一维 `intensity`

### smooth_signal_wavelet

```python
preprocess.filter.smooth_signal_wavelet(
  x: SpectrumBaseModule,
  wavelet: str = "db4",
  threshold_mode: str = "soft"
) -> np.ndarray
```

- 说明：小波阈值去噪；自动估计噪声并按 Donoho-Johnstone 阈值处理；重构长度与输入严格对齐（必要时截断/填充）
- 依赖：`pywt`（PyWavelets）
- 参数：`wavelet`（如 `db4`、`db8`、`haar`、`coif2`）、`threshold_mode`（`soft`/`hard`）

### 邻域搜索预处理（内部）

```python
preprocess.filter.smooth_ns_signal_pre(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 1
) -> tuple[np.ndarray, np.ndarray, np.ndarray]
```

- 说明：在 m/z 轴上做 kNN 搜索（`cKDTree`），返回邻域强度、距离与索引
- 行为：
  - `mz_list` 缺失或与强度长度不匹配时，回退为 `np.arange(N)`
  - `intensity` 中的 NaN 以邻域均值填充（邻域全 NaN 时使用全局中位数；仍 NaN 则置 0）
  - `k > N` 时自动将 `k` 调为 `N`
- 参数：`k >= 1`、`p >= 1`（Minkowski 距离参数；`p=1` 曼哈顿、`p=2` 欧氏）

### smooth_ns_signal_ma

```python
preprocess.filter.smooth_ns_signal_ma(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 1
) -> np.ndarray
```

- 说明：在 m/z 邻域上做等权平均平滑（行归一化）
- 参数：`k`（邻居数）、`p`（距离度量）

### smooth_ns_signal_gaussian

```python
preprocess.filter.smooth_ns_signal_gaussian(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 1,
  sd: float | None = None
) -> np.ndarray
```

- 说明：对邻域距离做高斯加权；指数在数值下溢前做 `clip`，行归一化避免除零
- 默认：`sd = median(max_row_distance) / 2`
- 参数：`k`、`p`、`sd`

### smooth_ns_signal_bi（双边）

```python
preprocess.filter.smooth_ns_signal_bi(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 2,
  sd_dist: float | None = None,
  sd_intensity: float | None = None
) -> np.ndarray
```

- 说明：同时考虑 m/z 距离与强度差的双边加权；适合保边（峰位）平滑
- 默认：
  - `sd_dist = median(max_row_distance) / 2`
  - `sd_intensity = scipy.stats.median_abs_deviation(intensity, scale="normal")`
- 数值稳定性：指数下溢前 `clip`，行归一化

### smooth_preprocess

```python
preprocess.filter.smooth_preprocess(data: SpectrumBaseModule) -> SpectrumBaseModule
```

- 说明：通用平滑前预处理，将负强度置零并清理引用，返回更新后的对象

## 快速上手示例

> 以下示例均为简化演示，图片位置留空。

### 1. 移动平均（MA）
```python
import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from module.ms_module import MS, SpectrumBaseModule
from module.ms_data_manager_imzml import MSDataManagerImzML
from preprocess.ms_preprocess import MSIPreprocessor
# Direct parameter settings
window_size = 7  # Window size for moving average
output_filename = "ma_denoised_plot.png"  # Output plot filename
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem

# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing (using directly set window size)
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="ma",
    window=window_size
)
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='MA'
)
# Save results
mz_den = np.asarray(denoised.mz_list, dtype=float)
inten_den = np.asarray(denoised.intensity, dtype=float)
out_base = f"{input_file_stem}_ma_denoised"
np.save(f"{out_base}_intensity.npy", inten_den)
np.save(f"{out_base}_mz.npy", mz_den)

print(f"Processing completed! Window size: {window_size}")
print(f"Plot saved as: {output_filename}")
print(f"Data saved as: {out_base}_intensity.npy and {out_base}_mz.npy")

```
![MA 滤波前后示例（待补图）]()
## 之后的示例代码，导入和保存部分与MA一致，将会被省略
### 2. 高斯平滑

```python
# Parameter settings
window_size = 7  # Window size for Gaussian filter
standard_deviation = 2.0  # Standard deviation for Gaussian filter
output_filename = "gaussian_denoised_plot.png"  # Output plot filename
# Check if window size is odd
if window_size % 2 == 0:
    print("Warning: Window size should be an odd number. Incrementing by 1.")
    window_size += 1
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Apply Gaussian noise reduction
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="gaussian",
    window=window_size,
    sd=standard_deviation
)
print(f"Denoising complete. Denoised spectrum length: {len(denoised.intensity)}")
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='seagreen',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Gaussian'
)
```

![高斯滤波前后示例（待补图）]()

### 3. Savitzky-Golay 平滑（保峰形）

```python
# Parameter settings
window_size = 11  # Window size for Savitzky-Golay filter (odd number, >=3)
polyorder = 3  # Polynomial order for Savitzky-Golay filter (must be less than window)
output_filename = "savgol_denoised_plot.png"  # Output plot filename
# Check window parameters
if window_size % 2 == 0:
    window_size += 1
    print(f"Warning: Window size must be odd. Adjusted to {window_size}.")
if window_size < 3:
    window_size = 3
    print(f"Warning: Window size must be at least 3. Adjusted to {window_size}.")
if polyorder >= window_size:
    polyorder = window_size - 1
    print(f"Warning: Polyorder must be less than window. Adjusted to {polyorder}.")
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Apply Savitzky-Golay noise reduction
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="savgol",
    window=window_size,
    polyorder=polyorder
)
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='darkorange',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Savitzky-Golay'
)

```

![SG 滤波前后示例（待补图）]()

### 4. 小波去噪

```python
# Parameter settings
wavelet_type = "db4"  # Wavelet type (e.g., 'db4', 'sym8', 'coif5')
threshold_mode = "hard"  # Thresholding mode ('soft' or 'hard')
output_filename = "wavelet_denoised_plot.png"  # Output plot filename
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]

print(f"Original spectrum loaded. Length: {len(sp.intensity)}")
print(f"Processing with Wavelet filter: wavelet={wavelet_type}, mode={threshold_mode}")
# Apply Wavelet noise reduction
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="wavelet",
    wavelet=wavelet_type,
    threshold_mode=threshold_mode
)
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='purple',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Wavelet'
)
```

![小波去噪前后示例（待补图）]()

### 5. 邻域搜索移动平均 (ma_ns)

```python
# Parameters
K = 7  # Number of neighbors
P = 1  # Minkowski distance parameter
OUTPUT_FILENAME = "ma_ns_denoised_plot.png"
FILE_PATH = "data\\neg-gz4.imzML"
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing
denoised_sp = MSIPreprocessor.noise_reduction(
    data=sp,
    method="ma_ns",
    window=K,  # For ns methods, window is used as k
    p=P
)
# Plotting
denoised_sp.plot(
    save_path=OUTPUT_FILENAME,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='MA-NS'
)
```

![MA_NS 滤波前后示例（待补图）]()

### 6. 邻域搜索高斯 (gaussian_ns)

```python
# Parameters
K = 7
P = 2
SD = None  # Let the function determine sd automatically
OUTPUT_FILENAME = "gaussian_ns_denoised_plot.png"
FILE_PATH = "data\\neg-gz4.imzML"
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing
denoised_sp = MSIPreprocessor.noise_reduction(
    data=sp,
    method="gaussian_ns",
    window=K,
    p=P,
    sd=SD
)
# Plotting
denoised_sp.plot(
    save_path=OUTPUT_FILENAME,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Gaussian-NS'
)
```
![Gaussian_NS 滤波前后示例（待补图）]()

### 7. 邻域搜索双边 (bi_ns)

```python
# Parameters
K = 7
P = 2
SD_DIST = None
SD_INTENSITY = None
OUTPUT_FILENAME = "bi_ns_denoised_plot.png"
FILE_PATH = "data\\neg-gz4.imzML"
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing
denoised_sp = MSIPreprocessor.noise_reduction(
    data=sp,
    method="bi_ns",
    window=K,
    p=P,
    sd=SD_DIST,
    sd_intensity=SD_INTENSITY
)
# Plotting
denoised_sp.plot(
    save_path=OUTPUT_FILENAME,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='BI-NS'
)
```

![Bilateral_NS 滤波前后示例（待补图）]()

### 6. 统一入口（推荐）：MSIPreprocessor.noise_reduction

```python
from preprocess.ms_preprocess import MSIPreprocessor

pre = MSIPreprocessor()
spec_ma = pre.noise_reduction(data=spec, method='ma', window=7)
spec_gauss = pre.noise_reduction(data=spec, method='gaussian', window=9, sd=2.5)
spec_sg = pre.noise_reduction(data=spec, method='savgol', window=9, polyorder=3)
spec_wt = pre.noise_reduction(data=spec, method='wavelet', wavelet='db4', threshold_mode='soft')

spec_ma_ns = pre.noise_reduction(data=spec, method='ma_ns', window=7, p=2)
spec_gauss_ns = pre.noise_reduction(data=spec, method='gaussian_ns', window=7, sd=None, p=2)
spec_bi_ns = pre.noise_reduction(data=spec, method='bi_ns', window=7, sd=None, sd_intensity=None, p=2)



## 参数说明与调参建议

- `window`（卷积/SG 窗口）
  - 正整数；偶数会自动转为最近的奇数（中心对齐）
  - 建议：5–15（根据峰宽与噪声强度调整）
- `coef`（自定义卷积核）
  - 传入时会自动归一化；长度决定有效窗口
- `sd`（高斯核标准差）
  - 经典卷积默认 `sd = window / 4`
  - 邻域高斯默认 `sd = median(max_row_distance) / 2`（数据自适应）
- `polyorder`（SG 多项式阶数）
  - 必须小于 `window`；建议 2–4
- `wavelet` 与 `threshold_mode`
  - 推荐：`db4` + `soft`；强噪声可试 `hard`
- `k`（邻域搜索近邻数）
  - `k >= 1`；过大易过平滑；一般 5–11
- `p`（Minkowski 距离参数）
  - `p=1` 曼哈顿、`p=2` 欧氏；更高的 `p` 更强调距离差异
- `sd_intensity`（双边强度差权重的尺度）
  - 默认为强度的 `median_abs_deviation`（鲁棒性好）

## 使用场景

- 基础降噪与快速平滑：`ma`、`gaussian`
- 保峰形与定量分析：`savgol`（适合尖锐峰且不希望峰位漂移）
- 强噪声/混合噪声：`wavelet`（多尺度鲁棒阈值）
- 不规则 m/z 采样或需要基于 m/z 距离的平滑：`ma_ns`、`gaussian_ns`
- 保边（峰位）且降噪：`bi_ns`（双边在尖峰附近更保真）

## 常见问题与排查

- 输入形状错误
  - 报错：`TypeError: x.intensity must be a 1D numpy array`
  - 解决：确保 `intensity` 为一维 `np.ndarray`
- 窗口参数无效
  - 报错：`ValueError: window must be a positive integer`
  - 解决：`window >= 3`，且最终会转为奇数
- SG 约束不满足
  - 情况：`polyorder >= window`；库会报错或过拟合
  - 解决：将 `polyorder` 调低（实现已自动校正为 `window-1`）
- 依赖缺失
  - 依赖：`numpy`、`scipy`（`stats`、`signal`、`spatial`）、`pywt`
  - 解决：按 `requirements.txt` 安装；或在 `gaussian` 缺 `scipy` 时用 `numpy` 实现
- 邻域搜索边界
  - 情况：`mz_list` 缺失或长度不匹配
  - 处理：自动回退 `mz_list = np.arange(N)`，并填充 `NaN` 强度
- 数值稳定性
  - 指数权重下溢：实现已 `clip` 到安全范围
  - 行归一化：避免除零；若行和为 0 则置为 1
- 统一入口方法值
  - 报错：`Unsupported smoothing method`
  - 解决：使用 `ma, gaussian, savgol, wavelet, ma_ns, gaussian_ns, bi_ns`

## 实践建议

- 预处理：先用 `smooth_preprocess` 将负强度清理
- 窗口选择：以目标峰宽为参考；过大窗口会拉平尖峰
- 自适应尺度：优先使用默认 `sd`/`sd_intensity`（基于数据鲁棒估计）
- 评估指标：对比峰高/峰宽/噪声底的变化，避免过平滑

## 参考

- `preprocess/filter.py`（滤波/去噪核心实现）
- `preprocess/ms_preprocess.py`（统一入口与参数分派）
- `module/ms_module.py`（`SpectrumBaseModule` 数据结构与可视化）