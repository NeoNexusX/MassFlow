# 质谱成像噪声消除算法使用指南
## 概述

`MSIPreprocessor.noise_reduction()` 方法现在支持四种去噪算法：

1. **移动平均** (`ma`) - 
2. **高斯平滑** (`gaussian`) -   
3. **Savitzky-Golay滤波** (`savgol`) 
4. **小波去噪** (`wavelet`) - 

### 使用方法

```python
from preprocess.ms_preprocess import MSIPreprocessor

# 基本使用
denoised_spectrum = MSIPreprocessor.noise_reduction(
    data=spectrum,
    method='savgol',
    window=11,        # 窗口大小（必须为奇数）
    polyorder=3       # 多项式阶数（必须小于窗口大小）
)
```

### 参数说明

- `window`: 窗口大小，必须为奇数且大于多项式阶数
- `polyorder`: 多项式阶数，通常取2-4，阶数越高保持细节越好但去噪效果可能降低

### 参数选择建议

- **窗口大小**: 
  - 小窗口(5-11): 保持更多细节，适合高分辨率数据
  - 大窗口(15-31): 更强的平滑效果，适合噪声较大的数据
  
- **多项式阶数**:
  - 2阶: 适合平滑的峰形
  - 3-4阶: 适合尖锐的峰形
  - 过高阶数可能引入振荡

## 性能比较

| 算法 | 计算速度 | 峰形保持 | 去噪效果 | 适用场景 |
|------|----------|----------|----------|----------|
| 移动平均 | 最快 | 一般 | 一般 | 快速预处理 |
| 高斯平滑 | 快 | 一般 | 好 | 一般去噪 |
| Savitzky-Golay | 中等 | 很好 | 好 | 保持峰形 |
| 小波去噪 | 较慢 | 好 | 很好 | 强噪声环境 |

---

## 四种算法参数选择指南

### 1) 移动平均（`ma`）
- 核心参数
  - `window`：窗口大小（建议使用奇数，内部会自动转为奇数）
  - `coef`：可选的自定义卷积核（不提供时为均匀权重）
- 选择建议
  - 峰宽与采样步长：若平均峰宽约为 `W` 个点，则可选 `window ≈ W/3 ~ W/2`，过大将导致峰形展宽或峰顶变平。
  - 数据分辨率高、峰密度大：用小窗口（5–11）；噪声大或峰稀疏：用大窗口（15–31）。
  - 需保留尖锐峰或肩峰：优先选择较小窗口，避免过度平滑。
- 快速配方
  - 常规：`window=11`
  - 高噪声：`window=21`
  - 保峰：`window=7`
- 边界与长度
  - 使用边缘值填充，输出长度与输入相同。

### 2) 高斯平滑（`gaussian`）
- 核心参数
  - `sd`：高斯核标准差（决定平滑强度）
  - `window`：核长度（内部通过移动平均实现卷积）
- 选择建议
  - FWHM（半高宽）关系：`FWHM = 2.355 × sd`
    - 若已知峰的典型 FWHM（单位：点数），可设 `sd ≈ FWHM / 2.355`
  - 采样步长已知（单位 m/z）：将 FWHM（m/z）除以步长得到点数，再套上式计算 `sd`。
  - `window` 尽量覆盖 ±3sd 区间：建议 `window ≈ int(6 × sd) + 1`（保证为奇数）
- 快速配方
  - 一般：`sd=2`，`window=11`
  - 更强平滑：`sd=3`，`window=19`
  - 保留细节：`sd=1.5`，`window=9`
- 注意
  - 高斯相比移动平均更“柔和”，对峰形展宽更小；但 `sd` 过大仍会平坦化峰顶。

### 3) Savitzky-Golay（`savgol`）
- 核心参数
  - `window`：窗口大小（必须为奇数且 > `polyorder`）
  - `polyorder`：多项式阶数（通常 2–4）
- 选择建议
  - 保峰优先：SG 滤波在去噪同时保持峰形与斜率，适合峰提取前的预处理。
  - 窗口与阶数的折中：
    - 低阶（`polyorder=2`）：更强平滑，边缘更稳健
    - 中高阶（`polyorder=3–4`）：保留更多细节，适合尖峰与复杂峰形
  - 典型组合：
    - 细节优先：`window=9–11`，`polyorder=3`
    - 更强平滑：`window=15–21`，`polyorder=2`
  - 过大窗口或过高阶数会引入振荡或伪峰，建议逐步增减观察。
- 快速配方
  - 通用：`window=11`，`polyorder=3`
  - 保峰：`window=9`，`polyorder=3`
  - 更强平滑：`window=17`，`polyorder=2`

### 4) 小波去噪（`wavelet`）
- 核心参数
  - `wavelet`：小波类型（如 `'db4'`、`'haar'`、`'coif2'`）
  - `threshold_mode`：阈值模式（`'soft'` 更平滑，`'hard'` 保留细节更强）
- 选择建议
  - 小波族选择：
    - `db4`/`db8`：通用型，适合大多数质谱数据（平滑且细节平衡）
    - `haar`：阶跃/突变信号；质谱中一般不推荐
    - `coif2`/`coif4`：更平滑的基函数，适合低噪或缓变峰形
  - 阈值模式：
    - `'soft'`：推荐默认，减少伪振铃与不连续
    - `'hard'`：保留更多高频细节，但可能引入毛刺
  - 阈值策略：目前实现使用通用（Donoho-Johnstone）阈值并对各细节层应用相同阈值，适合快速稳健去噪；如需更细粒度，可扩展为“按层自适应阈值”。
- 快速配方
  - 通用：`wavelet='db4'`，`threshold_mode='soft'`
  - 更强细节保留：`wavelet='db8'`，`threshold_mode='hard'`
  - 更平滑：`wavelet='coif2'`，`threshold_mode='soft'`
- 边界与长度
  - 使用 `'symmetric'` 边界条件；已在实现中处理重构后长度与输入一致（自动截断/填充）。

---

# 依赖与环境
- `ma`：仅依赖 `numpy`
- `gaussian`：需要 `scipy.stats`
- `savgol`：需要 `scipy.signal`
- `wavelet`：需要 `pywt`（PyWavelets）
- 安装示例（Windows）
  - `pip install numpy scipy pywavelets`

# 示例代码（不绘图，无指标）
```python
from pathlib import Path
import numpy as np
from module.ms_module import MS
from module.ms_data_manager_imzml import MSDataManagerImzML
from preprocess.ms_preprocess import MSIPreprocessor

file_path = str(Path("Dataset/neg-gz4.imzML"))
ms = MS()
ms_dm = MSDataManagerImzML(ms, filepath=file_path)
ms_dm.load_full_data_from_file()
spectrum = ms[0]

# 选择方法与参数（示例：移动平均）
denoised = MSIPreprocessor.noise_reduction(spectrum, method="ma", window=7)

# 保存结果
np.save("neg-gz4_ma_denoised_intensity.npy", np.array(denoised.intensity))
np.save("neg-gz4_ma_denoised_mz.npy", np.array(denoised.mz_list))
```

# 常见问题（FAQ）
- Q: 报错 `ImportError: No module named 'scipy'/'pywt'`？
  - A: 安装依赖或选择不需要依赖的 `ma` 方法。
- Q: `ValueError: Unsupported smoothing method`？
  - A: `method` 仅支持 `'ma'|'gaussian'|'savgol'|'wavelet'`。
- Q: `window` 是偶数或 `polyorder >= window` 怎么办？
  - A: 实现会自动将 `window` 调整为奇数，并将 `polyorder` 限制为 `< window`；但过高阶数可能导致振荡，建议手动检查。
- Q: 输出长度与输入不一致？
  - A: 当前实现对 MA/高斯/SG 保证长度一致；小波已做截断/填充处理匹配输入长度。
- Q: 去噪后峰顶变平或 SNR 降低？
  - A: 减小 `window` 或 `sd`，或尝试 `savgol`；对强噪声可用 `wavelet='db4', threshold_mode='soft'`。
- Q: 速度太慢？
  - A: 小波最慢；快速预览用 `ma` 或中等窗口的 `gaussian`/`savgol`。

---

# Guide to Noise Reduction Algorithms for Mass Spectrometry Imaging

## Overview

`MSIPreprocessor.noise_reduction()` supports four denoising methods:

1. Moving Average (`ma`)
2. Gaussian Smoothing (`gaussian`)
3. Savitzky–Golay (`savgol`)
4. Wavelet Denoising (`wavelet`)

### Usage

```python
from preprocess.ms_preprocess import MSIPreprocessor

denoised_spectrum = MSIPreprocessor.noise_reduction(
    data=spectrum,
    method='savgol',
    window=11,        # window size (must be odd)
    polyorder=3       # polynomial order (must be < window)
)
```

### Parameters

- `window`: Window size (must be odd and greater than `polyorder`)
- `polyorder`: Polynomial order (typically 2–4); higher orders preserve more details but may reduce smoothing effectiveness

### Parameter Selection Tips

- Window size:
  - Small (5–11): preserves more details; suitable for high-resolution data
  - Large (15–31): stronger smoothing; suitable for high noise
- Polynomial order:
  - Order 2: suitable for smooth peak shapes
  - Order 3–4: suitable for sharp or complex peaks
  - Excessive order may introduce oscillations

## Performance Comparison

| Method           | Speed | Peak Preservation | Denoising | Use Case           |
|------------------|-------|-------------------|-----------|--------------------|
| Moving Average   | Fastest | Average         | Average   | Quick preprocessing |
| Gaussian         | Fast  | Average           | Good      | General denoising   |
| Savitzky–Golay   | Medium | Very good        | Good      | Peak preservation   |
| Wavelet          | Slow  | Good              | Very good | Strong noise        |

---

## Parameter Selection Guide for Each Method

### 1) Moving Average (`ma`)
- Key parameters
  - `window`: window size (use odd; auto-adjusted to odd internally)
  - `coef`: optional custom convolution kernel (uniform weights if not provided)
- Tips
  - Peak width vs. sampling: if average peak width is `W` points, use `window ≈ W/3 ~ W/2`; overly large window may flatten peaks.
  - High resolution, dense peaks: small window (5–11); high noise or sparse peaks: large window (15–31).
  - Preserve sharp peaks or shoulders: prefer smaller window to avoid over-smoothing.
- Quick recipes
  - General: `window=11`
  - High noise: `window=21`
  - Peak preservation: `window=7`
- Boundary & length
  - Edge padding; output length equals input length.

### 2) Gaussian Smoothing (`gaussian`)
- Key parameters
  - `sd`: standard deviation of Gaussian kernel (controls smoothing strength)
  - `window`: kernel length (implemented via moving average convolution internally)
- Tips
  - FWHM relation: `FWHM = 2.355 × sd`
    - If typical peak FWHM (in points) is known, set `sd ≈ FWHM / 2.355`.
  - If sampling step (in m/z) is known: convert FWHM in m/z to points by dividing by step, then compute `sd`.
  - Cover ±3sd: use `window ≈ int(6 × sd) + 1` (odd).
- Quick recipes
  - General: `sd=2`, `window=11`
  - Strong smoothing: `sd=3`, `window=19`
  - Detail preservation: `sd=1.5`, `window=9`
- Notes
  - Gaussian is smoother than MA and causes less peak broadening; very large `sd` can still flatten peaks.

### 3) Savitzky–Golay (`savgol`)
- Key parameters
  - `window`: window size (must be odd and > `polyorder`)
  - `polyorder`: polynomial order (usually 2–4)
- Tips
  - Peak preservation: SG preserves peak shapes and slopes; suitable before peak picking.
  - Tradeoff between window and order:
    - Low order (`polyorder=2`): stronger smoothing, robust at boundaries
    - Medium–high order (`polyorder=3–4`): preserves more details; good for sharp/complex peaks
  - Typical combos:
    - Detail-first: `window=9–11`, `polyorder=3`
    - Strong smoothing: `window=15–21`, `polyorder=2`
  - Large window or high order may introduce oscillations or artifacts; adjust gradually.
- Quick recipes
  - General: `window=11`, `polyorder=3`
  - Peak preservation: `window=9`, `polyorder=3`
  - Strong smoothing: `window=17`, `polyorder=2`

### 4) Wavelet Denoising (`wavelet`)
- Key parameters
  - `wavelet`: wavelet family (e.g., `'db4'`, `'haar'`, `'coif2'`)
  - `threshold_mode`: `'soft'` (smoother) or `'hard'` (more details)
- Tips
  - Choosing wavelets:
    - `db4`/`db8`: general-purpose; balanced smoothing and details
    - `haar`: step-like signals; generally not recommended for MS
    - `coif2`/`coif4`: smoother basis, good for low-noise or slowly varying peaks
  - Threshold mode:
    - `'soft'`: recommended; fewer ringing/discontinuities
    - `'hard'`: retains more high-frequency details; may introduce spikes
  - Threshold strategy: current implementation uses a common Donoho–Johnstone threshold applied to all detail levels; good for quick robust denoising. Layer-wise adaptive thresholding can be added if needed.
- Quick recipes
  - General: `wavelet='db4'`, `threshold_mode='soft'`
  - More detail preservation: `wavelet='db8'`, `threshold_mode='hard'`
  - Smoother: `wavelet='coif2'`, `threshold_mode='soft'`
- Boundary & length
  - `'symmetric'` boundary; reconstruction length matches input via trim/pad.

---

# Dependencies & Environment
- `ma`: numpy only
- `gaussian`: scipy.stats
- `savgol`: scipy.signal
- `wavelet`: pywt (PyWavelets)
- Install (Windows):
  - `pip install numpy scipy pywavelets`

# Example (no plotting, no metrics)
```python
from pathlib import Path
import numpy as np
from module.ms_module import MS
from module.ms_data_manager_imzml import MSDataManagerImzML
from preprocess.ms_preprocess import MSIPreprocessor

file_path = str(Path("Dataset/neg-gz4.imzML"))
ms = MS()
ms_dm = MSDataManagerImzML(ms, filepath=file_path)
ms_dm.load_full_data_from_file()
spectrum = ms[0]

denoised = MSIPreprocessor.noise_reduction(spectrum, method="ma", window=7)

np.save("neg-gz4_ma_denoised_intensity.npy", np.array(denoised.intensity))
np.save("neg-gz4_ma_denoised_mz.npy", np.array(denoised.mz_list))
```

# FAQ
- Q: `ImportError: No module named 'scipy'/'pywt'`?
  - A: Install dependencies or choose `ma`, which requires only numpy.
- Q: `ValueError: Unsupported smoothing method`?
  - A: `method` supports only `'ma'|'gaussian'|'savgol'|'wavelet'`.
- Q: What if `window` is even or `polyorder >= window`?
  - A: Implementation auto-adjusts `window` to odd and clamps `polyorder < window`; however, high orders may cause oscillations—tune manually.
- Q: Output length mismatch?
  - A: MA/Gaussian/SG keep length equal; Wavelet handles trim/pad to match input.
- Q: Flattened peaks or lower SNR after denoising?
  - A: Reduce `window` or `sd`, or try `savgol`; for heavy noise use `wavelet='db4', threshold_mode='soft'`.
- Q: Too slow?
  - A: Wavelet is the slowest; use `ma` or moderate Gaussian/SG for quick preview.