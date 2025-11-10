# MassFlow 基线校正（baseline_correction.py）

本文档介绍 MassFlow 中的基线校正模块，重点说明 `preprocess/baseline_correction.py` 的核心接口及其在 `MSIPreprocessor.baseline_correction` 中的使用。内容包含 API 说明、算法原理、参数与调参建议、示例代码（来自 `test_Asls.py` 与 `test_snip.py`）、适用场景与常见问题。

- 相关模块
  - `preprocess/baseline_correction.py`（ASLS，SNIP 实现与统一入口）
  - `module/ms_module.py`（光谱数据结构）

## 概览

- 输入与输出
  - 输入：一维 `intensity` 数组（`np.ndarray`）或 `module.ms_module.SpectrumBaseModule`
  - 输出：二元组 `(corrected, baseline)`：
    - 当输入为 `np.ndarray`：返回 `(np.ndarray corrected, np.ndarray baseline)`
    - 当输入为 `SpectrumBaseModule`：返回 `(SpectrumBaseModule corrected_spectrum, np.ndarray baseline)`；`mz_list` 与坐标保留，`intensity` 替换为基线校正后的结果
- 算法类别
  - ASLS（Asymmetric Least Squares）：非对称最小二乘，鲁棒基线估计，保峰形
  - SNIP（Statistics-Sensitive Non-linear Iterative Peak-clipping）：统计敏感非线性迭代截峰，自适应早停，避免过度扣底
- 预防过校正
  - `baseline_scale`：对估计的基线乘一个系数（0–1），用于保留一定背景以避免过度扣底
  - 若希望忠实还原算法原始行为（不缩放基线），可设置 `baseline_scale=1.0`

### 函数关系图

```mermaid
graph LR
  A[MSIPreprocessor.baseline_correction(data, method, ...)] --> ASLS[asls_baseline(y, lam, p, niter)]
  A --> SNIP[snip_baseline(y, m, decreasing, epsilon)]
  ASLS --> B[estimated baseline]
  SNIP --> B
  B --> C[scaled_baseline = baseline_scale * baseline]
  C --> D[corrected = y - scaled_baseline; clip to non-negative]
```

## 核心 API

### MSIPreprocessor.baseline_correction

```python
preprocess.baseline_correction.MSIPreprocessor.baseline_correction(
  data: np.ndarray | SpectrumBaseModule,
  method: str = "asls",
  lam: float = 1e7,
  p: float = 0.01,
  niter: int = 15,
  baseline_scale: float = 0.8,
  m: int | None = None,
  decreasing: bool = True,
  epsilon: float = 1e-3
) -> tuple[np.ndarray | SpectrumBaseModule, np.ndarray]
```

- 说明：统一的基线校正入口。根据 `method` 分派到 ASLS 或 SNIP，并返回修正后的信号与估计基线。
- 支持的 `method`：
  - `"asls"`：建议用于光滑背景扣除，保留峰形
  - `"snip"`：迭代截峰，适合快速稳健扣底，带自适应早停
- 返回：二元组 `(corrected, baseline)`，类型根据输入决定（见“概览”）
- 异常：
  - `ValueError`（不支持的方法）
  - `TypeError`（输入类型不合法）

#### ASLS（内部实现）

- 目标：通过非对称权重抑制“高于基线”的峰值，使二阶差分平滑项主导基线形状
- 关键参数：
  - `lam`（平滑项权重）：建议范围 `1e4–1e8`；越大越平滑
  - `p`（非对称权重）：建议范围 `0.001–0.1`；越小越保守（更保峰）
  - `niter`（迭代次数）：典型 `5–30`；迭代更新权重
- 数值实现：
  - 优先使用 `scipy.sparse` 构造二阶差分矩阵 `D` 并解稀疏线性方程 `(W + λD'D)z = Wy`
  - 缺少 `scipy` 时退化为致密矩阵实现（速度较慢）

#### SNIP（内部实现）

- 目标：迭代地将信号与邻域平均进行“下截”，逐步抑制峰值，仅保留背景
- 关键参数：
  - `m`（最大半窗）：默认自动 `min(50, n//10)`；表示截峰邻域半宽
  - `decreasing`（迭代顺序）：`True` 表示从粗到细（`p=m..1`），更常用；`False` 表示从细到粗（`1..m`）
  - `epsilon`（相对变化早停阈值）：`1e-5–1e-3`；越大越早停，避免过度扣底
- 早停策略：
  - 在每个窗长迭代中，比较更新前后片段的相对变化，若低于 `epsilon` 则提前终止迭代

## 快速上手示例

### 1. ASLS 示例（使用明显基线漂移数据）

```python
import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Allow importing from project root
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from module.ms_module import SpectrumBaseModule
from preprocess.baseline_correction import MSIPreprocessor

def compute_metrics(mz: np.ndarray, y: np.ndarray, corrected: np.ndarray, baseline: np.ndarray):
    # Negative value ratio
    negative_ratio = float(np.mean(corrected < 0)) if corrected.size > 0 else 0.0
    # Correlation between original and corrected
    corr_all = float(np.corrcoef(y, corrected)[0, 1]) if y.size == corrected.size else float("nan")
    # Baseline smoothness (RMS of second derivative)
    d2 = np.diff(np.asarray(baseline, dtype=float), n=2)
    baseline_smoothness = float(np.sqrt(np.mean(d2 * d2)))
    # TIC retention ratio (consistent with test_baseline_visualization)
    tic_ratio = float(np.sum(corrected) / max(np.sum(y), 1e-12))
    return {
        "negative_ratio": negative_ratio,
        "corr_all": corr_all,
        "baseline_smoothness": baseline_smoothness,
        "tic_ratio": tic_ratio,
    }
# --- Parameters (ASLS) ---
lam = 5e7        # Smoothness weight, larger = smoother (reference baseline_visualization default)
p = 0.008        # Asymmetry factor, smaller = more conservative
niter = 25       # Iteration count
baseline_scale = 0.9  # Baseline scaling to avoid over-correction

# --- File and Output Settings ---
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
FILE_PATH = os.path.join(DATA_DIR, "gaussian_denoised_intensity.npy")  # Consistent with test_snip
DENOISED_MZ_FILE = os.path.join(DATA_DIR, "gaussian_denoised_mz.npy")
OUTPUT_FILENAME = "asls_baseline_correction_visualization_400-500.png"

# --- Data Loading ---
print(f"Loading data from {FILE_PATH}...")
intensity_data = np.load(FILE_PATH)
mz_data = np.load(DENOISED_MZ_FILE)

# Use first spectrum if multiple spectra are included
intensity = intensity_data[0] if intensity_data.ndim > 1 else intensity_data

# Create SpectrumBaseModule
sp = SpectrumBaseModule(mz_list=mz_data, intensity=intensity, coordinates=[0, 0, 0])
print(f"Loaded spectrum with {len(sp.mz_list)} data points.")
print("Data loaded successfully.")

# --- Baseline Correction (ASLS) ---
print(f"Applying ASLS baseline correction with lam={lam:.2e}, p={p:.4f}, niter={niter}, baseline_scale={baseline_scale:.2f}...")
corrected_sp, estimated_baseline = MSIPreprocessor.baseline_correction(
    data=sp,
    method="asls",
    lam=lam,
    p=p,
    niter=niter,
    baseline_scale=baseline_scale
)
# --- Metrics ---
metrics = compute_metrics(sp.mz_list, sp.intensity, corrected_sp.intensity, estimated_baseline)
print("\n=== Baseline Correction Metrics (ASLS) ===")
print(f"  Negative Ratio: {metrics['negative_ratio']:.4f}")
print(f"  Correlation with Original: {metrics['corr_all']:.4f}")
print(f"  Baseline Smoothness (RMS of 2nd derivative): {metrics['baseline_smoothness']:.6f}")
print(f"  TIC Retention Ratio: {metrics['tic_ratio']:.4f}")
print("========================================\n")

# --- Plotting ---
plt.figure(figsize=(16, 5))
plt.plot(sp.mz_list, sp.intensity, label='Original', color='steelblue', linewidth=1.0, alpha=0.9); plt.plot(sp.mz_list, estimated_baseline, label='Estimated Baseline', color='forestgreen', linewidth=1.0, alpha=0.9); plt.plot(corrected_sp.mz_list, corrected_sp.intensity, label='Corrected', color='darkorange', linewidth=1.0, alpha=0.9)
plt.title('ASLS Baseline Correction Visualization'); plt.xlabel("m/z"); plt.ylabel("Intensity"); plt.xlim(400, 450); plt.ylim(0, 2); plt.legend(); plt.grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig(OUTPUT_FILENAME, dpi=300); plt.close()
# --- Save Corrected Data ---
input_file_stem = Path(FILE_PATH).stem
out_base = f"{input_file_stem}_asls_corrected"
np.save(f"{out_base}_intensity.npy", np.asarray(corrected_sp.intensity, dtype=float))
np.save(f"{out_base}_mz.npy", np.asarray(corrected_sp.mz_list, dtype=float))
print(f"Corrected data saved to {out_base}_intensity.npy and {out_base}_mz.npy")
```

![ASLS 前后示例（待补图）]()

### 2. SNIP 示例（使用正常数据）

```python
import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from module.ms_module import SpectrumBaseModule
from preprocess.baseline_correction import MSIPreprocessor


def compute_metrics(mz: np.ndarray, y: np.ndarray, corrected: np.ndarray, baseline: np.ndarray):
    # Negative value ratio
    negative_ratio = float(np.mean(corrected < 0)) if corrected.size > 0 else 0.0
    # Correlation between original and corrected
    corr_all = float(np.corrcoef(y, corrected)[0, 1]) if y.size == corrected.size else float("nan")
    # Baseline smoothness (RMS of second derivative)
    d2 = np.diff(np.asarray(baseline, dtype=float), n=2)
    baseline_smoothness = float(np.sqrt(np.mean(d2 * d2)))
    # TIC retention ratio
    tic_y = np.sum(np.maximum(0, y))
    tic_corrected = np.sum(np.maximum(0, corrected))
    tic_ratio = float(tic_corrected / max(tic_y, 1e-12))
    return {
        "negative_ratio": negative_ratio,
        "corr_all": corr_all,
        "baseline_smoothness": baseline_smoothness,
        "tic_ratio": tic_ratio,
    }

# --- Parameters ---
# SNIP algorithm parameters
m = 30  # Maximum half-window for the filter. If None, it will be auto-determined.
decreasing = True  # Use decreasing window sizes (from m to 1), which is generally more effective.
epsilon = 1e-5  # Threshold for adaptive early stopping to prevent over-correction.
baseline_scale = 0.9  # Scale factor for the estimated baseline; 1.0 keeps original SNIP behavior.

# --- File and Output Settings ---
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
FILE_PATH = os.path.join(DATA_DIR, "neg-gz4_savgol_denoised_intensity.npy") # Use the denoised file as the main input
DENOISED_MZ_FILE = os.path.join(DATA_DIR, "neg-gz4_savgol_denoised_mz.npy")
OUTPUT_FILENAME = "snip_baseline_correction_visualization_400-500.png"
# --- Data Loading ---
print(f"Loading data from {FILE_PATH}...")
intensity_data = np.load(FILE_PATH)
# The mz file is now explicitly defined, no need to derive its path
mz_data = np.load(DENOISED_MZ_FILE)
# Use the first spectrum if the data contains multiple spectra
intensity = intensity_data[0] if intensity_data.ndim > 1 else intensity_data

# Create a SpectrumBaseModule instance directly, providing dummy coordinates
sp = SpectrumBaseModule(mz_list=mz_data, intensity=intensity, coordinates=[0, 0, 0])
print(f"Loaded spectrum with {len(sp.mz_list)} data points.")
print("Data loaded successfully.")
# --- Baseline Correction ---
if m is None:
    n = len(sp.intensity)
    m_auto = max(1, min(50, n // 10))
    print(f"SNIP 'm' is None, so it will be automatically determined as: {m_auto}")

print(f"Applying SNIP baseline correction with m={m}, decreasing={decreasing}, epsilon={epsilon}, baseline_scale={baseline_scale}...")
# The method returns a new spectrum object with the corrected intensity.
corrected_sp, estimated_baseline = MSIPreprocessor.baseline_correction(
    data=sp,
    method="snip",
    m=m,
    decreasing=decreasing,
    epsilon=epsilon,
    baseline_scale=baseline_scale
)
# The estimated baseline is now returned directly by the function.

# --- Metrics ---
metrics = compute_metrics(sp.mz_list, sp.intensity, corrected_sp.intensity, estimated_baseline)
print("\n=== Baseline Correction Metrics (SNIP) ===")
print(f"  Negative Ratio: {metrics['negative_ratio']:.4f}")
print(f"  Correlation with Original: {metrics['corr_all']:.4f}")
if metrics['baseline_smoothness'] is not None:
    print(f"  Baseline Smoothness (RMS of 2nd derivative): {metrics['baseline_smoothness']:.6f}")
print(f"  TIC Retention Ratio: {metrics['tic_ratio']:.4f}")
print("========================================\n")

# --- Plotting ---
plt.figure(figsize=(16, 5))
plt.plot(sp.mz_list, sp.intensity, label='Original', color='steelblue', linewidth=1.0, alpha=0.9); plt.plot(sp.mz_list, estimated_baseline, label='Estimated Baseline', color='forestgreen', linewidth=1.0, alpha=0.9); plt.plot(corrected_sp.mz_list, corrected_sp.intensity, label='Corrected', color='darkorange', linewidth=1.0, alpha=0.9)
plt.title('SNIP Baseline Correction Visualization'); plt.xlabel("m/z"); plt.ylabel("Intensity"); plt.xlim(400, 450); plt.ylim(0, 2); plt.legend(); plt.grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig(OUTPUT_FILENAME, dpi=300); plt.close()
# --- Save Corrected Data ---
input_file_stem = Path(FILE_PATH).stem
out_base = f"{input_file_stem}_snip_corrected"
np.save(f"{out_base}_intensity.npy", np.asarray(corrected_sp.intensity, dtype=float))
np.save(f"{out_base}_mz.npy", np.asarray(corrected_sp.mz_list, dtype=float))
print(f"Corrected data saved to {out_base}_intensity.npy and {out_base}_mz.npy")
```

![SNIP 前后示例（待补图）]()

## 参数说明与调参建议

- 通用
  - `baseline_scale`（基线缩放系数）：`0.0–1.0`；减小可避免过度扣底；忠实算法行为建议 `1.0`
- ASLS
  - `lam`（平滑权重）：控制基线曲率，越大越平滑；建议 `1e4–1e8`
  - `p`（非对称权重）：越小越保峰（降低高值点权重）；建议 `0.001–0.1`
  - `niter`（迭代次数）：`5–30`；增加迭代可稳定权重
- SNIP
  - `m`（最大半窗）：默认自动 `min(50, n//10)`；过大可能过度扣底
  - `decreasing`（窗序）：`True`（粗到细）通常更稳健；`False`（细到粗）适用于需要先处理局部细节
  - `epsilon`（早停阈值）：`1e-5–1e-3`；越大越早停，适合防止细节过度削减

## 使用场景

- 背景平滑扣除且保峰形：ASLS（适合定量与形状保持）
- 快速稳健的扣底：SNIP（适合含尖峰与宽峰混合的背景压制）
- 避免过度校正：降低 `baseline_scale` 或提高 `epsilon`（SNIP）

## 常见问题与排查

- 输入类型错误
  - 报错：`TypeError: data must be np.ndarray or SpectrumBaseModule`
  - 解决：确保传入 `np.ndarray` 或 `SpectrumBaseModule`
- 空数组或维度不匹配
  - 情况：长度为 0 或 `mz_list` 与 `intensity` 长度不一致
  - 解决：清洗数据并保证一维长度一致
- 稀疏依赖缺失（ASLS）
  - 情况：环境缺少 `scipy.sparse`；会自动退化为致密实现
  - 影响：计算速度较慢；建议安装 `scipy`
- 过度扣底
  - 症状：校正后峰形过低，TIC 丢失明显
  - 解决：降低 `lam`（ASLS）、减小 `m` 或增大 `epsilon`（SNIP）、提高 `baseline_scale`

## 实践建议

- 基线缩放优先级：优先设置 `baseline_scale`（如 `0.6–0.9`）以避免过度扣底
- 指标评估：可参考示例中的四个指标：
  - 负值比例（`negative_ratio`）
  - 与原始信号的相关性（`corr_all`）
  - 基线光滑度（`baseline_smoothness`，二阶差分 RMS）
  - TIC 保留率（`tic_ratio`）
- 路径管理：Windows 环境使用 `\\` 分隔；NPY 读取时确保强度与 m/z 对齐

## 参考

- `preprocess/baseline_correction.py`（基线校正核心实现：ASLS、SNIP）
- `module/ms_module.py`（`SpectrumBaseModule` 数据结构与可视化）