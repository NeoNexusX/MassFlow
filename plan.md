# Peak Reference Filtering Plan

## 背景

当前 flat pipeline 中，peak alignment 会将每条谱投影到同一个 `reference` m/z 轴上。问题是 `reference_computer` 生成的 reference 可能非常长，后续 `peak_align_flat` 会生成形状为：

```text
n_spectra x reference.size
```

的 dense aligned matrix。如果 reference 中包含大量低支持度、低强度或冗余峰，那么对齐阶段会产生不必要的内存、I/O 和计算开销。

因此 peak filtering 的主要目标应当是：

1. 在正式 align 之前尽量缩短 reference。
2. 保持 `peak_align` 的职责单一，只做投影和可选统计。
3. 将滤峰策略放在独立的 `peak_filter_helper.py` 中，避免把业务规则耦合进 align 内核。

## 总体设计

推荐拆成两层过滤/统计：

```text
peak_pick
  -> reference_computer 生成候选 reference，并基于 reference counts 做预过滤
  -> peak_align 使用过滤后的 reference 正式对齐
  -> peak_align 可选产出 alignment stats
  -> peak_filter_helper 基于 stats 独立实现更丰富的滤峰策略
```

核心原则：

- `reference_computer` 可以使用 counts 做预过滤，因为这是减少 reference 长度、降低 align 成本的最佳位置。
- `Preprocessor.peak_align` 可以接收 `min_count` / `min_frequency` 作为 reference 预过滤参数，但最终 `peak_align_flat` / align numba kernel 不接收这些滤峰规则，避免对齐职责膨胀。
- `peak_align` 可以维护和 reference 等长的统计数组，例如 `counts`、`intensity_sum`、`intensity_max`，供后续分析或二次过滤使用。
- `peak_filter_helper.py` 独立消费 `reference + stats`，返回 `filtered_reference + keep_mask`。

## 模块职责

### `peak_align_helper.py`

继续负责：

- domain/resolution 估计。
- reference 构建。
- raw peaks 到 reference 的 binning。
- 调用 numba align kernel 执行对齐。

需要扩展：

- `merge_peaks` 同步合并 counts。
- `bin_peaks_parallel` 返回 reference 及对应 counts。
- `reference_computer` 增加可选预过滤参数。

不建议加入：

- 同位素过滤逻辑。
- 强度过滤策略。
- 复杂业务规则。

### `peak_filter_helper.py`

负责独立滤峰策略：

- 按出现次数过滤。
- 按出现频率过滤。
- 按平均强度过滤。
- 按最大强度过滤。
- 未来可扩展同位素峰过滤。

它只消费统计结果，不主动执行 alignment，也不关心数据来自 reference 阶段还是 align 阶段。

### `peak_align_numba.py`

继续提供基础匹配/对齐 primitive：

- `search_nearest_jit`
- `search_nearest_flat_parallel_jit`
- `align_spectra_flat_parallel_jit`
- 现有 accumulate/binning 相关 kernel

可扩展：

- 支持 `merge_peaks` 返回 merged counts 的 numba kernel。
- 支持 align 阶段返回 batch-level stats 的 kernel。

不建议：

- 在 numba align kernel 中直接执行过滤策略。

### 可选：`peak_filter_numba.py`

如果后续强度统计、同位素过滤或大规模 mask 计算较重，可以新增：

```text
src/massflow/preprocess/numba/peak_filter_numba.py
```

用于放置纯计算 kernel，例如：

- reference stats reduce。
- isotope candidate search。
- keep mask 计算。

第一阶段可以先不新增，除非 `peak_filter_helper.py` 中的 Python/Numpy 实现成为瓶颈。

## 第一层：Reference 阶段 Counts 预过滤

### 当前问题

当前 `reference_computer` 的流程大致是：

```text
estimate_domain_parallel
  -> bin_peaks_parallel
  -> merge_peaks
  -> reference
```

`bin_peaks_parallel` 内部已经维护了：

```python
counts = np.zeros(domain.size, dtype=np.int64)
```

这个 counts 表示候选 domain bin 被多少条谱/峰支持。但当前 `merge_peaks` 只返回 merged m/z reference，没有同步返回 merged counts，导致最终 reference 和 counts 不再一一对应。

### 修改目标

将 merge 后 reference 和 counts 绑定起来：

```python
merged_reference, merged_counts = merge_peaks(...)
```

合并规则：

- merged m/z 使用 counts 加权平均，保持当前 reference 合并逻辑。
- merged counts 使用 sum 合并。
- NaN 或无效 peaks 不参与合并。
- 输出的 `merged_reference[i]` 必须与 `merged_counts[i]` 一一对应。

### API 建议

`merge_peaks`：

```python
def merge_peaks(
    mz_list: NDArray,
    counts: NDArray,
    tolerance: float,
    tol_method: str = "abs",
) -> tuple[NDArray, NDArray]:
    ...
```

`bin_peaks_parallel`：

```python
def bin_peaks_parallel(...) -> tuple[NDArray, NDArray]:
    ...
    return reference, reference_counts
```

`reference_computer`：

```python
def reference_computer(
    data_manager: MSDataManager,
    reference: Optional[NDArray] = None,
    binfun: str = "median",
    binratio: float = 2.0,
    tolerance: Optional[float] = None,
    units: str = "ppm",
    batch_size: int = 256,
    max_threads: int = 2,
    min_count: Optional[int] = None,
    min_frequency: Optional[float] = None,
    return_counts: bool = False,
) -> tuple[NDArray, float] | tuple[NDArray, float, NDArray]:
    ...
```

### 过滤规则

`min_count`：

```python
keep = reference_counts >= min_count
```

`min_frequency`：

```python
keep = np.minimum(reference_counts, n_spectra) / n_spectra >= min_frequency
```

其中 `n_spectra` 可由 `len(data_manager.ms)` 得到，或在 flat iteration 中累计。
因为 merge 后 counts 采用求和合并，理论上可能超过 `n_spectra`，用于 frequency 时需要封顶到 1.0。

如果同时设置：

```python
keep = count_mask & frequency_mask
```

### 语义说明

Reference 阶段 counts 表示：

> 有多少条谱或候选峰支持这个 reference m/z。

它适合用于 align 前的粗过滤，因为它计算早、成本低、能显著减少 reference 长度。

但它不是最终 alignment hit count。原因：

- reference 构建阶段使用 domain bin。
- merge 后 counts 需要重新聚合。
- align 阶段有自己的 nearest/tolerance 匹配规则。
- align 阶段可能因 reference 间距或 tolerance 产生不同命中结果。

因此 reference counts 主要用于预过滤，不建议作为最终质量评估的唯一依据。

## 第二层：Align 阶段 Stats

### 目标

正式 align 时，可选维护和 reference 等长的统计数组：

```python
counts = np.zeros(reference.size, dtype=np.int64)
intensity_sum = np.zeros(reference.size, dtype=np.float64)
intensity_max = np.zeros(reference.size, dtype=np.float64)
```

这些统计应当严格按最终 align 规则计算，因此它们更适合用于：

- 对齐质量评估。
- 过滤策略调参。
- 二次滤峰。
- 可视化或报告。

### Counts 语义

Align 阶段 counts 应定义为：

> 有多少条 spectrum 在该 reference m/z 上成功匹配到至少一个 raw peak。

不要简单使用：

```python
aligned[:, i] > 0
```

因为真实 intensity 可能为 0，或者经过处理后出现 0。counts 应在“匹配成功”时更新，而不是根据输出强度值推断。

### 多线程注意事项

`align_spectra_flat_parallel_jit` 使用 `prange`。不能在每个 spectrum 并行循环中直接写共享全局：

```python
counts[i] += 1
intensity_sum[i] += value
intensity_max[i] = max(...)
```

这会有 race condition。

更稳的方案：

1. 每个 batch 的 align kernel 返回 batch-level stats。
2. Python 层在 batch 之间聚合全局 stats。
3. 或者 numba 返回 per-spectrum/per-batch stats，再做 reduce。

推荐第一阶段实现 batch-level stats：

```python
aligned, batch_counts, batch_intensity_sum, batch_intensity_max = align_spectra_flat_parallel_with_stats_jit(...)
```

batch 内部避免 race 的方式：

- 简单方案：不在 `prange` 内更新全局 stats，而是在生成 `aligned` 后由另一个 serial/reduction kernel 根据 match mask 统计。
- 更优方案：每个 spectrum 写一行 boolean/intensity stats，然后 reduce；但会增加内存。
- 折中方案：按 spectrum 输出局部 matched indices，再聚合；实现复杂度较高。

第一版可以先做简单可靠的方案，后续再优化。

### 是否默认启用

不建议默认启用 align stats。因为它会增加计算和内存开销。

建议增加参数：

```python
collect_stats: bool = False
```

或者单独方法：

```python
peak_aligner_with_stats(...)
```

避免影响当前 pipeline 默认结果。

## `peak_filter_helper.py` 设计

### 数据结构

建议定义：

```python
from dataclasses import dataclass
import numpy as np

@dataclass(slots=True)
class PeakReferenceStats:
    reference: np.ndarray
    counts: np.ndarray | None = None
    intensity_sum: np.ndarray | None = None
    intensity_max: np.ndarray | None = None
    n_spectra: int | None = None

    @property
    def frequency(self) -> np.ndarray | None:
        if self.counts is None or self.n_spectra is None or self.n_spectra <= 0:
            return None
        return self.counts / float(self.n_spectra)

    @property
    def mean_intensity(self) -> np.ndarray | None:
        if self.counts is None or self.intensity_sum is None:
            return None
        out = np.zeros_like(self.intensity_sum, dtype=np.float64)
        mask = self.counts > 0
        out[mask] = self.intensity_sum[mask] / self.counts[mask]
        return out
```

### 过滤函数

```python
def filter_reference_peaks(
    stats: PeakReferenceStats,
    *,
    min_count: int | None = None,
    min_frequency: float | None = None,
    min_mean_intensity: float | None = None,
    min_max_intensity: float | None = None,
    top_k: int | None = None,
    remove_isotopes: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    ...
```

返回：

```python
filtered_reference, keep_mask
```

其中：

- `filtered_reference = stats.reference[keep_mask]`
- `keep_mask.shape == stats.reference.shape`

### 基础过滤规则

`min_count`：

```python
keep &= stats.counts >= min_count
```

`min_frequency`：

```python
keep &= stats.frequency >= min_frequency
```

`min_mean_intensity`：

```python
keep &= stats.mean_intensity >= min_mean_intensity
```

`min_max_intensity`：

```python
keep &= stats.intensity_max >= min_max_intensity
```

`top_k`：

- 可以按 `mean_intensity`、`max_intensity` 或 `counts` 排序。
- 第一版建议不要默认启用，避免隐式改变 m/z coverage。

### 同位素峰过滤

同位素过滤建议作为第二阶段功能，先不要和 counts/强度过滤混在一起。

后续可设计：

```python
def filter_isotope_peaks(
    reference: np.ndarray,
    intensity_score: np.ndarray,
    *,
    charge_states: tuple[int, ...] = (1, 2),
    mz_tolerance: float = 0.01,
    intensity_ratio_range: tuple[float, float] | None = None,
) -> np.ndarray:
    ...
```

初步策略：

- 对每个 m/z，查找 `+1.0033548378 / z` 附近的峰。
- 如果后者强度/支持度明显更低，可标记为 isotope candidate。
- 保留 monoisotopic peak，删除或标记 isotope peak。

注意：

- 同位素过滤容易误删，第一版建议只提供可选功能。
- 默认关闭。
- 需要单独测试。

## Pipeline 集成建议

### Reference 预过滤集成

在 `PreprocessorAPI.peak_align` 中可增加参数：

```python
min_count: int | None = None
min_frequency: float | None = None
```

传入 `_prepare_peak_align_task`，最终交给 `reference_computer`。

调用链：

```text
Preprocessor.peak_align(min_frequency=...)
  -> _prepare_peak_align_task
  -> reference_computer(min_frequency=...)
  -> filtered reference
  -> FlatPreprocess.peak_align_flat(reference=filtered_reference)
```

这样正式 align 使用的 reference 已经缩短。

### Align stats 集成

建议不要第一版就改 `FlatBatchResult`，否则影响范围较大。

可以后续增加：

```python
@dataclass(slots=True)
class PeakAlignStats:
    counts: np.ndarray
    intensity_sum: np.ndarray
    intensity_max: np.ndarray
    n_spectra: int
```

并在更高层保存或返回。

如果需要在 pipeline 中跨 batch 聚合，可以由 `Preprocessor` 维护：

```python
self._peak_align_stats_accumulator
```

但这会让 pipeline 状态变复杂，不建议第一阶段实现。

第一阶段优先完成 reference 预过滤，解决 reference 太长的主要问题。

## 实施顺序

### Phase 1：Reference counts 预过滤

目标：在 align 前缩短 reference。

任务：

1. 修改 `merge_peaks`，返回 `(merged_reference, merged_counts)`。
2. 修改对应 numba kernel，使 counts 随 m/z merge 同步聚合。
3. 修改 `bin_peaks_parallel`，返回 `(reference, reference_counts)`。
4. 修改 `reference_computer`，支持：
   - `min_count`
   - `min_frequency`
   - `return_counts`
5. 修改 `PreprocessorAPI.peak_align` 和 `_prepare_peak_align_task`，把过滤参数传入 reference 计算。
6. 增加测试：
   - merge 后 counts 与 reference 长度一致。
   - 相邻 peaks merge 时 counts 正确求和。
   - `min_count` 能删除低支持 reference。
   - `min_frequency` 能按 `n_spectra` 比例删除低支持 reference。
   - 不设置过滤参数时，reference 输出与当前行为一致。

验收标准：

- 默认参数下 pipeline 结果不变。
- 设置 `min_count` / `min_frequency` 时 reference 变短。
- 对齐输出的 `lengths` 变为过滤后的 reference 长度。

### Phase 2：Peak filter helper

目标：建立独立滤峰策略模块。

任务：

1. 定义 `PeakReferenceStats`。
2. 实现 `filter_reference_peaks`。
3. 支持 counts/frequency/mean/max intensity 基础规则。
4. 增加纯 numpy 单元测试，不依赖 imzML 数据文件。

验收标准：

- 输入 stats，输出 `filtered_reference` 和 `keep_mask`。
- 多个过滤条件组合时按 AND 语义工作。
- 缺少必要 stats 时给出明确错误。

### Phase 3：Align stats

目标：按最终 align 规则产出统计信息。

任务：

1. 设计 `PeakAlignStats`。
2. 增加可选 stats collection API。
3. 在 batch 粒度返回：
   - counts
   - intensity_sum
   - intensity_max
4. 避免 `prange` race condition。
5. 增加合成 flat 数据测试，验证 stats 与 align 匹配结果一致。

验收标准：

- `collect_stats=False` 时当前 align 输出不变。
- `collect_stats=True` 时 stats 长度等于 reference 长度。
- counts 表示“每个 reference 被多少条 spectrum 命中”，不是 `aligned > 0` 的简单计数。

### Phase 4：同位素过滤

目标：可选移除 isotope candidate。

任务：

1. 设计 isotope filter 参数。
2. 基于 m/z 差值和强度/支持度 score 标记 isotope。
3. 默认关闭。
4. 增加针对合成同位素峰的测试。

验收标准：

- 默认不改变 reference。
- 开启后只过滤明确的 isotope candidate。
- 不影响 counts/frequency 基础过滤。

## 风险点

1. **Counts 语义混淆**

   Reference counts 和 alignment counts 不完全相同。文档和变量名要明确：

   - `reference_counts`: reference 构建阶段的支持度。
   - `alignment_counts`: 最终 align 匹配命中数。

2. **Merge 后 counts 错位**

   修改 `merge_peaks` 时必须保证 m/z 和 counts 同步 merge，同步过滤。

3. **多线程 stats race**

   align 阶段如果在 `prange` 中更新共享 counts/sum/max，会产生不稳定结果。必须避免。

4. **默认行为改变**

   新增过滤参数必须默认关闭。默认情况下 reference 和 alignment 输出应保持当前行为。

5. **同位素误删**

   同位素过滤不应第一阶段默认启用，需要单独验证。

## 推荐默认参数

第一版建议：

```python
min_count = None
min_frequency = None
return_counts = False
collect_stats = False
remove_isotopes = False
```

这样默认 pipeline 结果不变，用户显式开启过滤时才缩短 reference。

## 最终目标

理想状态下，用户可以这样使用：

```python
processed = (
    Preprocessor(dm)
    .peak_pick(method="quantile", snr=2.0)
    .peak_align(
        units="ppm",
        binfun="median",
        binratio=2.0,
        min_frequency=0.01,
    )
    .start()
)
```

其中：

- `reference_computer` 根据 `min_frequency` 先缩短 reference。
- `peak_align_flat` 只对过滤后的 reference 做正式对齐。
- 未来如果需要更复杂规则，可以通过 `peak_filter_helper.py` 独立组合 stats 和过滤逻辑。
