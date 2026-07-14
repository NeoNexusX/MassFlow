# MassFlow

本文介绍 MassFlow 中的同位素峰过滤（isotope filtering）模块，重点包括：

- 流水线级入口 `Preprocessor.filter_peaks_by_isotopes`
- 扁平批处理入口 `FlatPreprocess.isotope_filter_flat`
- `preprocess/helper/isotope_filter_helper.py` 中的逐谱过滤函数

## 概述

同位素峰过滤用于在峰提取之后识别同一同位素簇中的后续同位素峰，并只保留最低 m/z 的起始峰以及未被分配到同位素簇的峰。

- 输入必须是 **processed + centroided** 数据：每条谱具有自己的峰 m/z 数组。
- 该操作按谱独立执行，不会跨像素合并峰，也不会创建共享 m/z 轴。
- 输出会压缩 m/z、强度和 `lengths`，被识别为后续同位素峰的数据点会被真正删除。
- 推荐顺序为 `peak_pick → filter_peaks_by_isotopes → peak_align`。

### 流程图

```mermaid
graph LR
  A["Preprocessor.filter_peaks_by_isotopes"] --> B["FlatPreprocess.isotope_filter_flat"]
  B --> C["filter_peaks_by_isotopes"]
  C --> D["逐谱构建 keep mask"]
  D --> E["压缩 m/z、intensity 和 lengths"]
  E --> F["保留单同位素峰和未分配峰"]
```

## 核心 API

### Preprocessor.filter_peaks_by_isotopes（流水线级）

```python
massflow.preprocess.api.PreprocessorAPI.filter_peaks_by_isotopes(
    *,
    mz_tolerance: float = 10.0,
    units: Literal["da", "ppm"] = "ppm",
    min_charge: int = 1,
    max_charge: int = 3,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
) -> Self
```

- 描述：在每条 centroided 谱中识别同位素簇，删除簇内除最低 m/z 起始峰以外的峰。
- 任务类型：batch 任务；使用 `Preprocessor.batch_size` 分批读取和写出。
- 输入：应先经过 `peak_pick`，或本身已经是 processed centroid 数据。
- 输出：调用 `.start()` 后返回新的数据管理器；谱数量和坐标保持不变，每条谱的峰数可能减少。
- 参数：
  - `mz_tolerance`：同位素峰匹配容差，默认 `10.0`。
  - `units`：容差单位，`"ppm"` 或 `"da"`，默认 `"ppm"`。
  - `min_charge`、`max_charge`：搜索的电荷范围。搜索顺序从 `max_charge` 递减到 `min_charge`。
  - `min_isotope_count`：一个候选簇被接受所需的最少峰数，包含起始峰。
  - `max_isotope_count`：一个候选簇最多检查的峰数，包含起始峰。
  - `use_decreasing_model`：是否启用后续同位素峰强度递减约束。
  - `start_intensity_check`：从第几个同位素索引开始检查强度递减，默认 `2`。

示例：

```python
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import Preprocessor

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()

filtered = (
    Preprocessor(
        dm,
        batch_size=512,
        output_path="data/isotope_filtered.imzML",
    )
    .peak_pick(method="mad", snr=2.0)
    .filter_peaks_by_isotopes(
        mz_tolerance=10.0,
        units="ppm",
        min_charge=1,
        max_charge=3,
        min_isotope_count=3,
    )
    .start()
)

filtered.close()
dm.close()
```

### filter_peaks_by_isotopes（辅助函数级）

```python
preprocess.helper.isotope_filter_helper.filter_peaks_by_isotopes(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    *,
    mz_tolerance: float = 10.0,
    units: Literal["da", "ppm"] = "ppm",
    min_charge: int = 1,
    max_charge: int = 1,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]
```

- 返回值：`(filtered_mz, filtered_intensity, filtered_lengths)`。
- 说明：辅助函数的 `max_charge` 默认值为 `1`；流水线 API 的默认值为 `3`，直接调用 helper 时需要注意这一差异。
- 输入 m/z 应按每条谱升序排列，因为内部使用 `np.searchsorted` 搜索容差窗口。

### FlatPreprocess.isotope_filter_flat（扁平批处理级）

```python
massflow.preprocess.flat_pre_fun.FlatPreprocess.isotope_filter_flat(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    *,
    mz_tolerance: float = 10.0,
    units: Literal["da", "ppm"] = "ppm",
    min_charge: int = 1,
    max_charge: int = 3,
    min_isotope_count: int = 3,
    max_isotope_count: int = 10,
    use_decreasing_model: bool = True,
    start_intensity_check: int = 2,
) -> FlatBatchResult
```

- 描述：流水线中每个 flat batch 的入口。
- 输出：返回压缩后的 m/z、强度和谱长度；坐标由 pipeline 原样传递。

## 算法说明

当前实现使用固定的碳同位素质量差：

```text
ISOTOPE_MASS_DIFFERENCE = 1.003354838 Da
```

对于起始峰 `M` 和电荷 `z`，第 `k` 个候选同位素峰位置为：

```text
expected_mz = M + k × 1.003354838 / z
```

逐谱处理过程如下：

1. 从最低 m/z 的未分配峰开始。
2. 按 `max_charge → min_charge` 的顺序测试电荷状态。
3. 对每个期望 m/z，在容差窗口内选择误差最小的峰。
4. 遇到缺失峰、已分配峰或不满足强度模型时停止扩展。
5. 若簇大小达到 `min_isotope_count`，保留起始峰并删除其余同位素峰。
6. 未组成有效同位素簇的峰保持不变。

### 强度递减模型

当 `use_decreasing_model=True` 时，只有 `isotope_index >= start_intensity_check` 的候选峰才需要满足“不高于前一个同位素峰”。默认 `start_intensity_check=2`，因此：

- `M+1` 可以高于 `M`；
- 从 `M+2` 开始要求强度不再上升。

## 参数限制

- `mz_tolerance` 必须为有限正数。
- `ppm` 容差不能超过 `100 ppm`。
- `da` 容差不能超过 `0.1 Da`。
- 电荷范围必须满足 `1 <= min_charge <= max_charge`。
- 同位素峰数必须满足 `2 <= min_isotope_count <= max_isotope_count`。
- `start_intensity_check` 必须为非负整数值。

## 使用注意事项

- 不要在 `peak_align` 之后调用本算法。共享 m/z 轴会被识别为 continuous 数据并触发 `ValueError`。
- 该实现只依据同位素质量间距和可选强度趋势判断，不使用分子式、理论同位素丰度或空间相关性。
- 较大的 `max_charge` 和 `max_isotope_count` 会增加每个起始峰的搜索工作量。
- 容差过宽可能把邻近非同位素峰并入簇；容差过窄可能漏掉真实同位素峰。

## 参考

- `preprocess/helper/isotope_filter_helper.py`：输入校验与数组压缩入口
- `preprocess/numba/isotope_filter_numba.py`：逐谱 keep mask 和 Numba 并行实现
- `preprocess/flat_pre_fun.py`：`FlatPreprocess.isotope_filter_flat`
- `preprocess/api.py`：`PreprocessorAPI.filter_peaks_by_isotopes`
- `preprocess/pipeline/planner.py`：峰提取、同位素过滤、参考轴构建与峰对齐的执行顺序

