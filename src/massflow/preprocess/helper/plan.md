# Adduct Filtering Implementation Plan

> **For agentic workers:** 先按本计划实现。不要顺手重构无关模块；不要跑完整 pytest。默认只做静态检查，除非用户另行明确要求运行测试。

**Goal:** 参考现成 LC/MS adduct 注释算法，实现 MassFlow 的加合离子过滤：从 aligned/shared m/z feature 轴中识别同一 neutral mass 的不同 ion/adduct forms，保留代表峰，过滤冗余 adduct peaks。

**Architecture:** 不采用简单 `mz + delta` 的 per-spectrum 删除。第一版采用 CAMERA `findAdducts` / MZmine Ion Identity Networking 共同使用的 neutral-mass hypothesis 思路：对每个 m/z feature 和每条 adduct rule 反推 neutral mass，把相近 neutral mass hypothesis 聚类，再用 MSI 的空间强度 profile correlation 做确认。过滤本身沿用 `stat_filter_helper.py` 的 dataset/shared-mask 写回路径。

**Tech Stack:** Python 3.12, NumPy, Numba where useful, existing `MSDataManagerImzML.flat_generator`, existing dataset-scope preprocessing pipeline.

---

## 调研结论

### CAMERA / Bioconductor 是最适合第一版参考的算法

官方包说明：CAMERA 是 mass spectrometry annotation 工具，包含 rule-based isotope/adduct annotation 和 EIC correlation tagging。官方 `findAdducts` 文档说明它对 `xsAnnotate` 的 pseudospectra 注释 adducts/fragments，参数包括 `ppm`、`mzabs`、`multiplier`、`polarity`、`rules`、`max_peaks`；细节中明确：对每个由 `groupFWHM` / `groupCorr` 生成的 pseudospectra group 计算所有可能 adduct 并映射到 peaks，至少两个 adduct 匹配时才能计算 possible molecule mass，且推荐先跑 isotope annotation。

已核对 CAMERA 1.68.0 官方源码包：

- `R/fct_findAdducts.R`
  - `massDiffMatrix(m, rules)` 使用公式：
    `neutral_mass = (abs(charge) * mz - massdiff) / nmol`
  - `annotateGrp(...)` 对一个 pseudospectrum group 内的 peaks 生成 neutral-mass hypotheses。
  - `createHypothese(...)` 将相近 neutral masses 聚类，cluster 中至少有 2 个 candidates 才进入后续判断。
  - 后续有 isotope charge consistency、quasimolecular ion、OID causality、IPS score 等过滤。
- `inst/rules/primary_adducts_pos.csv`
  - positive primary rules: `[M+H]+`, `[M+Na]+`, `[M+K]+`, `[M+NH4]+`
- `inst/rules/primary_adducts_neg.csv`
  - negative primary rules: `[M-H]-`, `[M-2H+Na]-`, `[M-2H+K]-`, `[M+Cl]-`

结论：第一版应参考 CAMERA 的 rule table + neutral mass clustering，而不是写一个没有来源的固定 offset 删除。

### MZmine 提供了更现代的网络化实现，但第一版不必完整照搬

已核对 MZmine `master` 官方源码：

- `IonNetworkingTask.java`
  - 要求 feature list 已经有 correlation groups；没有 groups 会提示先跑 correlate grouping。
  - 在每个 group 内只比较 correlated rows。
  - 对 row pair 调用 `library.searchRows(rowA, rowB, mzTolerance)`，匹配后合并为 ion identity networks。
- `SearchableIonLibrary.java`
  - 按 charge 和 molecule count 对 ion types 分组并按质量排序，做 early-exit 搜索。
  - row pair 的匹配本质上也是把两个 observed m/z 通过不同 ion types 反推到一致 neutral mass。
- `IonType.java`
  - neutral mass 公式同 CAMERA 等价：
    `neutral_mass = (mz * abs(charge) - totalMass) / molecules`
  - 从 neutral mass 正推 m/z：
    `mz = (neutral_mass * molecules + totalMass) / abs(charge)`
- `id_cliquems/.../ComputeAdduct.java` 和 `AdductAnnotationCliqueMS.java`
  - CliqueMS 路线会先用 isotope/clique groups，再对 candidate neutral masses 打分并输出 top annotations。
  - 这个路线更完整，但对当前 MassFlow 第一版过重。

结论：MZmine 说明 adduct annotation 不应只靠 m/z 差值；它至少需要 feature group/correlation。MassFlow 的 MSI 版本可以用“ion image spatial correlation”替代 LC-MS 的 EIC/RT correlation。

### 对 MassFlow 的实现选择

MassFlow 当前没有 LC-MS RT/EIC pseudospectrum 概念；如果在单个 pixel spectrum 内直接删 `mz + delta`，误删风险高，且不符合 CAMERA/MZmine 的算法边界。第一版应做 dataset-scope shared m/z 过滤：

1. 要求输入为 aligned/shared m/z 数据。
2. 把每个 m/z feature 当作一个 ion feature。
3. 用 CAMERA primary rules 为每个 feature 反推 neutral mass hypotheses。
4. 对 neutral mass 相近、且来自至少两个不同 peaks/rules 的 hypothesis cluster 建 adduct group。
5. 用各 peak across-pixels intensity profile 的相关性确认 group。
6. 每个 group 只保留代表峰，过滤其余 adduct peaks。

这和现有 `filter_peaks_by_stats` 的 shared peak mask 写回路径一致。

## 正确性与溯源门槛

这个计划不能承诺“在所有 MSI 数据上算法一定有效”；它能承诺的是：实现不得脱离可溯源算法，且必须产出足够证据供审核判断是否正确有效。GLM5.2 实现前必须满足以下门槛：

1. **Source parity**
   - 默认 adduct rule 常数必须逐条对应 CAMERA `primary_adducts_pos.csv` / `primary_adducts_neg.csv`，在代码注释或测试中保留 rule 名称、`nmol`、`charge`、`massdiff`、`quasi`、`ips` 来源。
   - neutral mass 公式必须同时等价于 CAMERA `massDiffMatrix(...)` 和 MZmine `IonType.getMass(...)`：
     `neutral_mass = (abs(charge) * mz - massdiff) / nmol`。
   - candidate group 必须保留 CAMERA/MZmine 的核心条件：同一 neutral mass hypothesis、至少两个不同 peaks/rules；不得实现成简单固定 m/z 差值删除。
2. **Adaptation boundary**
   - “spatial profile correlation 替代 LC-MS EIC/RT correlation”是 MassFlow/MSI 适配层，不是 CAMERA/MZmine 原始算法本身。
   - 如果实现者不能解释该适配如何用 MSI ion-image across-pixel intensity profile 支撑同源 ion feature 判断，就不要把它写成默认强结论；必须保留为可调参数，并在报告中输出 correlation。
3. **No source, no default**
   - 不要新增未溯源的 adduct rules、score、intensity heuristic、经验阈值或删除策略。
   - 若确需加入项目特有参数，必须命名为显式用户参数，并在 docstring 中标注“MassFlow adaptation”，不能写成 CAMERA/MZmine 规则。
4. **Auditability**
   - public API 必须能输出 annotation/report table，至少包含：`neutral_mass`、`representative_peak`、`removed_peak`、`representative_mz`、`removed_mz`、`representative_rule`、`removed_rule`、`spatial_correlation`、`rule_source`。
   - 没有 report table 的实现视为不可审核；即使过滤结果看起来合理，也不能通过审核。
5. **Validation gates**
   - 单元测试必须复现 CAMERA primary rules 下的 neutral mass 公式。
   - synthetic toy case 必须构造已知 neutral mass 的 `[M+H]+`、`[M+Na]+`、`[M+K]+`、`[M+NH4]+`，验证只保留代表峰。
   - negative control 必须覆盖：无 neutral-mass rule match 的相近峰不能删除；有 neutral-mass rule match 但 spatial correlation 低于阈值的峰不能删除。
   - 小型真实或 fixture artifact 检查必须确认输出 shared m/z axis 只删除 report 中列出的 peaks，且 kept/removed index 与 report 一致。

## 需要修改的文件

- `src/massflow/preprocess/helper/adduct_filter_helper.py`
  - 放规则表、neutral mass hypothesis、spatial profile correlation、keep mask 构建、annotation report、public dataset filter。
- `src/massflow/preprocess/numba/adduct_filter_numba.py`
  - 第一版只放必要的 candidate profile aggregate kernel；不要写 per-spectrum `mz + delta` filter。
- `src/massflow/preprocess/api.py`
  - 增加 `filter_peaks_by_adducts(...)` dataset-scope API。
- `src/massflow/preprocess/pipeline/planner.py`
  - 将 adduct filter 放在 `peak_align` 后面；如果用户没有先 peak-align，API 应提示该任务需要 shared m/z axis。
- `tests/test_adduct_filter.py`
  - 添加小型 shared m/z toy tests。

## 核心数据结构

### Task 1: 定义 CAMERA primary rule table

**Files:**
- Modify: `src/massflow/preprocess/helper/adduct_filter_helper.py`

- [ ] **Step 1: 替换当前 Markdown 任务说明为 Python 模块**

```python
"""Helpers for CAMERA-like adduct filtering on shared m/z feature axes."""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.numba.stat_filter_numba import filter_peaks_flat_parallel_jit
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.adduct_filter_helper")

Polarity = Literal["positive", "negative"]
```

- [ ] **Step 2: 增加 rule dataclass**

```python
@dataclass(frozen=True, slots=True)
class AdductRule:
    """CAMERA-style adduct rule."""

    name: str
    nmol: int
    charge: int
    massdiff: float
    score: float
    quasi: bool
    priority: int
```

- [ ] **Step 3: 增加 source-backed primary rules**

数值来自 CAMERA 1.68.0 `inst/rules/primary_adducts_pos.csv` 和 `inst/rules/primary_adducts_neg.csv`。

```python
PRIMARY_ADDUCT_RULES_POSITIVE: tuple[AdductRule, ...] = (
    AdductRule("[M+H]+", 1, 1, 1.007276, 1.0, True, 0),
    AdductRule("[M+Na]+", 1, 1, 22.989218, 1.0, True, 1),
    AdductRule("[M+K]+", 1, 1, 38.963158, 1.0, True, 2),
    AdductRule("[M+NH4]+", 1, 1, 18.033823, 1.0, True, 3),
)

PRIMARY_ADDUCT_RULES_NEGATIVE: tuple[AdductRule, ...] = (
    AdductRule("[M-H]-", 1, -1, -1.007276, 1.0, True, 0),
    AdductRule("[M-2H+Na]-", 1, -1, 20.974666, 0.5, False, 1),
    AdductRule("[M-2H+K]-", 1, -1, 36.948606, 0.5, False, 2),
    AdductRule("[M+Cl]-", 1, -1, 34.969402, 1.0, True, 3),
)
```

- [ ] **Step 4: 增加 rule resolver**

```python
def _resolve_rules(
    polarity: Polarity,
    rules: tuple[AdductRule, ...] | None,
) -> tuple[AdductRule, ...]:
    if rules is not None:
        return rules
    if polarity == "positive":
        return PRIMARY_ADDUCT_RULES_POSITIVE
    if polarity == "negative":
        return PRIMARY_ADDUCT_RULES_NEGATIVE
    raise ValueError("polarity must be either 'positive' or 'negative'.")
```

### Task 2: 实现 neutral mass hypothesis clustering

**Files:**
- Modify: `src/massflow/preprocess/helper/adduct_filter_helper.py`

- [ ] **Step 1: 增加 hypothesis dataclass**

```python
@dataclass(frozen=True, slots=True)
class AdductHypothesis:
    peak_index: int
    rule_index: int
    neutral_mass: float
```

- [ ] **Step 2: 增加 neutral mass 公式**

必须和 CAMERA / MZmine 公式一致。

```python
def _neutral_mass(mz: float, rule: AdductRule) -> float:
    return (abs(rule.charge) * mz - rule.massdiff) / rule.nmol
```

- [ ] **Step 3: 为 shared m/z axis 生成 hypotheses**

```python
def _build_hypotheses(
    mz_axis: np.ndarray,
    rules: tuple[AdductRule, ...],
) -> list[AdductHypothesis]:
    hypotheses: list[AdductHypothesis] = []
    for peak_index, mz in enumerate(mz_axis):
        for rule_index, rule in enumerate(rules):
            hypotheses.append(
                AdductHypothesis(
                    peak_index=peak_index,
                    rule_index=rule_index,
                    neutral_mass=_neutral_mass(float(mz), rule),
                )
            )
    hypotheses.sort(key=lambda item: item.neutral_mass)
    return hypotheses
```

- [ ] **Step 4: 增加 CAMERA-style neutral mass tolerance**

CAMERA `createHypothese` 对 candidate neutral masses 用 `2 * ppm * mean_mass + mzabs` 量级聚类。这里保留同样语义。

```python
def _neutral_mass_tolerance(
    neutral_mass: float,
    *,
    ppm: float,
    mzabs: float,
) -> float:
    return 2.0 * ppm * 1e-6 * neutral_mass + mzabs
```

- [ ] **Step 5: 聚类 hypotheses**

```python
def _cluster_hypotheses(
    hypotheses: list[AdductHypothesis],
    *,
    ppm: float,
    mzabs: float,
) -> list[list[AdductHypothesis]]:
    clusters: list[list[AdductHypothesis]] = []
    current: list[AdductHypothesis] = []

    for hypothesis in hypotheses:
        if not current:
            current = [hypothesis]
            continue

        center = float(np.mean([item.neutral_mass for item in current]))
        tolerance = _neutral_mass_tolerance(center, ppm=ppm, mzabs=mzabs)
        if abs(hypothesis.neutral_mass - center) <= tolerance:
            current.append(hypothesis)
        else:
            clusters.append(current)
            current = [hypothesis]

    if current:
        clusters.append(current)

    return clusters
```

### Task 3: 从 clusters 生成候选 adduct groups

**Files:**
- Modify: `src/massflow/preprocess/helper/adduct_filter_helper.py`

- [ ] **Step 1: 增加 group dataclass**

```python
@dataclass(frozen=True, slots=True)
class AdductGroup:
    neutral_mass: float
    peak_indices: tuple[int, ...]
    rule_indices: tuple[int, ...]
    representative_peak: int
    spatial_correlations: tuple[float, ...]
```

`rule_indices` 和 `spatial_correlations` 必须与 `peak_indices` 同序，便于后续 report 精确说明每个 removed peak 对应的 rule 和 correlation。

- [ ] **Step 2: 增加 report dataclass**

```python
@dataclass(frozen=True, slots=True)
class AdductAnnotation:
    neutral_mass: float
    representative_peak: int
    removed_peak: int
    representative_mz: float
    removed_mz: float
    representative_rule: str
    removed_rule: str
    spatial_correlation: float
    rule_source: str
```

- [ ] **Step 3: 实现代表峰选择**

优先保留 quasi/main ion；tie 用 rule priority；再 tie 用 summed intensity。

```python
def _choose_representative_peak(
    cluster: list[AdductHypothesis],
    rules: tuple[AdductRule, ...],
    summed_intensity: np.ndarray,
) -> int:
    best = min(
        cluster,
        key=lambda item: (
            not rules[item.rule_index].quasi,
            rules[item.rule_index].priority,
            -float(summed_intensity[item.peak_index]),
            item.peak_index,
        ),
    )
    return best.peak_index
```

- [ ] **Step 4: 实现 cluster filtering**

```python
def _choose_rule_per_peak(
    cluster: list[AdductHypothesis],
    rules: tuple[AdductRule, ...],
) -> dict[int, int]:
    selected: dict[int, int] = {}
    for item in cluster:
        current = selected.get(item.peak_index)
        if current is None:
            selected[item.peak_index] = item.rule_index
            continue
        current_rule = rules[current]
        item_rule = rules[item.rule_index]
        if (not item_rule.quasi, item_rule.priority) < (
            not current_rule.quasi,
            current_rule.priority,
        ):
            selected[item.peak_index] = item.rule_index
    return selected


def _candidate_groups_from_clusters(
    clusters: list[list[AdductHypothesis]],
    rules: tuple[AdductRule, ...],
    summed_intensity: np.ndarray,
    *,
    require_quasi: bool,
) -> list[AdductGroup]:
    groups: list[AdductGroup] = []
    for cluster in clusters:
        rule_by_peak = _choose_rule_per_peak(cluster, rules)
        peak_indices = tuple(sorted(rule_by_peak))
        rule_indices = tuple(rule_by_peak[peak_index] for peak_index in peak_indices)
        if len(peak_indices) < 2 or len(set(rule_indices)) < 2:
            continue
        if require_quasi and not any(rules[item.rule_index].quasi for item in cluster):
            continue

        representative = _choose_representative_peak(cluster, rules, summed_intensity)
        groups.append(
            AdductGroup(
                neutral_mass=float(np.mean([item.neutral_mass for item in cluster])),
                peak_indices=peak_indices,
                rule_indices=rule_indices,
                representative_peak=representative,
                spatial_correlations=tuple(np.nan for _ in peak_indices),
            )
        )
    return groups
```

### Task 4: 实现 shared m/z dataset profile 统计

**Files:**
- Modify: `src/massflow/preprocess/helper/adduct_filter_helper.py`

- [ ] **Step 1: 验证输入必须是 shared m/z**

```python
def _validate_shared_mz_data(data_manager: MSDataManagerImzML) -> None:
    if data_manager.ms.meta.continuous is not True:
        raise ValueError("Adduct filtering requires aligned/shared m/z data. Run peak_align first.")
    if data_manager.ms.shared_mz_list is None:
        raise ValueError("Adduct filtering requires a shared m/z axis.")
```

- [ ] **Step 2: 计算 summed intensity**

复用 `stat_filter_helper.compute_peak_filter_profile(...)` 最省事；如果直接 import 会产生循环或职责不清，再复制 `_aggregate_batches(...)` 的最小逻辑。

```python
def _compute_summed_intensity(data_manager: MSDataManagerImzML) -> np.ndarray:
    summed = None
    for mz_axis, intensity_flat, lengths, _ in data_manager.flat_generator(include_mz=True):
        batch_sum = np.zeros(mz_axis.size, dtype=np.float64)
        offsets = np.concatenate(([0], np.cumsum(lengths, dtype=np.int64)))
        for spectrum_index in range(lengths.size):
            start = offsets[spectrum_index]
            end = offsets[spectrum_index + 1]
            batch_sum[: int(lengths[spectrum_index])] += intensity_flat[start:end]
        summed = batch_sum if summed is None else summed + batch_sum
    if summed is None:
        raise ValueError("Cannot compute adduct groups from empty data.")
    return summed
```

- [ ] **Step 3: 计算候选 group 内 peak profile correlation**

第一版只对 candidate groups 内的 peaks 计算相关性，不构建全 feature x feature correlation matrix。

```python
def _collect_candidate_profiles(
    data_manager: MSDataManagerImzML,
    peak_indices: tuple[int, ...],
) -> dict[int, np.ndarray]:
    profiles = {peak_index: [] for peak_index in peak_indices}
    for _mz_axis, intensity_flat, lengths, _ in data_manager.flat_generator(include_mz=False):
        offsets = np.concatenate(([0], np.cumsum(lengths, dtype=np.int64)))
        for spectrum_index in range(lengths.size):
            start = offsets[spectrum_index]
            length = int(lengths[spectrum_index])
            for peak_index in peak_indices:
                value = intensity_flat[start + peak_index] if peak_index < length else 0.0
                profiles[peak_index].append(float(value))
    return {key: np.asarray(values, dtype=np.float64) for key, values in profiles.items()}
```

```python
def _pearson_correlation(a: np.ndarray, b: np.ndarray) -> float:
    a_centered = a - np.mean(a)
    b_centered = b - np.mean(b)
    denom = np.linalg.norm(a_centered) * np.linalg.norm(b_centered)
    if denom == 0.0:
        return 0.0
    return float(np.dot(a_centered, b_centered) / denom)
```

- [ ] **Step 4: 用 correlation 确认 groups**

```python
def _confirm_groups_by_correlation(
    data_manager: MSDataManagerImzML,
    groups: list[AdductGroup],
    *,
    min_spatial_correlation: float,
) -> list[AdductGroup]:
    confirmed: list[AdductGroup] = []
    for group in groups:
        profiles = _collect_candidate_profiles(data_manager, group.peak_indices)
        representative = profiles[group.representative_peak]
        correlated_peaks: list[int] = []
        correlated_rules: list[int] = []
        correlations: list[float] = []
        for peak_index, rule_index in zip(group.peak_indices, group.rule_indices, strict=True):
            if peak_index == group.representative_peak:
                corr = 1.0
            else:
                corr = _pearson_correlation(representative, profiles[peak_index])
            if peak_index == group.representative_peak or corr >= min_spatial_correlation:
                correlated_peaks.append(peak_index)
                correlated_rules.append(rule_index)
                correlations.append(corr)
        if len(correlated_peaks) >= 2 and len(set(correlated_rules)) >= 2:
            confirmed.append(
                AdductGroup(
                    neutral_mass=group.neutral_mass,
                    peak_indices=tuple(correlated_peaks),
                    rule_indices=tuple(correlated_rules),
                    representative_peak=group.representative_peak,
                    spatial_correlations=tuple(correlations),
                )
            )
    return confirmed
```

后续性能优化可以把 `_collect_candidate_profiles` 改成一次 scan 收集所有 candidate peaks，或写 Numba pair-aggregate kernel。第一版先保持正确和可审。

### Task 5: 生成 keep mask 并写出 filtered dataset

**Files:**
- Modify: `src/massflow/preprocess/helper/adduct_filter_helper.py`

- [ ] **Step 1: 从 confirmed groups 生成 keep mask**

```python
def _build_adduct_keep_mask(
    mz_axis_size: int,
    groups: list[AdductGroup],
) -> np.ndarray:
    keep_mask = np.ones(mz_axis_size, dtype=bool)
    for group in groups:
        for peak_index in group.peak_indices:
            if peak_index != group.representative_peak:
                keep_mask[peak_index] = False
    return keep_mask
```

- [ ] **Step 2: 生成 annotation report rows**

```python
def _build_adduct_annotations(
    groups: list[AdductGroup],
    mz_axis: np.ndarray,
    rules: tuple[AdductRule, ...],
) -> list[AdductAnnotation]:
    annotations: list[AdductAnnotation] = []
    for group in groups:
        representative_offset = group.peak_indices.index(group.representative_peak)
        representative_rule = rules[group.rule_indices[representative_offset]]
        for peak_index, rule_index, correlation in zip(
            group.peak_indices,
            group.rule_indices,
            group.spatial_correlations,
            strict=True,
        ):
            if peak_index == group.representative_peak:
                continue
            removed_rule = rules[rule_index]
            annotations.append(
                AdductAnnotation(
                    neutral_mass=group.neutral_mass,
                    representative_peak=group.representative_peak,
                    removed_peak=peak_index,
                    representative_mz=float(mz_axis[group.representative_peak]),
                    removed_mz=float(mz_axis[peak_index]),
                    representative_rule=representative_rule.name,
                    removed_rule=removed_rule.name,
                    spatial_correlation=float(correlation),
                    rule_source="CAMERA primary adduct rules",
                )
            )
    return annotations
```

- [ ] **Step 3: 可选写出 CSV audit report**

```python
def _write_adduct_report(
    annotations: list[AdductAnnotation],
    report_path: str | Path | None,
) -> None:
    if report_path is None:
        return

    import csv

    fieldnames = [
        "neutral_mass",
        "representative_peak",
        "removed_peak",
        "representative_mz",
        "removed_mz",
        "representative_rule",
        "removed_rule",
        "spatial_correlation",
        "rule_source",
    ]
    with Path(report_path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for annotation in annotations:
            writer.writerow(
                {
                    "neutral_mass": annotation.neutral_mass,
                    "representative_peak": annotation.representative_peak,
                    "removed_peak": annotation.removed_peak,
                    "representative_mz": annotation.representative_mz,
                    "removed_mz": annotation.removed_mz,
                    "representative_rule": annotation.representative_rule,
                    "removed_rule": annotation.removed_rule,
                    "spatial_correlation": annotation.spatial_correlation,
                    "rule_source": annotation.rule_source,
                }
            )
```

- [ ] **Step 4: 复制 shared-mask 写回逻辑**

参考 `stat_filter_helper._filter_data_by_mask(...)`，不要改 `stat_filter_helper.py`。

```python
def _filter_data_by_mask(
    data_manager: MSDataManagerImzML,
    keep_mask: np.ndarray,
    *,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
) -> MSDataManagerImzML:
    filtered_data = MSDataManagerImzML(
        filepath=None if output_path is None else str(output_path),
        temp_dir=None if temp_dir is None else str(temp_dir),
        max_threads=data_manager.max_threads,
        mz_dtype=data_manager.mz_dtype,
        intensity_dtype=data_manager.intensity_dtype,
    )
    filtered_data.copy_meta(data_manager)
    filtered_data.ms.meta.continuous = True
    filtered_data.ms.meta.processed = False

    filtered_mz = None
    for mz_axis, intensity_flat, lengths, coordinates in data_manager.flat_generator(
        include_mz=True,
    ):
        if filtered_mz is None:
            filtered_mz = mz_axis[keep_mask]
        filtered_intensity, filtered_lengths = filter_peaks_flat_parallel_jit(
            intensity_flat,
            lengths,
            keep_mask,
        )
        filtered_data.write_flat_batch(
            intensity_flat=filtered_intensity,
            lengths=filtered_lengths,
            coordinates=coordinates,
            mz_axis=filtered_mz,
        )

    filtered_data.close_writer()
    filtered_data.load_head_data()
    return filtered_data
```

- [ ] **Step 5: public API**

```python
def filter_peaks_by_adducts(
    data_manager: MSDataManagerImzML,
    *,
    polarity: Polarity = "positive",
    rules: tuple[AdductRule, ...] | None = None,
    ppm: float = 5.0,
    mzabs: float = 0.015,
    require_quasi: bool = True,
    min_spatial_correlation: float = 0.7,
    output_path: str | Path | None = None,
    temp_dir: str | Path | None = None,
    report_path: str | Path | None = None,
) -> MSDataManagerImzML:
    """Filter adduct peaks using CAMERA-like neutral-mass hypotheses."""
    _validate_shared_mz_data(data_manager)
    if ppm <= 0 or not np.isfinite(ppm):
        raise ValueError("ppm must be a finite positive number.")
    if mzabs < 0 or not np.isfinite(mzabs):
        raise ValueError("mzabs must be a finite non-negative number.")
    if not -1 <= min_spatial_correlation <= 1:
        raise ValueError("min_spatial_correlation must be in [-1, 1].")

    selected_rules = _resolve_rules(polarity, rules)
    mz_axis = np.asarray(data_manager.ms.shared_mz_list)
    summed_intensity = _compute_summed_intensity(data_manager)

    hypotheses = _build_hypotheses(mz_axis, selected_rules)
    clusters = _cluster_hypotheses(hypotheses, ppm=ppm, mzabs=mzabs)
    groups = _candidate_groups_from_clusters(
        clusters,
        selected_rules,
        summed_intensity,
        require_quasi=require_quasi,
    )
    confirmed_groups = _confirm_groups_by_correlation(
        data_manager,
        groups,
        min_spatial_correlation=min_spatial_correlation,
    )
    annotations = _build_adduct_annotations(confirmed_groups, mz_axis, selected_rules)
    _write_adduct_report(annotations, report_path)
    keep_mask = _build_adduct_keep_mask(mz_axis.size, confirmed_groups)
    return _filter_data_by_mask(
        data_manager,
        keep_mask,
        output_path=output_path,
        temp_dir=temp_dir,
    )
```

### Task 6: 接入 pipeline API

**Files:**
- Modify: `src/massflow/preprocess/api.py`
- Modify: `src/massflow/preprocess/pipeline/planner.py`

- [ ] **Step 1: 在 `api.py` import helper**

```python
from massflow.preprocess.helper.adduct_filter_helper import (
    filter_peaks_by_adducts,
    AdductRule,
    Polarity,
)
```

- [ ] **Step 2: 在 `PreprocessorAPI` 增加 dataset task**

放在 `filter_peaks_by_isotopes(...)` 后或 `filter_peaks_by_stats(...)` 附近均可；这是 dataset-scope task。

```python
    def filter_peaks_by_adducts(
        self,
        *,
        polarity: Polarity = "positive",
        rules: tuple[AdductRule, ...] | None = None,
        ppm: float = 5.0,
        mzabs: float = 0.015,
        require_quasi: bool = True,
        min_spatial_correlation: float = 0.7,
        output_path: str | None = None,
        report_path: str | None = None,
    ) -> Self:
        """Register CAMERA-like adduct filtering on aligned/shared m/z data."""
        return self._register_task(
            "filter_peaks_by_adducts",
            scope="dataset",
            apply_fn=filter_peaks_by_adducts,
            polarity=polarity,
            rules=rules,
            ppm=ppm,
            mzabs=mzabs,
            require_quasi=require_quasi,
            min_spatial_correlation=min_spatial_correlation,
            output_path=output_path,
            report_path=report_path,
        )
```

- [ ] **Step 3: 在 planner 顺序中放到 peak_align 后**

```python
_OPERATION_ORDER: dict[str, int] = {
    "baseline_correction": 10,
    "noise_reduction": 20,
    "normalization": 30,
    "peak_pick": 40,
    "filter_peaks_by_isotopes": 50,
    "reference_computer": 60,
    "peak_align": 70,
    "filter_peaks_by_adducts": 80,
}
```

原因：adduct filter 第一版要求 shared m/z feature axis。若用户链式调用 `.peak_pick().filter_peaks_by_adducts()` 而没有 `.peak_align()`，运行时应由 `_validate_shared_mz_data` 报错。

### Task 7: 添加聚焦测试

**Files:**
- Create: `tests/test_adduct_filter.py`

默认不要运行 pytest。测试代码先写好，后续用户允许时只跑该文件。

- [ ] **Step 1: 测 neutral mass formula**

```python
import numpy as np

from massflow.preprocess.helper.adduct_filter_helper import (
    PRIMARY_ADDUCT_RULES_POSITIVE,
    _neutral_mass,
    _build_hypotheses,
    _cluster_hypotheses,
    _candidate_groups_from_clusters,
)


def test_neutral_mass_formula_matches_camera_rules() -> None:
    h_rule = PRIMARY_ADDUCT_RULES_POSITIVE[0]
    na_rule = PRIMARY_ADDUCT_RULES_POSITIVE[1]

    assert np.isclose(_neutral_mass(101.007276, h_rule), 100.0)
    assert np.isclose(_neutral_mass(122.989218, na_rule), 100.0)
```

- [ ] **Step 2: 测 CAMERA-like cluster 能把 H/Na/K/NH4 归为一个候选 group**

```python
def test_candidate_groups_from_primary_positive_rules() -> None:
    mz_axis = np.array(
        [
            101.007276,  # [M+H]+ for neutral mass 100
            118.033823,  # [M+NH4]+ for neutral mass 100
            122.989218,  # [M+Na]+ for neutral mass 100
            138.963158,  # [M+K]+ for neutral mass 100
            250.0,
        ],
        dtype=np.float64,
    )
    summed = np.array([100.0, 20.0, 40.0, 10.0, 500.0], dtype=np.float64)

    hypotheses = _build_hypotheses(mz_axis, PRIMARY_ADDUCT_RULES_POSITIVE)
    clusters = _cluster_hypotheses(hypotheses, ppm=5.0, mzabs=0.015)
    groups = _candidate_groups_from_clusters(
        clusters,
        PRIMARY_ADDUCT_RULES_POSITIVE,
        summed,
        require_quasi=True,
    )

    target_groups = [group for group in groups if set(group.peak_indices) == {0, 1, 2, 3}]
    assert len(target_groups) == 1
    assert target_groups[0].representative_peak == 0
```

- [ ] **Step 3: 测 keep mask 只保留代表峰**

```python
from massflow.preprocess.helper.adduct_filter_helper import AdductGroup, _build_adduct_keep_mask


def test_build_adduct_keep_mask_removes_non_representatives() -> None:
    group = AdductGroup(
        neutral_mass=100.0,
        peak_indices=(0, 1, 2, 3),
        rule_indices=(0, 1, 2, 3),
        representative_peak=0,
        spatial_correlations=(1.0, 0.98, 0.97, 0.99),
    )

    keep = _build_adduct_keep_mask(5, [group])

    np.testing.assert_array_equal(
        keep,
        np.array([True, False, False, False, True]),
    )
```

- [ ] **Step 4: 测 negative controls**

```python
from massflow.preprocess.helper.adduct_filter_helper import _pearson_correlation


def test_low_spatial_correlation_does_not_support_group_confirmation() -> None:
    representative = np.array([10.0, 20.0, 30.0, 40.0])
    unrelated = np.array([40.0, 30.0, 20.0, 10.0])

    assert _pearson_correlation(representative, unrelated) < 0.7
```

```python
def test_unrelated_mz_peak_is_not_removed_by_neutral_mass_group() -> None:
    mz_axis = np.array([101.007276, 122.989218, 250.0], dtype=np.float64)
    summed = np.array([100.0, 40.0, 500.0], dtype=np.float64)

    hypotheses = _build_hypotheses(mz_axis, PRIMARY_ADDUCT_RULES_POSITIVE)
    clusters = _cluster_hypotheses(hypotheses, ppm=5.0, mzabs=0.015)
    groups = _candidate_groups_from_clusters(
        clusters,
        PRIMARY_ADDUCT_RULES_POSITIVE,
        summed,
        require_quasi=True,
    )

    keep = _build_adduct_keep_mask(mz_axis.size, groups)

    assert keep[2]
```

- [ ] **Step 5: 后续如需集成测试**

用一个 fake/shared-mz data manager 或小 imzML fixture 验证 `filter_peaks_by_adducts(...)` 写出后 shared m/z 轴减少。不要在第一版计划里跑完整数据集。

### Task 8: 静态检查

**Files:**
- No source edits unless checks reveal a real issue in files touched above.

- [ ] **Step 1: 只做 py_compile**

```bash
.venv/bin/python -m py_compile \
  src/massflow/preprocess/helper/adduct_filter_helper.py \
  src/massflow/preprocess/api.py \
  src/massflow/preprocess/pipeline/planner.py \
  tests/test_adduct_filter.py
```

Expected: no output and exit code 0.

- [ ] **Step 2: 检查 whitespace**

```bash
git diff --check
```

Expected: no output and exit code 0.

- [ ] **Step 3: 不默认运行 pytest**

如果用户明确允许 targeted pytest，再运行：

```bash
.venv/bin/python -m pytest tests/test_adduct_filter.py -q
```

## 审核重点

- 算法必须是 CAMERA/MZmine-style neutral mass hypothesis，不是简单 `mz + delta` 删除。
- 默认规则表必须注明来源于 CAMERA primary adduct rules，并用测试锁住公式和数值。
- 任何未能溯源到 CAMERA/MZmine 的策略必须标注为 MassFlow/MSI adaptation；不能写成来源算法。
- 必须有 annotation/report table，否则无法确认 removed peaks 与 source-backed rule match 对应关系。
- adduct filter 第一版应是 dataset-scope，放在 `peak_align` 之后；不要接成 isotope filter 那样的 per-spectrum batch task。
- 需要有 spatial profile correlation 或等价的 group confirmation；没有确认步骤会偏离 CAMERA/MZmine 的核心前提。
- negative controls 必须证明“只有 m/z rule match 不够删除”，还需要 spatial profile confirmation。
- 不要创建全 feature x feature correlation 矩阵；只对 neutral-mass candidate groups 内的 peaks 取 profile 做相关性。
- 不要运行完整 pytest。

## 参考来源

- Bioconductor CAMERA package page: https://bioconductor.org/packages/release/bioc/html/CAMERA.html
- CAMERA `findAdducts` documentation: https://rdrr.io/bioc/CAMERA/man/findAdducts-methods.html
- CAMERA source package inspected: `CAMERA_1.68.0.tar.gz`, especially `R/fct_findAdducts.R` and `inst/rules/primary_adducts_*.csv`
- MZmine Ion Networking source: https://github.com/mzmine/mzmine/blob/master/mzmine-community/src/main/java/io/github/mzmine/modules/dataprocessing/id_ion_identity_networking/ionidnetworking/IonNetworkingTask.java
- MZmine ion type formula source: https://github.com/mzmine/mzmine/blob/master/mzmine-community/src/main/java/io/github/mzmine/datamodel/identities/iontype/IonType.java
- MZmine CliqueMS adduct source: https://github.com/mzmine/mzmine/blob/master/mzmine-community/src/main/java/io/github/mzmine/modules/dataprocessing/id_cliquems/cliquemsimplementation/AdductAnnotationCliqueMS.java
