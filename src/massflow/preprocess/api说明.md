**Centroid = 峰已经被提取成离散 m/z-intensity 点
 Profile = 原始连续谱**

# Preprocessor API 参数参考

**需要用户传参数：**

1. ## **`peak_pick` — 峰检测（batch 或 dataset 范围）**

| 参数          | 类型      | 默认值     | 可传值 / 取值范围                                            |
| ------------- | --------- | ---------- | ------------------------------------------------------------ |
| `method`      | `str`     | `"diff"`   | python 后端：`"sd"`、`"mad"`、`"quantile"`、`"diff" `                                            cardinal 后端：`"diff"`、`"sd"`、`"mad"`、`"quantile"`、`"filter"`、`"cwt"` |
| `snr`         | `float`   | `2.0`      | 正数（≤0 或 NaN 时跳过 SNR 过滤）                            |
| `return_type` | `str`     | `"height"` | `"height"` 或 `"area"`                                       |
| `width`       | `int`     | `5`        | ≥1 整数（小于 1 会被钳为 1）如果减少一点参数，该参数可以不要 |
| `backend`     | `Backend` | `"python"` | `"python"`（batch）或 `"cardinal"`（dataset，调 R Cardinal） |

---

## 2. `normalization` — 归一化（batch 范围）

| 参数            | 类型                 | 默认值        | 可传值 / 取值范围                                  |
| --------------- | -------------------- | ------------- | -------------------------------------------------- |
| `method`        | `str`                | `"tic_numba"` | **仅** `"tic_numba"`、`"rms_numba"`、`"ref_numba"` |
| `scale`         | `float \| None`      | `None`        | 非负有限数或 `None`（输出整体缩放因子）一般为1     |
| `mz_flat`       | `np.ndarray \| None` | `None`        | 1D ndarray，**`ref_numba` 必填**                   |
| `ref`           | `float \| None`      | `None`        | float，**`ref_numba` 必填**（参考 m/z）            |
| `ref_tolerance` | `float`              | `0.1`         | 非负有限数（参考 m/z 容差，`ref_numba` 用）        |

## 3. `peak_align` — 峰对齐（dataset 范围）

| 参数        | 类型            | 默认值     | 可传值 / 取值范围                                            |
| ----------- | --------------- | ---------- | ------------------------------------------------------------ |
| `tolerance` | `float \| None` | `None`     | 正数或 `None`（None 时按 `binratio` 自动估计；`units="ppm"` 时单位为 ppm） |
| `units`     | `str`           | `"ppm"`    | `"ppm"` / `"relative"`（同义）或 `"Da"` / `"absolute"`（同义）（[peak_align_helper.py:26-28](helper/peak_align_helper.py#L26-L28)） |
| `binfun`    | `str`           | `"median"` | `"median"`、`"min"`、`"max"`、`"mean"`（[peak_align_helper.py:115-121](helper/peak_align_helper.py#L115-L121)；其它字符串会回退为 mean） |
| `binratio`  | `float`         | `2.0`      | 正数（容差/分辨率比例）                                      |
| `backend`   | `Backend`       | `"python"` | `"python"` 或 `"cardinal"`（cardinal 后端不支持 `min_count` / `min_frequency`，传入会 raise，[api.py:156-157](api.py#L156-L157)） |

## 4. `noise_reduction` — 平滑降噪（batch 范围）

| 参数        | 类型            | 默认值             | 可传值 / 取值范围                                            |
| ----------- | --------------- | ------------------ | ------------------------------------------------------------ |
| `method`    | `str`           | `"gaussian_numba"` | **仅** `"savgol_numba"`、`"gaussian_numba"`、`"ma_numba"`（[flat_pre_fun.py:81-83](flat_pre_fun.py#L81-L83)） |
| `window`    | `int`           | `5`                | 正整数                                                       |
| `sd`        | `float \| None` | `None`             | 正数或 `None`（gaussian_numba 用：高斯标准差）               |
| `polyorder` | `int`           | `3`                | 非负整数（savgol_numba 用：多项式阶）                        |
| `deriv`     | `int`           | `0`                | 非负整数（savgol_numba 用：导数阶，0=平滑）                  |
| `delta`     | `float`         | `1.0`              | 正数（savgol_numba 用：采样间距）                            |



**不需要用户传参数：**

##  5.`baseline_correction` — 基线校正（batch 范围）

| 参数             | 类型          | 默认值           | 可传值 / 取值范围                                   |
| ---------------- | ------------- | ---------------- | --------------------------------------------------- |
| `method`         | `str`         | `"locmin_numba"` | **仅** `"locmin_numba"` 或 `"snip_numba"`           |
| `smooth`         | `str`         | `"none"`         | `"none"` 或 `"loess"`（仅 locmin_numba 使用）       |
| `span`           | `float`       | `0.1`            | `(0, 1]`，loess 平滑窗占比（locmin_numba 用）       |
| `upper`          | `bool`        | `False`          | `True` / `False`（locmin_numba 用，反向 envelope）  |
| `width`          | `int`         | `5`              | 正整数，局部极值窗口（locmin_numba 用）             |
| `niter`          | `int`         | `15`             | 正整数，迭代次数（locmin_numba 用）                 |
| `baseline_scale` | `float`       | `1.0`            | `(0, 1]`，基线扣除强度缩放                          |
| `m`              | `int \| None` | `None`           | 正整数或 `None`（snip_numba 的滤波半宽，None 自动） |
| `decreasing`     | `bool`        | `True`           | `True` / `False`（snip_numba：递减/递增窗策略）     |

---

## **Centroid数据一般只做2和3**

---

## 关键约束速查

- **后端互斥参数**：`peak_align` 的 `min_count` / `min_frequency` 只能 python 后端用；`peak_pick` 的 `relheight` 只能 python 后端用。
- **`ref_numba` 归一化必填**：`mz_flat` 和 `ref` 缺一不可。
- **所有方法名都是字符串字面量**，传入会 `.strip().lower()` 后校验，错值会抛 `ValueError`。
- **scope**：`baseline_correction` / `noise_reduction` / `normalization` 是 batch 范围；`peak_align` 始终是 dataset 范围；`peak_pick` 在 python 下是 batch、cardinal 下是 dataset。
