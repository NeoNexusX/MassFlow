# MassFlow MSI Zarr 格式规范

版本：`1.1`

## 1. 文件结构

v1.1 使用两个并列的布局数据组：`/ion_image`（ion-major）与 `/spectra`（pixel-major）。
一份文件按写入模式取以下三种合法结构之一：

```text
# 1) ion-only（IonImageDataManagerZarr 默认）
sample.zarr/
├── metadata/                       # 两种布局共用
├── axes/
│   ├── coordinates                 # (n_pixels, 3), <u4，像素坐标 [x, y, z]
│   └── mz                          # (n_ions,)，共享 m/z 轴
├── ion_image/
│   ├── .zattrs                     # row_axis="ion", encoding="continuous"
│   ├── intensity                   # ion-major，shape=(n_ions * n_pixels,)
│   └── offsets                     # [0, n_pixels, 2*n_pixels, ...]
└── stats/
    └── mean_spectrum               # 可选，source="/ion_image"
```

```text
# 2) spectra-only（MSDataManagerZarr）
sample.zarr/
├── metadata/                       # 两种布局共用
├── axes/
│   ├── coordinates                 # (n_pixels, 3), <u4，像素坐标 [x, y, z]
│   └── mz                          # continuous 模式必选；processed 模式不存在
├── spectra/
│   ├── .zattrs                     # row_axis="pixel", encoding=continuous|processed
│   ├── intensity                   # pixel-major，shape=(n_pixels * n_mz,) 或 (sum(lengths),)
│   ├── offsets                     # [0, n_mz, 2*n_mz, ...] 或 cumulative lengths
│   └── mz                          # 仅 processed 模式，逐点 m/z
└── stats/
    └── tic                         # 可选，source="/spectra"
```

```text
# 3) continuous combined（IonImageDataManagerZarr(include_spectra=True)）
sample.zarr/
├── metadata/                       # 两种布局共用
├── axes/                           # 两种布局共用
│   ├── coordinates                 # (n_pixels, 3)，顺序与 /spectra 的行一致
│   └── mz                          # (n_ions,)，continuous 共享 m/z 轴
├── ion_image/
│   ├── .zattrs                     # row_axis="ion", encoding="continuous"
│   ├── intensity                   # ion-major，shape=(n_ions * n_pixels,)
│   └── offsets                     # [0, n_pixels, 2*n_pixels, ...]
├── spectra/
│   ├── .zattrs                     # row_axis="pixel", encoding="continuous"
│   ├── intensity                   # pixel-major，shape=(n_pixels * n_ions,)
│   └── offsets                     # [0, n_ions, 2*n_ions, ...]
└── stats/
    ├── tic                         # (n_pixels,)，source="/spectra"
    └── mean_spectrum               # (n_ions,)，source="/ion_image"
```

v1.1 不再有 `/data` 组。结构判定只依据合法组合：`/ion_image`、`/spectra`
至少存在一个，且 `/data` 必须不存在。v1.0 明确不受支持，也不提供兼容读取。

---

## 2. 根属性

所有根属性均为必需，位于 `/.zattrs`。

```json
{
  "format": "massflow.msi_zarr",
  "format_version": "1.1",
  "write_state": "complete",
  "coordinate_order": ["x", "y", "z"],
  "coordinate_base": 1,
  "spatial_shape": [height, width]
}
```

| 属性 | 类型 | 含义 |
| --- | --- | --- |
| `format` | string | 固定为 `massflow.msi_zarr`，标识 MassFlow MSI Zarr 格式 |
| `format_version` | string | 该文件对应的 schema 版本，固定为 `1.1` |
| `write_state` | string | `"writing"` 或 `"complete"`；只有 `"complete"` 可被 reader 打开 |
| `coordinate_order` | string[] | 坐标含义顺序，固定为 `["x", "y", "z"]` |
| `coordinate_base` | integer | 坐标原点基准，固定为 `1`（one-based） |
| `spatial_shape` | integer[] | 空间维度 `[height, width]`，单位为像素 |

`format_version` 指 MassFlow MSI Zarr schema 版本，不是 Zarr 库的存储格式版本。

---

## 3. `axes/` 规范

`axes/` 被两份布局数据组共用，不重复保存。`metadata/` 同样共用。

### 3.1 `axes/coordinates`

```text
shape = (n_pixels, 3)
dtype = <u4
```

约束：

- 所有坐标为 one-based，与 `coordinate_base` 一致。
- 每行对应一个采集像素，顺序与 `/spectra` 的行严格对应；在 ion-only 文件中
  与 `/ion_image` 中每个 ion 向量内的第 `p` 个值对应。
- ion-major 写入时自动去重；坐标重复时保留最后写入的像素并记录 warning。
  combined 文件中 `/spectra` 与 `/ion_image` 使用完全相同的有效像素集合与顺序。
- 二维数据 `z` 固定为 `1`。
- `1 ≤ x ≤ width` 且 `1 ≤ y ≤ height`。

坐标顺序和基准只记录在根属性，不在数组 attrs 中重复。

### 3.2 `axes/mz`

```text
shape = (n_mz,)
dtype = <f8 或 <f4
attrs.unit = "m/z"
```

存在条件：

- **pixel-major continuous**：必须存在；
- **ion-major continuous**：必须存在；
- **combined continuous**：必须存在，同时作为 `/ion_image` 的 ion 轴与
  `/spectra` 每行谱的列轴；
- **pixel-major processed**：不存在（逐点 m/z 位于 `/spectra/mz`）。

---

## 4. 布局数据组规范

### 4.1 布局组属性

每个布局数据组的 `.zattrs` 均为必需。

`/ion_image/.zattrs`：

```json
{
  "row_axis": "ion",
  "encoding": "continuous"
}
```

`/spectra/.zattrs`：

```json
{
  "row_axis": "pixel",
  "encoding": "continuous"
}
```

或 `"encoding": "processed"`。

| 属性 | 允许值 | 含义 |
| --- | --- | --- |
| `row_axis` | `"pixel"` / `"ion"` | 一维数组按哪条逻辑轴组织为行 |
| `encoding` | `"continuous"` / `"processed"` | m/z 编码方式 |

约束：

- `row_axis="ion"` 必须配合 `encoding="continuous"`，且只能出现在 `/ion_image`；
- `row_axis="pixel"` 允许 `continuous` 或 `processed`，且只能出现在 `/spectra`；
- **combined 文件只支持 continuous**；spectra-only 仍支持 continuous/processed。

不保存 `n_rows`、`total_points`、`intensity_dtype`、`mz_dtype` 等可从数组
shape/dtype 直接推导的属性，避免形成不一致的第二事实来源。

### 4.2 `intensity` 与 `offsets`

两份布局共用同一套强度数组约定，只是行组织不同：

- `intensity`：一维，按 `offsets` 划分为 `n_rows` 行。默认 dtype 为 `float32`
  （写入）；读取接受任何浮点类型。
- `offsets`：`(n_rows + 1,)`，`<i8`，满足
  `offsets[0] == 0`、`offsets[i+1] >= offsets[i]`、`offsets[-1] == len(intensity)`。

即使所有行等长，`offsets` 也必须落盘，使 reader 永远使用同一条行切片路径。
空行由相等的相邻 offset 表达；reader 可读空行，writer 默认拒绝写入空谱。

#### `/ion_image`

```text
intensity.shape = (n_ions * n_pixels,)
offsets         = [0, n_pixels, 2*n_pixels, ..., n_ions*n_pixels]
```

第 `i` 个 ion vector：

```python
start, end = ion_offsets[i : i + 2]
mz = axes_mz[i]
values = ion_intensity[start:end]
```

`values[j]` 对应 `axes/coordinates[j]`，而不是 `(H, W)` 中的第 `j` 个稠密位置。
生成 ion image：

```python
image = np.zeros((height, width), dtype=values.dtype)
x = coordinates[:, 0]
y = coordinates[:, 1]
image[y - 1, x - 1] = values
```

#### `/spectra`（continuous）

```text
intensity.shape = (n_pixels * n_ions,)
offsets         = [0, n_ions, 2*n_ions, ..., n_pixels*n_ions]
```

第 `p` 个像素谱：

```python
start, end = spectra_offsets[p : p + 2]
mz = axes_mz[:]
intensity = spectra_intensity[start:end]
coordinate = axes_coordinates[p]
```

#### `/spectra`（processed）

```text
intensity.shape = (sum(lengths),)
mz.shape        = (sum(lengths),)
offsets         = cumulative lengths
```

第 `p` 条谱：

```python
start, end = spectra_offsets[p : p + 2]
mz = spectra_mz[start:end]
intensity = spectra_intensity[start:end]
```

`spectra/mz` 的每一行切片与对应 intensity 切片长度相同。

---

## 5. 三种合法形态

| 形态 | `/ion_image` | `/spectra` | `encoding` |
| --- | --- | --- | --- |
| ion-only | ✓ | — | continuous |
| spectra-only | — | ✓ | continuous 或 processed |
| combined | ✓ | ✓ | continuous（两份） |

判定规则：`/ion_image`、`/spectra` 至少存在一个，`/data` 必须不存在；不再识别
v1.0 文件。

---

## 6. `metadata/` 规范

### 6.1 `metadata/.zattrs`

所有标量元数据字段，由 `MetaDataZarr` 写入。以下是标准化字段名：

| 字段 | 类型 | 含义 |
| --- | --- | --- |
| `name` | string | 数据集名称 |
| `version` | number | 元数据版本号 |
| `spectrum_count_num` | integer | 谱总数（有效像素数） |
| `max_count_of_pixels_x` | integer | 空间宽度 |
| `max_count_of_pixels_y` | integer | 空间高度 |
| `pixel_size_x` | number | X 方向像素大小（μm） |
| `pixel_size_y` | number | Y 方向像素大小（μm） |
| `absolute_position_offset_x` | number | 绝对值坐标偏移 X |
| `absolute_position_offset_y` | number | 绝对值坐标偏移 Y |
| `instrument_model` | string | 仪器型号 |
| `centroid_spectrum` | boolean | 是否为 Centroided 谱 |
| `profile_spectrum` | boolean | 是否为 Profile 谱 |
| `continuous` | boolean | continuous 编码标志 |
| `processed` | boolean | processed 编码标志 |
| `ms1_spectrum` | boolean | 含 MS1 谱 |
| `msn_spectrum` | boolean | 含 MSn 谱 |

其他源 metadata 字段会原样保留。

### 6.2 `metadata/mask`

可选，存在时为二维二值数组，shape 为 `(height, width)`，dtype 由 Zarr 序列化
决定。值为 1 的位置表示有实际采集数据，0 表示空洞或未采集区域。

---

## 7. `stats/` 规范

v1.1 中 `tic` 与 `mean_spectrum` 不再彼此互斥：

- **combined 文件二者必选**，即使调用方没有显式打开统计量开关；
- **单布局文件二者按各自开关可选**：ion-only 的 `mean_spectrum` 可选、
  spectra-only 的 `tic` 可选。

### 7.1 `stats/mean_spectrum`

```text
shape = (n_ions,)
dtype = float64
attrs.source = "/ion_image"
attrs.reduced_axis = "pixel"
attrs.denominator = "sampled_pixels"
```

定义：

```text
mean_spectrum[i] = sum(ion_image 第 i 行) / n_pixels
```

分母是 `len(axes/coordinates)`，不是 `height * width`。如果需要包含空间空洞的
平均值，必须用不同的 `denominator` 明确标识。

### 7.2 `stats/tic`

```text
shape = (n_pixels,)
dtype = float64
attrs.source = "/spectra"
attrs.reduced_axis = "mz"
```

定义：

```text
tic[p] = sum(spectra 的第 p 行)
```

第 `p` 个元素对应 `axes/coordinates[p]`，顺序与 `/spectra/offsets` 的行顺序
严格一致。TIC 是求和而非平均，不挂 `denominator` 属性。

continuous 与 processed 编码均可生成；TIC 只依赖每条谱的 intensity 切片，不需要
共享 m/z 网格。combined 中 TIC 在批写入阶段按 pixel-major 行直接求和累积，
finalize 时不重读任何一份完整 intensity。

坐标去重时（重复坐标保留最后一条谱），`tic` 保留最后一条谱的 TIC；`mean_spectrum`
同样基于 compact 后的有效像素计算。

---

## 8. 读取约定

- ion-only 文件将 `/ion_image` 作为主布局；`row_axis` 返回 `"ion"`。
- spectra-only 文件将 `/spectra` 作为主布局；`row_axis` 返回 `"pixel"`。
- combined 文件将 `/ion_image` 作为 Python reader 的主布局；`/spectra` 只作为
  前端直接读取的 schema 契约，reader 不把它绑定到任何 pixel-major 读取 API。
- reader 不接受 `format_version="1.0"`，不查找 `/data`，无 fallback/迁移/别名。
