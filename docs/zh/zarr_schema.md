# MassFlow MSI Zarr 格式规范

版本：`1.0`

## 1. 文件结构

```text
sample.zarr/
├── metadata/                       # 实验、仪器、谱类型等语义元数据
│   ├── .zattrs                     # 所有标量元数据字段
│   └── mask                        # 可选，二维二值空间掩膜
│
├── axes/
│   ├── coordinates                 # (n_pixels, 3), <u4，像素坐标 [x, y, z]
│   └── mz                          # 可选，共享 m/z 轴，continuous 模式必选
│
├── data/
│   ├── intensity                   # 一维，所有谱/离子的强度数据
│   ├── offsets                     # (n_rows + 1), <i8，行起始位置索引
│   └── mz                          # 可选，processed 模式的逐点 m/z
│
└── stats/                          # 可选，可删除、可重算的统计量
    └── mean_spectrum
```

`data/` 是唯一规范数据位置。`stats/` 及其子数组不是正确读取文件所必需的。

---

## 2. 根属性

所有根属性均为必需，位于 `/.zattrs`。

```json
{
  "format": "massflow.msi_zarr",
  "format_version": "1.0",
  "write_state": "complete",
  "coordinate_order": ["x", "y", "z"],
  "coordinate_base": 1,
  "spatial_shape": [height, width]
}
```

| 属性 | 类型 | 含义 |
| --- | --- | --- |
| `format` | string | 固定为 `massflow.msi_zarr`，标识 MassFlow MSI Zarr 格式 |
| `format_version` | string | 该文件对应的 schema 版本，固定为 `1.0` |
| `write_state` | string | `"writing"` 或 `"complete"`；只有 `"complete"` 可被 reader 打开 |
| `coordinate_order` | string[] | 坐标含义顺序，固定为 `["x", "y", "z"]` |
| `coordinate_base` | integer | 坐标原点基准，固定为 `1`（one-based） |
| `spatial_shape` | integer[] | 空间维度 `[height, width]`，单位为像素 |

`format_version` 指 MassFlow MSI Zarr schema 版本，不是 Zarr 库的存储格式版本。

---

## 3. `axes/` 规范

### 3.1 `axes/coordinates`

```text
shape = (n_pixels, 3)
dtype = <u4
```

约束：

- 所有坐标为 one-based，与 `coordinate_base` 一致。
- 每行对应一个采集像素，顺序与 `data/` 谱顺序严格对应。
- 写入时自动去重；坐标重复时保留最后写入的像素并记录 warning。
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
- **pixel-major processed**：不存在。

`axes/mz` 必须严格对应 ion 行编号或 continuous 谱中的列编号。

---

## 4. `data/` 规范

### 4.1 `data` 属性

所有 data 属性均为必需，位于 `/data/.zattrs`。

```json
{
  "row_axis": "pixel",
  "encoding": "continuous"
}
```

| 属性 | 允许值 | 含义 |
| --- | --- | --- |
| `row_axis` | `"pixel"` / `"ion"` | 一维数组按哪条逻辑轴组织为行 |
| `encoding` | `"continuous"` / `"processed"` | m/z 编码方式 |

不保存 `n_rows`、`total_points`、`intensity_dtype`、`mz_dtype` 等可从数组 shape/dtype 直接推导的属性，避免形成不一致的第二事实来源。

`encoding` 不只依赖节点存在性推断。reader 应同时检查属性与结构一致性：

- `encoding="continuous"` 时必须有 `axes/mz` 且不得存在 `data/mz`；
- `encoding="processed"` 时必须有 `data/mz` 且不得存在 `axes/mz`；
- `row_axis="pixel"` 允许 `continuous` 或 `processed`；
- `row_axis="ion"` 必须配合 `encoding="continuous"`。

### 4.2 `data/intensity`

所有规范强度值组成的一维数组。逻辑上按 `offsets` 划分为 `n_rows` 行。

实现要求：

- 默认 dtype：`float32`（写入）；读取接受任何浮点类型。
- chunk 是一维 point chunk，以目标字节数控制上限。
- compressor 和 filter 不写入格式规范，由 writer 配置决定。

### 4.3 `data/offsets`

```text
shape = (n_rows + 1,)
dtype = <i8
```

必须满足：

```text
offsets[0] == 0
offsets[i + 1] >= offsets[i]
offsets[-1] == len(data/intensity)
```

即使所有行等长，也必须落盘。这样 reader 永远使用同一条行切片路径。

空行在结构上由相等的相邻 offset 表达。reader 必须能够读取空行；writer 默认拒绝写入空谱。

### 4.4 `data/mz`

仅用于 pixel-major processed：

```text
shape = data/intensity.shape
dtype = <f8 或 <f4
attrs.unit = "m/z"
```

每一行的 m/z 切片与对应 intensity 切片长度相同。

---

## 5. 三种规范形态

### 5.1 Pixel-major Continuous

```text
data.attrs.row_axis = "pixel"
data.attrs.encoding = "continuous"

axes/coordinates.shape = (n_pixels, 3)
axes/mz.shape          = (n_mz,)
data/intensity.shape   = (n_pixels * n_mz,)
data/offsets           = [0, n_mz, 2*n_mz, ...]
data/mz                = 不存在
```

第 `p` 条谱：

```python
start, end = offsets[p : p + 2]
mz = axes_mz[:]
intensity = data_intensity[start:end]
```

### 5.2 Pixel-major Processed

```text
data.attrs.row_axis = "pixel"
data.attrs.encoding = "processed"

axes/coordinates.shape = (n_pixels, 3)
axes/mz                 = 不存在
data/intensity.shape    = (sum(lengths),)
data/mz.shape           = (sum(lengths),)
data/offsets             = cumulative lengths
```

第 `p` 条谱：

```python
start, end = offsets[p : p + 2]
mz = data_mz[start:end]
intensity = data_intensity[start:end]
```

### 5.3 Ion-major Continuous

```text
data.attrs.row_axis = "ion"
data.attrs.encoding = "continuous"

axes/coordinates.shape = (n_pixels, 3)
axes/mz.shape          = (n_ions,)
data/intensity.shape   = (n_ions * n_pixels,)
data/offsets           = [0, n_pixels, 2*n_pixels, ...]
data/mz                = 不存在
```

第 `i` 个 ion vector：

```python
start, end = offsets[i : i + 2]
mz = axes_mz[i]
values = data_intensity[start:end]
```

`values[j]` 对应 `axes/coordinates[j]`，而不是 `(H, W)` 中的第 `j` 个稠密位置。

生成 ion image：

```python
image = np.zeros((height, width), dtype=values.dtype)
x = coordinates[:, 0]
y = coordinates[:, 1]
image[y - 1, x - 1] = values
```

这种表示只保存实际采集像素，不为未采集位置写零。

---

## 6. `metadata/` 规范

### 6.1 `metadata/.zattrs`

所有标量元数据字段，由 `MetaDataZarr` 写入。以下是标准化字段名：

| 字段 | 类型 | 含义 |
| --- | --- | --- |
| `name` | string | 数据集名称 |
| `version` | number | 元数据版本号 |
| `spectrum_count_num` | integer | 谱总数 |
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

可选，存在时为二维二值数组，shape 为 `(height, width)`，dtype 由 Zarr 序列化决定。值为 1 的位置表示有实际采集数据，0 表示空洞或未采集区域。

---

## 7. `stats/` 规范

### 7.1 `stats/mean_spectrum`

仅当写入时指定 `compute_mean_spectrum=True` 时存在。

```text
shape = (n_mz,)
dtype = float64
attrs.reduced_axis = "pixel"
attrs.denominator = "sampled_pixels"
attrs.source = "/data"
```

定义：

```text
mean_spectrum[i] = sum(spectrum_p[i]) / n_pixels
```

分母是 `len(axes/coordinates)`，不是 `height * width`。如果需要包含空间空洞的平均值，必须用不同的 `denominator` 明确标识。

Processed spectra 没有共享 m/z 网格，除非先执行 binning/alignment，否则不能生成该数组。
