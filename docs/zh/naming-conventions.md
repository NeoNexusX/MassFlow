# MassFlow 命名习惯与规范

本文以当前源码中的类、模块和存储布局为准。

## Python 命名

- 模块、函数、参数和变量使用 `snake_case`：`ms_data_manager_imzml.py`、`load_head_data`、`write_flat_batch`、`mz_axis`。
- 类使用 `PascalCase`，保留已有缩写：`MSDataManagerImzML`、`MSDataManagerZarr`、`IonImageDataManagerH5`、`MetaDataImzMl`。
- 常量使用 `ALL_CAPS`：`ZARR_MODE_CONTINUOUS`、`ZARR_FORMAT_VERSION`。
- 内部状态和 helper 使用单下划线前缀：`_reader`、`_queue`、`_require_reader`。
- 使用公开导出路径，例如 `from massflow.data_manager import MSDataManagerImzML`。

## 类与职责

| 职责 | 当前类名 |
| --- | --- |
| 单条谱 | `Spectrum`、`SpectrumImzML`、`SpectrumZarr` |
| 谱集合 | `MassSpectrumSet` |
| 像素坐标 | `PixelCoordinates` |
| 谱数据管理 | `MSDataManager`、`MSDataManagerImzML`、`MSDataManagerZarr` |
| 离子图管理 | `IonImageDataManager`、`IonImageDataManagerBinary`、`IonImageDataManagerH5`、`IonImageDataManagerZarr` |
| 元数据 | `MetaDataBase`、`MetaDataImzMl`、`MetaDataZarr` |
| 预处理 | `Preprocessor`、`PreprocessorAPI`、`FlatPreprocess` |

源码按职责分为 `module/`、`data_manager/`、`preprocess/`、`msi_zarr/`、`segmentation/` 和 `tools/` 等目录。

## 操作与参数

操作名使用当前 API 中的动词形式：

- `load_head_data()`：加载元数据和惰性谱占位对象。
- `extract_metadata()`：将元数据读入管理器的模型。
- `flat_generator()`、`matrix_generator()`：生成批次数据。
- `write_flat_batch()`：写入扁平批次；写出器还提供 `add_spectrum()`。
- `read_ion_image()` / `plot_ion_image()`：通过可视化 helper 读取或显示离子图。
- `get_ion_image()`：离子图管理器的读取操作。
- `finalize()`：完成支持该操作的写出器；`close()`：释放资源。

统一使用 `mz`、`intensity`、`lengths`、`coordinates`、`mz_axis` 表达数据含义。跨 API 复制示例前先核对签名，例如底层 Zarr 写出器使用 `add_spectrum(mz, intensity, coordinates)`。

## 元数据字段

公开字段不带 `meta_` 前缀：`meta.pixel_size_x`、`meta.continuous`、`meta.profile_spectrum`。描述符和属性会把值同步到 `_meta` 字典。打印结果在标签前加 `meta_`，打印标签不是 Python 属性名。

通过 `to_dict()`、`get()`、`items()`、`keys()`、`values()` 读取元数据。类名和字段含义见[元数据](./metadata.md)。

`MassSpectrumSet` 提供 `add_spectrum()`、`update_spectrum_with_coord()`、`update_spectrum_with_index()`、迭代、索引和 `len()`；队列属于内部实现。

## 存储命名

- Zarr 组名为 `axes`、`metadata`、`spectra`、`ion_image`、`stats`；数组名包括 `intensity`、`offsets`、`coordinates`、`mz`，详见 [Zarr 格式规范](./zarr_schema.md)。
- H5 离子图的数据集为 `mz_axis`、`coordinates`、`ion_images`；维度和计数保存在根属性中。
- 编码名为 `continuous`、`processed`；写出器可接受 `auto` 推断。编码表达 m/z 的组织方式。
- 支持读写模式的管理器使用 `mode="r"` 或 `mode="w"`，不要把编码名与文件访问模式混用。

## 示例与校验

示例使用 `source_dm`、`processed_dm`、`mz_axis`、`coordinates` 等明确名称，显式给出文件路径，默认参数以源码签名为准。

公开输入和 I/O 边界使用异常表达错误，断言表达内部不变量。扩展 API 时沿用已有命名，避免为同一操作新增平行别名。
