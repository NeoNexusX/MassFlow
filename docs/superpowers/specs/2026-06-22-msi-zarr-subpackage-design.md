# MSI Zarr 子包设计

## 目标

将当前位于 `massflow.tools` 的 MSI Zarr reader、writer 和 schema 定义移动到专用子包 `massflow.msi_zarr`，使格式实现拥有清晰边界。

本次仅调整模块组织和 import 路径，不改变：

- Zarr schema；
- reader、writer 的运行行为；
- Binary、H5、Zarr 的读写和性能；
- data manager、绘图、分割及预处理行为；
- 公共类名和常量名。

## 目录结构

```text
src/massflow/msi_zarr/
├── __init__.py
├── schema.py
├── reader.py
└── writer.py
```

### `schema.py`

集中定义 reader 和 writer 共享的格式契约：

- `ZARR_FORMAT_NAME`
- `ZARR_FORMAT_VERSION`
- `ZARR_MODE_CONTINUOUS`
- `ZARR_MODE_PROCESSED`
- `ZARR_MODE_AUTO`
- `ZarrStorageMode`
- `ZarrWriteMode`
- `ZarrRowAxis`

reader 和 writer 均依赖 `schema.py`，消除 reader 从 writer 导入 schema 定义的反向依赖。

### `reader.py`

容纳当前 `massflow.tools.zarr_reader` 的 `MSIZarrReader` 及 reader 私有辅助函数。

reader 只负责：

- 打开和验证 MassFlow MSI Zarr store；
- 暴露 metadata、维度和 dtype；
- 读取 spectrum、ion row、batch 和统计量；
- 校验已落盘结构。

### `writer.py`

容纳当前 `massflow.tools.zarr_writer` 的 `MSIZarrWriter`、`CoordinateIndex` 及 writer 私有辅助函数。

writer 继续负责：

- pixel-major 和 ion-major 写入；
- 坐标验证与去重；
- chunk 计算；
- continuous/processed 编码处理；
- finalize 和写后结构校验。

当前坐标逻辑只由 writer 使用，因此不单独创建 `coordinates.py`。

### `__init__.py`

统一导出子包公共 API：

```python
from massflow.msi_zarr import (
    MSIZarrReader,
    MSIZarrWriter,
    ZARR_MODE_AUTO,
    ZARR_MODE_CONTINUOUS,
    ZARR_MODE_PROCESSED,
    ZarrRowAxis,
    ZarrStorageMode,
    ZarrWriteMode,
)
```

writer 私有辅助函数不从 `__init__.py` 导出。当前测试直接使用的 `prepare_coordinate_index` 从 `massflow.msi_zarr.writer` 导入。

## 保留位置

以下内容本次不移动：

- `MetaDataZarr`：仍位于 `massflow.module.ms_meta_data_zarr`；
- `IonMajorTempStore`：仍位于 `massflow.tools.ion_major_temp_store`；
- `MSDataManagerZarr` 和 `IonImageDataManagerZarr`：仍位于 `massflow.data_manager`；
- schema 文档：仍位于 `docs/zh/zarr_schema.md`。

这样可以避免把纯模块整理扩大为 data manager 或 metadata 重构。

## Import 迁移

仓库内消费者统一改为从 `massflow.msi_zarr` 导入公共 API：

- `MSDataManagerZarr`
- `IonImageDataManagerZarr`
- `IonImageZarrDataset`
- 相关测试

不保留 `massflow.tools.zarr_reader` 和 `massflow.tools.zarr_writer` 兼容模块。它们属于当前尚未稳定发布的重构内容，保留转发文件只会形成两套入口。

`massflow.data_manager` 现有的 Zarr mode 常量导出保持不变，避免影响从 data manager 顶层导入这些常量的调用方。

## 依赖方向

```text
massflow.data_manager ───────┐
massflow.morch ──────────────┼──> massflow.msi_zarr
tests ───────────────────────┘

massflow.msi_zarr.reader ────> massflow.msi_zarr.schema
massflow.msi_zarr.writer ────> massflow.msi_zarr.schema
massflow.msi_zarr.reader ────> MetaDataZarr
massflow.msi_zarr.writer ────> MetaDataZarr, IonMajorTempStore
```

`massflow.msi_zarr` 不依赖任何 data manager，防止循环依赖。

## 验证

按当前仓库要求不运行 pytest。完成迁移后执行：

1. `compileall` 检查 `src` 和 `tests`；
2. `git diff --check`；
3. 全仓搜索，确认不存在旧 `massflow.tools.zarr_reader`、`massflow.tools.zarr_writer` import；
4. continuous 和 processed spectra Zarr round-trip；
5. ion-major Zarr 写出、重开、ion image/vector/block 读取；
6. Binary、H5、Zarr 写出与读取性能对比；
7. imzML、Binary、H5、Zarr 的 `plot_ion_image`、`plot_tic_image`、`plot_umap_kmeans_image` 冒烟验证。

## 非目标

本次不处理此前审核发现的：

- Zarr chunk 布局性能问题；
- `coordinate_base` 兼容性；
- AUTO 模式单谱编码判断；
- metadata 扩展字段透传；
- `apply_region_mask` 性能回归。

这些问题需要独立修改和独立验证，不能混入本次纯目录重构。
