# MassFlow 命名习惯与规范


## 总则
- 统一采用“清晰、可读、含义明确”的命名，避免歧义与过度缩写。
- 保持风格一致：文件/函数/变量用 `snake_case`，类用 `PascalCase`。
- 专有缩写（如 `MSI`、`HDF5`、`ZYS`）按约定保留大写并嵌入类名中。
- 名称以职责为先，使用动词开头的方法名表达操作语义（如 `load_*`、`write_*`、`inspect_*`）。

## 文件与目录
- 目录使用小写单词，必要时用下划线分隔：`module/`。
- Python 文件使用 `snake_case`：`ms_module.py`、`ms_data_manager_imzml.py`、`msi_data_manager_msi.py`、`msi_data_manager_zys.py`。
- 示例文件：`example.ipynb`；不在示例中放业务逻辑。

## 模块与包
- 公开模块导出在 `module/__init__.py` 中维护（按需）。
- 相对导入保持简洁：`from .ms_module import MS, SpectrumBaseModule`。

## 类命名（PascalCase）
- 领域模型类：`MSI`。
- 组合/管理类：`MSIDataManager`、`MSIDataManagerZYS`。
- 基元/数据片类：`MSIBaseModule`。
- 规则：
  - 以领域名或职责名为主体。
  - 专有缩写保持大写嵌入，如 `MSI`，形成 `MSIDataManager` 而非 `MsiDataManager`。

## 方法与函数（snake_case）
- 动词开头，体现操作与对象：`load_full_data_from_file`、`load_meta_from_file`、`allocate_data_from_meta`、`add_msi_slice`、`get_msi_by_mz`、`plot_msi`、`inspect_data`、`write2local`。
- 辅助/内部方法可以下划线前缀标记：`_inspect_hdf5_structure`、`_write_meta_data`、`_create_datasets`。
- 带范围或过滤的参数以含义清晰的词组命名：`target_mz_range`、`display_threshold_percent`。

## 变量命名（snake_case）
- 普通变量：`mz`、`msroi`、`base_mask`、`coords`、`num_pixels`、`selected_channels`。
- 布尔变量以条件/状态语义命名：`need_base_mask`。
- 常量使用全大写加下划线：示例中为 `FILE_PATH`。
- 避免单字母除非循环/数学上下文（如 `i`）。

## 属性与元数据
- 私有元数据字段统一前缀：`_meta_*`，如 `_meta_name`、`_meta_version`、`_meta_mask`、`_meta_mz_num`、`_meta_storage_mode`、`_meta_need_base_mask`。
- 对应公开属性统一使用 `meta_*` 形式，并通过 `@property` 暴露：`meta_name`、`meta_version`、`meta_mask`、`meta_mz_num`、`meta_storage_mode`、`meta_need_base_mask`。
- 数据矩阵属性：公开 `data`，私有存储 `_data`。
- 更新元数据使用统一入口：`update_metadata()` 自动收集 `_meta_*` 字段。
- 队列访问统一：`add_msi_slice()`、`get_queue()`、`__len__()`、`__iter__()`、`__getitem__()`。

## 常量与枚举
- 使用全大写下划线：`FILE_PATH`、`FILE_PATH_ZYS`（示例中）。
- 模式字符串以小写明确：`mode in {"merge", "split"}`；写入前设置 `meta_storage_mode`。

## 缩写与专有名词
- `MSI`（Mass Spectrometry Imaging）：在类名/方法名中保持大写。
- `HDF5`：方法/变量中保持规范英文缩写（如 `_inspect_hdf5_structure`）。
- `ZYS`：来源数据的专有缩写可用于派生类名 `MSIDataManagerZYS` 与相关方法。
- `mz`、`msroi` 等领域术语保持与数据集一致，不做额外转换。

## 私有/内部实现约定
- 内部方法/字段使用 `_` 前缀：`_create_datasets`、`_write_meta_data`、`_metadata`、`_queue`。
- 仅通过公开方法/属性访问内部状态，避免直接访问私有字段。
- 断言与错误信息：
  - 使用 `assert` 与明确消息：`assert self.meta_mz_num > 0, "meta_mz_num must be greater than 0"`。
  - 参数/文件类型检查：`assert self.filepath.endswith('.mat'), "Error: filepath is not a .mat file."`。

## HDF5/数据集命名
- 组名：`mz_{mz_value:.5f}`（五位小数）。
- 数据集：
  - `mz`：对应单个 m/z 值。
  - `msroi`：二维图像矩阵（可 gzip 压缩）。
  - 元数据统一使用 `meta_*` 数据集名（与公开属性一致）。
- 元数据写入：
  - 公开键使用无下划线前缀形式：`meta_version`、`meta_mz_num`。
  - 类型：`version`/`num` 用数值类型（`np.float32`），字符串用 `h5py.string_dtype('utf-8')`。

## 参数与类型注解
- 类型注解按需添加，关键接口与返回值优先：`def get_msi(self) -> MSI`、`h5_data_zys: Optional[h5py.File]`。
- 默认值表达常见场景，提升可用性：`target_mz_range=None`、`storage_mode='split'`、`need_base_mask=False`。

## 示例一致性
- 示例中的实例/变量遵循相同风格：`msi_dm`、`msi_dm_zys`、`target_mz_range=[101, 102]`。
- 主入口使用惯例：
  ```python
  if __name__ == "__main__":
      ...
  ```

## 命名反例（避免）
- 混用大小写风格：如 `MsiDataManager`（错误，应为 `MSIDataManager`）。
- 文件使用驼峰或大写：如 `MsiModule.py`（错误，应为 `msi_module.py`）。
- 无语义的缩写：如 `tmp1`、`val2`（避免，使用语义化名称）。

## 落地建议
- 新增类/方法时先确定职责和数据语义，再落定命名。
- 领域缩写统一登记（如 `MSI`、`HDF5`、`ZYS`），在类名中保持一致。
- 涉及元数据的属性/数据集必须采用 `meta_*` 前缀，并维护 `_meta_*` 私有字段与属性映射。
- HDF5 写入/读取时严格遵循组/数据集命名与类型约定，避免后续解析不一致。

---
如需将该规范扩展到测试、日志、异常分类等，可在本文件增补相应章节并给出代码片段示例。