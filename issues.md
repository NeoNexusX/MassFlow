# MassFlow 代码审阅报告

- **审阅对象**：`dev` 分支 `a70b5b3`，另含未纳入版本控制的 `scripts/imzml_to_zarr.py`
- **审阅日期**：2026-09-26
- **范围**：`src/massflow` 全部模块（约 2 万行）、`tests/`、中英文 `docs/`、`main.py`、`scripts/`、`.github/workflows/`、`pyproject.toml`
- **方法**：逐文件阅读代码；可疑点写独立的小脚本实测复现（关键复现代码见[附录](#8-附录关键问题复现代码)）；文档中的函数签名默认值用脚本与 `inspect.signature` 自动比对。

**标记说明**

| 标记 | 含义 |
|---|---|
| P0 | 正常使用即可触发：数据被写错、流程崩溃，或污染用户环境 |
| P1 | 明显的逻辑错误或显著的性能问题 |
| P2 | 边界条件、API 一致性问题、潜在隐患 |
| P3 | 小问题、代码质量 |
| 实测 | 已用脚本复现 |
| 读码 | 根据代码推断，未单独复现 |

**问题统计**

| 分类 | 数量 | P0 | P1 | P2 | P3 |
|---|---|---|---|---|---|
| 1. 严重问题 | 12 | 7 | 5 | 0 | 0 |
| 2. 逻辑错误与隐患 | 28 | 0 | 3 | 15 | 10 |
| 3. 性能 | 10 | 0 | 2 | 7 | 1 |
| 4. 文档与代码不一致 | 11 | 0 | 3 | 4 | 4 |
| 5. 测试、CI 与工程化 | 3 | 0 | 2 | 1 | 0 |
| 6. 代码质量与小问题 | 10 | 0 | 0 | 0 | 10 |
| **合计** | **74** | **7** | **15** | **27** | **25** |

**测试现状**（本地运行 `pytest`）：97 通过、21 失败、14 错误、64 跳过。35 个失败/错误中有 34 个、以及全部 64 个跳过，都是因为测试数据不在仓库里。唯一真实的失败是 `test_mean_spectrum_is_absent_by_default`，而且需要改的是测试本身（见 TEST-03）。

**建议优先处理**：BUG-01～BUG-07。它们在正常使用中就会触发，而且改动都不大。

## 目录

1. [严重问题](#1-严重问题)
2. [逻辑错误与隐患](#2-逻辑错误与隐患)
3. [性能](#3-性能)
4. [文档与代码不一致](#4-文档与代码不一致)
5. [测试、CI 与工程化](#5-测试ci-与工程化)
6. [代码质量与小问题](#6-代码质量与小问题)
7. [实测确认没有问题的部分](#7-实测确认没有问题的部分)
8. [附录：关键问题复现代码](#8-附录关键问题复现代码)

---

## 1. 严重问题

- [ ] **BUG-01 [P0][实测] imzML 写出器在 auto 模式下会把不同谱的 m/z 写成同一条**
  - 位置：[stream_imzml_writer.py:378-393](src/massflow/tools/stream_imzml_writer.py#L378-L393)
  - 问题：缓存用 `(长度, 首值, 尾值)` 判断两条 m/z 是否相同。两条 processed 谱只要这三项一样，第二条就会复用第一条的 m/z 偏移，而且没有任何报错。元数据没有标明 continuous 或 processed 时，`MSDataManagerImzML.add_spectrum` 和 `swap_matrix_data_out2disk` 都会走这个模式。
  - 复现：依次写入 `[100, 150.2, 200]` 和 `[100, 170.5, 200]`，读回的第二条是 `[100, 150.2, 200]`。
  - 建议：用完整数组比较（或对全部内容做哈希）作为缓存键。

- [x] **BUG-02 [P0][实测] `fastmath=True` 让 NaN/Inf 检查失效**
  - 位置：[baseline_correction_numba.py:6](src/massflow/preprocess/numba/baseline_correction_numba.py#L6)、[:32](src/massflow/preprocess/numba/baseline_correction_numba.py#L32)、[:61](src/massflow/preprocess/numba/baseline_correction_numba.py#L61)；[peak_align_numba.py:6](src/massflow/preprocess/numba/peak_align_numba.py#L6)、[:176](src/massflow/preprocess/numba/peak_align_numba.py#L176)、[:231](src/massflow/preprocess/numba/peak_align_numba.py#L231)
  - 问题：fastmath 允许 LLVM 假定数据里没有 NaN/Inf，`np.isfinite`/`np.isnan` 可能被优化掉，和 NaN 的比较结果也变成未定义。
  - 复现：`_median_of_finite([1, nan, 3, 5])` 返回 `4.0`，正确值是 `3.0`；`estimate_domain_stats_flat_parallel_jit` 遇到 NaN 时 max 变成 `nan`。此外，`search_nearest_jit` 用 `np.nan` 做比较哨兵（[:213](src/massflow/preprocess/numba/peak_align_numba.py#L213)），`scalar_diff_jit` 在分母为 0 时返回 NaN，两者也都依赖 NaN 的语义。
  - 建议：依赖 NaN 语义的函数去掉 fastmath，或者改用不含 `nnan`/`ninf` 的细粒度 fastmath 标志。

- [x] **BUG-03 [P0][实测] 连续调用 `Preprocessor` 且不传 `output_path` 时，第二步会读坏自己的输入**
  - 位置：[planner.py:42-43](src/massflow/preprocess/pipeline/planner.py#L42-L43)
  - 问题：默认输出路径写死为 `temp_dir/preprocess_output.<后缀>`。第二次调用时，输入和输出是同一个文件，写出器以 `"w"` 打开时会先把它截断。代码中也没有检查"输出路径等于输入路径"的情况。Zarr 输出同理（读码）。
  - 复现：先 `step1 = Preprocessor(dm).filter_peaks_by_isotopes().start()`，再 `Preprocessor(step1).peak_align().start()`，报 `could not broadcast input array from shape (0,) into shape (40,)`。
  - 建议：默认文件名加上 uuid；路径解析后如果与输入相同，直接报错。

- [ ] **BUG-04 [P0][实测] H5 离子图在 finalize 时崩溃**
  - 位置：[ion_image_data_manager_h5.py:206-211](src/massflow/data_manager/ion_image_data_manager_h5.py#L206-L211)、[:428](src/massflow/data_manager/ion_image_data_manager_h5.py#L428)
  - 问题：`dense[:, flat_idx] = sparse` 默认所有像素都有效，只要有一个坐标越界就报 shape mismatch。另外，判断坐标从 0 还是从 1 开始时只看 x 的最小值（第 208 行）：如果 y 含 0 而 x 不含 0，程序会把 x、y 整体减 1，y 就变成 -1，同样崩溃。同一组坐标在 Binary 后端下输出正确。
  - 建议：散射时只取有效列（`sparse[:, valid]`），并和其他后端统一坐标基准（见 LOGIC-06）。

- [ ] **BUG-05 [P0][实测] 离子图写出器为每个临时 batch 保留一个打开的文件描述符**
  - 位置：[ion_major_temp_store.py:128-132](src/massflow/tools/ion_major_temp_store.py#L128-L132)
  - 问题：`read_ion_block` 按 batch 缓存 fd，直到 `close()` 才关闭，batch 数一旦超过进程的打开文件上限就会失败。Binary、H5、Zarr 三种 ion-image 输出都走这条路径。
  - 复现：把软上限设为 256（macOS 上 launchd 的默认值，从 Finder/Dock 启动的 Jupyter 会继承），写入 300 个 batch 后调用 `finalize()`，报 `OSError: [Errno 24] Too many open files`。按 `batch_size=256` 计算，约 6.4 万像素就会触发；Linux 常见的默认上限 1024 对应约 26 万像素。
  - 建议：每读完一个 band 就关闭 fd，或者只用一个容量很小的 fd LRU，或者把临时 batch 合并写入同一个文件。

- [ ] **BUG-06 [P0][实测] import 时修改宿主程序的全局日志和警告配置**
  - 位置：[logger.py:273](src/massflow/tools/logger.py#L273)、[:118-119](src/massflow/tools/logger.py#L118-L119)、[:136-138](src/massflow/tools/logger.py#L136-L138)、[:150-151](src/massflow/tools/logger.py#L150-L151)
  - 问题：import 任意一个 massflow 模块都会产生以下副作用：
    - 在当前工作目录创建 `logs/`；
    - 删除 root logger 上已有的 handler，而且清理循环边遍历边删除，实测会漏删；
    - 开启 `captureWarnings`，并替换全局的 `warnings.formatwarning`；
    - `massflow` logger 设为 `propagate=False`，用户无法通过 root logger 收到它的日志。
  - 建议：库内只使用 `logging.getLogger(__name__)` 并挂一个 `NullHandler`；`init_app()` 改为由用户显式调用。

- [ ] **BUG-07 [P0][实测] 峰提取后只要有像素没有峰，输出为 `spectra_zarr` 就会崩溃**
  - 位置：[msi_zarr/writer.py:159-160](src/massflow/msi_zarr/writer.py#L159-L160)
  - 问题：Zarr 写出器拒绝长度为 0 的谱。MSI 数据中组织外的全零背景像素很常见，这些像素在峰提取后长度为 0，于是 `Preprocessor(..., output_format="spectra_zarr").peak_pick().start()` 直接失败。同样的数据输出为 imzML 则正常。用 `scripts/imzml_to_zarr.py` 转换含空谱的 processed imzML 时，也会在这一步失败（读码）。
  - 复现：见附录 B7，报 `ValueError: Writer does not accept empty spectra.`。
  - 建议：允许空谱（offsets 可以表示长度为 0 的行），同时调整 auto 编码推断中对长度为 0 的处理。

- [ ] **BUG-08 [P1][实测] SNIP 基线的结果随 batch 大小变化**
  - 位置：[baseline_correction_numba.py:300-304](src/massflow/preprocess/numba/baseline_correction_numba.py#L300-L304)
  - 问题：默认窗口 `m` 用整个 batch 的总点数计算，而不是单条谱的长度。
  - 复现：同一条 300 点的谱，batch=1 和 batch=64 算出的基线最大相差 0.17。
  - 建议：在 jit 函数内按每条谱的长度计算默认 `m`。

- [ ] **BUG-09 [P1][实测/读码] `ref` 归一化会把整条谱清零，默认参考 m/z 还随 batch 变化**
  - 位置：[normalization_numba.py:84-85](src/massflow/preprocess/numba/normalization_numba.py#L84-L85)、[:99](src/massflow/preprocess/numba/normalization_numba.py#L99)、[:127](src/massflow/preprocess/numba/normalization_numba.py#L127)；[flat_pre_fun.py:255-256](src/massflow/preprocess/flat_pre_fun.py#L255-L256)
  - 问题：
    - 找不到参考峰时函数直接 return，而输出数组是用 `np.zeros` 初始化的，这条谱就变成全 0（实测）。
    - 没有传 `ref` 时，程序取当前 flat m/z 数组中间的那个值。对 processed 数据，这个值每个 batch 都不一样，各 batch 的归一化基准不一致（读码）。
  - 建议：找不到参考峰时保留原值或写入 NaN 并告警；`ref` 为 None 时直接报错，或者在整个数据集上只确定一次。

- [ ] **BUG-10 [P1][实测/读码] `start()` 会改写已注册的任务**
  - 位置：[runner.py:75-80](src/massflow/preprocess/pipeline/runner.py#L75-L80)、[planner.py:137-138](src/massflow/preprocess/pipeline/planner.py#L137-L138)
  - 问题：runner 把算出来的 reference/tolerance 写回了 `Preprocessor._tasks` 里的同一个任务对象（实测），在 Jupyter 里第二次调用 `start()` 时会沿用上一次的参考轴。planner 还会覆盖 dataset 任务的 `output_path`，用户显式传入的值也会被盖掉（读码）。
  - 建议：规划前对任务做 `copy.deepcopy`，运行过程中产生的结果只写进副本。

- [ ] **BUG-11 [P1][实测] 文件名带点时路径解析出错，四处规则也不一致**
  - 位置：[ms_data_manager.py:70-73](src/massflow/data_manager/ms_data_manager.py#L70-L73)；对比 [stream_imzml_writer.py:204](src/massflow/tools/stream_imzml_writer.py#L204)、[export_imzml.py:361-364](src/massflow/tools/export_imzml.py#L361-L364)、[scripts/imzml_to_zarr.py:113](scripts/imzml_to_zarr.py#L113)
  - 问题：DataManager 会循环剥掉所有扩展名，`sample.v2.imzML` 被解析成 `sample.imzML`，结果找不到文件；另外三处都只剥掉一层。
  - 建议：只去掉已知后缀（`.imzML`、`.ibd`、`.zarr`，不区分大小写），四处共用同一个函数。

- [ ] **BUG-12 [P1][实测/读码] 生成器会永久改掉 DataManager 的线程数**
  - 位置：[ms_data_manager_imzml.py:270-276](src/massflow/data_manager/ms_data_manager_imzml.py#L270-L276)
  - 问题：调用 `flat_generator(max_threads=2)` 之后，`dm.max_threads` 从 8 变成了 2（实测）。画图函数和 `reference_computer` 默认传 2，dataset executor 传 `output_max_workers`，所以用户的输入 DataManager 会被悄悄降速。如果 `output_max_workers=1`，还会直接报 "Single-thread mode is not supported"（读码）。
  - 建议：生成器内部使用局部的线程数，不要写回实例属性。

## 2. 逻辑错误与隐患

- [ ] **LOGIC-01 [P1][读码] Zarr 转 imzML 时元数据被丢弃，MassFlow 写出的 imzML 都没有像素尺寸**
  - 位置：[ms_data_manager_imzml.py:260-265](src/massflow/data_manager/ms_data_manager_imzml.py#L260-L265)、[stream_imzml_writer.py:648-688](src/massflow/tools/stream_imzml_writer.py#L648-L688)
  - 问题：`copy_meta` 要求源元数据和目标是同一类型，从 `MetaDataZarr` 复制到 `MetaDataImzMl` 时只打一条警告、什么都不复制；而反方向的 `MetaDataZarr.copy` 支持跨类型，两边不对称。另外，`ImzMLWriter.from_ms_meta` 本身就不写像素尺寸（IMS:1000046/47）和极性，所以经 MassFlow 写出的 imzML 都会丢掉像素尺寸。
  - 建议：按字段名跨类型复制；模板补上 pixel size 和 polarity。

- [ ] **LOGIC-02 [P1][读码] 谱类型未知时一律写成 centroid，重新读入后 profile 步骤会被丢弃**
  - 位置：[stream_imzml_writer.py:659](src/massflow/tools/stream_imzml_writer.py#L659)
  - 问题：只要 `profile_spectrum` 不是 True，就写 centroid。之后读入这个文件时，planner 会把它当 centroid 数据，丢掉 baseline、noise 和 peak_pick 任务。
  - 建议：类型未知时不写谱类型的 cvParam，或者在写出前先推断。

- [ ] **LOGIC-03 [P2][读码] 写出的 imzML 不完全符合规范**
  - 位置：[stream_imzml_writer.py:30-33](src/massflow/tools/stream_imzml_writer.py#L30-L33)、[:68](src/massflow/tools/stream_imzml_writer.py#L68)、[:528](src/massflow/tools/stream_imzml_writer.py#L528)、[:541-552](src/massflow/tools/stream_imzml_writer.py#L541-L552)、[:706-711](src/massflow/tools/stream_imzml_writer.py#L706-L711)
  - 问题：
    - `cvList` 没有声明 `IMS`，文件里却大量使用 `cvRef="IMS"`；
    - `ms level` 写成了 0；
    - flat 写出路径忽略压缩配置，XML 里却照样声明压缩方式（默认不压缩，所以目前不会触发）；
    - `__exit__` 遇到异常时会留下半成品文件。

- [ ] **LOGIC-04 [P2][实测] 连续读取路径返回只读数组**
  - 位置：[read_matrix.py:148-152](src/massflow/tools/read_matrix.py#L148-L152)、[read_flat.py:37-41](src/massflow/tools/read_flat.py#L37-L41)
  - 问题：磁盘 dtype 与输出 dtype 相同时，程序直接返回 `frombuffer` 的只读视图，而碎片化读取路径返回的是可写数组。用户对结果做原地运算时会报错，而且报不报错取决于数据在磁盘上是否连续。
  - 建议：统一返回可写数组，例如用 `readinto` 读进预分配的数组。

- [ ] **LOGIC-05 [P1][实测] H5 数据集无法 pickle，在 macOS 上用多进程 DataLoader 会失败**
  - 位置：[ion_image_h5_dataset.py:85-89](src/massflow/morch/ion_image_h5_dataset.py#L85-L89)
  - 问题：数据集被访问过一次之后会持有 `h5py.File`，此时 `pickle.dumps` 报 `h5py objects cannot be pickled`。macOS 默认用 spawn 启动 worker，需要 pickle 数据集。Binary 版本没有这个问题。
  - 建议：实现 `__getstate__`，把 `_file` 和 `_dset` 丢掉。

- [ ] **LOGIC-06 [P2][读码] 判断坐标从 0 还是从 1 开始，代码里有 6 套不同的规则**
  - 位置（括号内是各处的判断方式）：
    - [ion_image_data_manager_binary.py:608](src/massflow/data_manager/ion_image_data_manager_binary.py#L608)（x 或 y 含 0）
    - [ion_image_data_manager_h5.py:208](src/massflow/data_manager/ion_image_data_manager_h5.py#L208)（只看 x 的最小值）
    - [morch/ion_image_binary_dataset.py:63](src/massflow/morch/ion_image_binary_dataset.py#L63)（只看 x）
    - [plot_ion_image.py:175-177](src/massflow/tools/plot_ion_image.py#L175-L177) 与 [:480](src/massflow/tools/plot_ion_image.py#L480)
    - [region_mask.py:78](src/massflow/segmentation/region_mask.py#L78)
    - [umap_kmeans.py:692-695](src/massflow/segmentation/umap_kmeans.py#L692-L695)
  - 问题：同一份数据在不同模块里可能被当成不同的基准。[region_mask.py:66-79](src/massflow/segmentation/region_mask.py#L66-L79) 没有越界检查，坐标为 -1 时会静默取到最后一行或最后一列；[:165](src/massflow/segmentation/region_mask.py#L165) 读取的 `_coordinate_base` 属性在任何类里都不存在。
  - 建议：imzML 规范和 Zarr schema 都规定坐标从 1 开始。可以统一成一个函数，或者直接固定为从 1 开始。

- [ ] **LOGIC-07 [P2][读码] processed 数据的离子图取"最近峰"时不设容差**
  - 位置：[plot_ion_image.py:146-150](src/massflow/tools/plot_ion_image.py#L146-L150)
  - 问题：`mz_tolerance=None` 时，每个像素都取离目标最近的峰，即使相距几十 Da 也会被画进离子图。
  - 建议：没给容差时使用一个默认的 ppm 容差，超出范围的像素填 0。

- [ ] **LOGIC-08 [P2][读码] Cardinal 后端的若干问题**
  - 位置：[api.py:261](src/massflow/preprocess/api.py#L261)、[:266](src/massflow/preprocess/api.py#L266)；[adapter.py:260-261](src/massflow/r_preprocess/adapter.py#L260-L261)
  - 问题：
    - `peak_align` 的 `min_frequency` 默认值是 0.01，所以"只有 Python 后端支持"的警告每次都会打印；
    - `cardinal_runner` 直接读 `data_manager.imzml_filepath`，输入是 Zarr 时会报 AttributeError；
    - R 端会读取整个文件，`target_locs` 被静默忽略。

- [ ] **LOGIC-09 [P2][读码] R 环境初始化的问题**
  - 位置：[environment.py:78-82](src/massflow/r_preprocess/environment.py#L78-L82)、[:109-111](src/massflow/r_preprocess/environment.py#L109-L111)、[:42-43](src/massflow/r_preprocess/environment.py#L42-L43)；[r_preprocess/\_\_init\_\_.py:5-12](src/massflow/r_preprocess/__init__.py#L5-L12)
  - 问题：
    - 单例先给 `_instance` 赋值再初始化，初始化失败后再调用，拿到的是一个半初始化的对象；
    - rpy2 没有安装时提示 "This error should not occur"，而没有提示安装 `--extra r`；
    - `LC_ALL=C` 改的是整个 Python 进程的环境变量；
    - `set_default_r_home` 先赋值、后校验路径。

- [ ] **LOGIC-10 [P2][读码] 没有注册任务且输出为 imzML 时，`output_path` 被静默忽略**
  - 位置：[preprocessor.py:140-149](src/massflow/preprocess/preprocessor.py#L140-L149)
  - 问题：程序直接返回输入的 DataManager，不写任何文件，也不给提示。

- [ ] **LOGIC-11 [P2][读码] 各模块支持的数据类型不一致**
  - 位置：[plot_ion_image.py:28-31](src/massflow/tools/plot_ion_image.py#L28-L31)、[region_mask.py:32-37](src/massflow/segmentation/region_mask.py#L32-L37)
  - 问题：流水线可以输出 `spectra_zarr`，UMAP/KMeans 也接受 `MSDataManagerZarr`，但 `read_ion_image`、`plot_ion_image`、`plot_tic_image` 和 `apply_region_mask` 都不支持它。结果是 `SegmentationPipeline` 用 `MSDataManagerZarr` 能跑完 UMAP，到 `apply_region_mask` 这一步才报错。

- [ ] **LOGIC-12 [P2][读码] 输入数据可能被意外修改或删除**
  - 位置：[umap_kmeans.py:938-943](src/massflow/segmentation/umap_kmeans.py#L938-L943)、[region_mask.py:203](src/massflow/segmentation/region_mask.py#L203)、[ms_data_manager.py:247-248](src/massflow/data_manager/ms_data_manager.py#L247-L248)、[ms_data_manager_imzml.py:254-258](src/massflow/data_manager/ms_data_manager_imzml.py#L254-L258)
  - 问题：
    - UMAP 的 `save_matrix="zarr"` 以 `r+` 模式把结果写进源 MSI Zarr；
    - Binary 后端的 region mask 在输出路径等于输入路径时，会先截断源文件再读取；
    - `filepath=None` 的 imzML DataManager 被垃圾回收时（`__del__`）会自动删除文件，只保存了路径字符串的调用方会发现文件不见了。

- [ ] **LOGIC-13 [P2][读码] 构造 DataManager 时会覆盖传入的 `ms` 的元数据**
  - 位置：[ms_data_manager_imzml.py:76](src/massflow/data_manager/ms_data_manager_imzml.py#L76)；Zarr 版本在 [ms_data_manager_zarr.py:65](src/massflow/data_manager/ms_data_manager_zarr.py#L65)，问题相同。

- [ ] **LOGIC-14 [P2][读码] `units="relative"` 时 tolerance 的含义不明确**
  - 位置：[api.py:290](src/massflow/preprocess/api.py#L290)、[peak_align_helper.py:41-43](src/massflow/preprocess/helper/peak_align_helper.py#L41-L43)、[:449](src/massflow/preprocess/helper/peak_align_helper.py#L449)
  - 问题：只有 `units == "ppm"` 时才会乘以 1e-6，而 `_normalize_units` 把 `"relative"` 和 `"ppm"` 当成同一类。所以传入 `units="relative", tolerance=10` 时，10 会被当成相对值 10（即 1000%）。

- [ ] **LOGIC-15 [P2][读码] ion-image 的 `from_data_manager` 没有处理元数据中宽高为 0 的情况**
  - 位置：[ion_image_data_manager_zarr.py:246-247](src/massflow/data_manager/ion_image_data_manager_zarr.py#L246-L247)、[ion_image_data_manager_binary.py:236-237](src/massflow/data_manager/ion_image_data_manager_binary.py#L236-L237)、[ion_image_data_manager_h5.py:337-338](src/massflow/data_manager/ion_image_data_manager_h5.py#L337-L338)
  - 问题：宽高直接取 `int(meta.max_count_of_pixels_x/y)`，而这两个字段的默认值是 0。宽高为 0 时，Zarr 会把所有坐标判为越界，Binary 会写出宽度为 0 的文件（离子图全空），H5 则会崩溃。`MSDataManagerZarr` 对 `<= 0` 做了保护，这三个后端没有。

- [ ] **LOGIC-16 [P2][实测/读码] `MassSpectrumSet`、`PixelCoordinates`、`Spectrum` 的 API 缺陷**
  - 位置：[mass_spectrum_set.py:229-236](src/massflow/module/mass_spectrum_set.py#L229-L236)、[:271-277](src/massflow/module/mass_spectrum_set.py#L271-L277)；[pixel_coordinates.py:39-58](src/massflow/module/pixel_coordinates.py#L39-L58)、[:90](src/massflow/module/pixel_coordinates.py#L90)、[:92-111](src/massflow/module/pixel_coordinates.py#L92-L111)；[spectrum.py:150-155](src/massflow/module/spectrum.py#L150-L155)
  - 问题：
    - `ms[x, y, z, w]` 静默返回 None（实测）；
    - `ms[x, y]` 读取时优先取 z=0，写入时却取第一个 z，两者不一致；
    - numpy 整数在构造函数里被拒绝，在 y/z 的 setter 里被静默忽略（实测）；
    - y 和 z 的报错信息都写成了 "X coordinate"；
    - `Spectrum` 定义了 `__eq__` 却没有定义 `__hash__`，因此不可哈希（实测）。

- [ ] **LOGIC-17 [P2][读码] DataManager 基类的边界处理**
  - 位置：[ms_data_manager.py:79-81](src/massflow/data_manager/ms_data_manager.py#L79-L81)、[:182-193](src/massflow/data_manager/ms_data_manager.py#L182-L193)、[:202](src/massflow/data_manager/ms_data_manager.py#L202)、[:214-215](src/massflow/data_manager/ms_data_manager.py#L214-L215)
  - 问题：
    - `_init_target_locs` 抛出的 `ValueError` 不带任何信息；
    - `ms` 为空时（例如 `target_locs` 没有选中任何像素），`detect_sorted` 调用 `randint(0, -1)` 会崩溃；而且它只随机抽一条谱来判断是否有序，结果不确定；
    - `create_ms_meta_mask` 的 `is None` 检查永远不成立（字段默认值是 0）；
    - 坐标用 int16 存储，超过 32767 会溢出。

- [ ] **LOGIC-18 [P2][读码] 同一参数在各层的默认值不一致**

  | 参数 | helper / 底层函数 | API 或另一处 | 位置 |
  |---|---|---|---|
  | isotope `max_charge` | 1 | 3 | [isotope_filter_helper.py:59](src/massflow/preprocess/helper/isotope_filter_helper.py#L59) |
  | smoother `polyorder` | 2 | 3 | [noise_reduction_helper.py:41](src/massflow/preprocess/helper/noise_reduction_helper.py#L41) |
  | `noise_reduction_flat` 的 `method` | `ma_numba` | `gaussian_numba` | [flat_pre_fun.py:79](src/massflow/preprocess/flat_pre_fun.py#L79) |
  | locmin `niter` | 3 | 15 | [baseline_correction_numba.py:319](src/massflow/preprocess/numba/baseline_correction_numba.py#L319) |
  | `pipeline_stall_timeout` | 60（`InputMetaData`） | 180（`Preprocessor`） | [types.py:62](src/massflow/preprocess/pipeline/types.py#L62)、[preprocessor.py:59](src/massflow/preprocess/preprocessor.py#L59) |

  - 影响：通过 `Preprocessor` API 调用时不受影响；直接调用 helper 或自己构造 `InputMetaData` 时，行为会和 API 不同。

- [ ] **LOGIC-19 [P3][读码] 任务顺序被静默重排，重复注册的任务会执行两次**
  - 位置：[planner.py:20-30](src/massflow/preprocess/pipeline/planner.py#L20-L30)、[:58-77](src/massflow/preprocess/pipeline/planner.py#L58-L77)
  - 问题：planner 按固定顺序重排任务，用户注册的顺序被忽略，只在 info 日志里能看出来；同一个任务注册两次会执行两次。另外，`Preprocessor` 的类型注解只写了 `MSDataManagerImzML`（[preprocessor.py:44](src/massflow/preprocess/preprocessor.py#L44)），实际上 Zarr 输入也支持。

- [ ] **LOGIC-20 [P3][读码] 流水线的停滞检测可能误杀正常运行的任务**
  - 位置：[flat_executor.py:446-464](src/massflow/preprocess/pipeline/flat_executor.py#L446-L464)
  - 问题：只有在队列存取时才记录一次进度。如果单个 batch 的处理时间超过 `pipeline_stall_timeout`，就会被判定为停滞并中止。例如缓存为空时首次 JIT 编译全部内核，或者参考轴极长时做对齐，都可能超时。

- [ ] **LOGIC-21 [P3][读码] 噪声估计与 Cardinal/matter 的细节差异**
  - 位置：[est_noise_numba.py:381-395](src/massflow/preprocess/numba/est_noise_numba.py#L381-L395)、[:505-544](src/massflow/preprocess/numba/est_noise_numba.py#L505-L544)
  - 问题：`overlap == 0` 时，`dynamic=True` 会被静默忽略；`nbins > 1` 时，matter 用 `spline()` 做三次样条插值，这里用的是局部三次 Lagrange 插值，结果会不同。

- [ ] **LOGIC-22 [P3][读码] 局部极大值检测的边界规则在两处不一致**
  - 位置：[peak_pick_numba.py:42](src/massflow/preprocess/numba/peak_pick_numba.py#L42)（`i > n - r`）；[baseline_correction_numba.py:229](src/massflow/preprocess/numba/baseline_correction_numba.py#L229)（`i >= n - r`）
  - 问题：峰提取在 `i = n - r` 处会接受一个窗口被截断的"峰"。

- [ ] **LOGIC-23 [P3][读码] 写入和读取两端的 API 不一致**
  - 位置：[msi_zarr/writer.py:126-131](src/massflow/msi_zarr/writer.py#L126-L131)、[ion_image_data_manager.py:51-52](src/massflow/data_manager/ion_image_data_manager.py#L51-L52)、[ion_image_data_manager_zarr.py:46](src/massflow/data_manager/ion_image_data_manager_zarr.py#L46)、[ion_image_data_manager_h5.py:66](src/massflow/data_manager/ion_image_data_manager_h5.py#L66)
  - 问题：
    - `MSIZarrWriter.add_spectrum(intensity, coordinates, mz)` 和 `ImzMLWriter.add_spectrum(mz, intensity, coords)` 的参数顺序不同；
    - 区间为空时，基类的 `get_ion_block` 返回形状 `(0, 0)`，H5 返回 `(0, n_pixels)`；
    - Zarr 后端直接 `del merge_workers`，H5 的 `_merge_workers` 从未被使用。

- [ ] **LOGIC-24 [P3][读码] 资源清理不完整**
  - 位置：[ms_data_manager_zarr.py:268-276](src/massflow/data_manager/ms_data_manager_zarr.py#L268-L276)、[runner.py:52-86](src/massflow/preprocess/pipeline/runner.py#L52-L86)
  - 问题：`MSDataManagerZarr.close()` 在还没有 finalize 时直接关闭，既不报错也不警告，留下一个 `write_state="writing"`、无法读取的 store。`run_segments` 中途出错时不关闭中间的 DataManager，临时文件要等垃圾回收时才会被删除。

- [ ] **LOGIC-25 [P3][读码] 潜在缺陷：目前没有调用方，将来启用时会出问题**
  - 位置：[spectrum.py:97-116](src/massflow/module/spectrum.py#L97-L116)、[:214-220](src/massflow/module/spectrum.py#L214-L220)；[spectrum_imzml.py:101-114](src/massflow/module/spectrum_imzml.py#L101-L114)；[ms_meta_data.py:191-194](src/massflow/module/ms_meta_data.py#L191-L194)
  - 问题：
    - 使用 `shared_mz_list` 时，`sort_by_mz = True` 只排序了 `_mz_list`，getter 却返回共享的 m/z 轴，于是 m/z 与强度错位；
    - `SpectrumImzML.swap_out2disk` 在持有 reader 时不写出修改，只清空数据，之后会从原文件重新读回原始数据；
    - `crop_range(mode="new")` 生成的副本保留了 `_reader`，调用 `clear_data()` 后会重新加载未裁剪的数据；
    - `update_meta()` 会把通过 `_set` 写入的自定义键全部改成 None。

- [ ] **LOGIC-26 [P3][读码] 谱图对比指标框的计算方式不对**
  - 位置：[plot.py:172-201](src/massflow/tools/plot.py#L172-L201)
  - 问题：两条谱的强度是按下标逐点比较的。例如对比 profile 谱和峰提取后的谱时，两者的 m/z 轴完全不同，这样算出的相关性没有意义。此外，指标框忽略 `mz_range`，显示的是整条谱的范围。

- [ ] **LOGIC-27 [P3][读码] 默认把文件写进当前工作目录**
  - 位置：`logs/`（[logger.py:101](src/massflow/tools/logger.py#L101)）、`./temp`（[preprocessor.py:51](src/massflow/preprocess/preprocessor.py#L51)）、`output/<uuid>`（[ion_image_data_manager_binary.py:163](src/massflow/data_manager/ion_image_data_manager_binary.py#L163)、[ion_image_data_manager_h5.py:57](src/massflow/data_manager/ion_image_data_manager_h5.py#L57)、[ion_image_data_manager_zarr.py:52](src/massflow/data_manager/ion_image_data_manager_zarr.py#L52)）、`./tic_image_output`（[plot_ion_image.py:728](src/massflow/tools/plot_ion_image.py#L728)）、`./segmentation_output`、`./umap_kmeans_output`
  - 建议：库函数不应默认写入当前工作目录。确实需要默认值时，放到 `tempfile.gettempdir()` 下，或者要求调用方显式传入路径。

- [ ] **LOGIC-28 [P3][读码] processed 输入默认不检查 m/z 是否升序**
  - 位置：[preprocessor.py:60](src/massflow/preprocess/preprocessor.py#L60)、[validation.py:81-82](src/massflow/preprocess/pipeline/validation.py#L81-L82)
  - 问题：`validate_processed_mz` 默认是 False。而对齐时的 `searchsorted`、参考轴构建和峰面积积分都假定 m/z 是升序的，m/z 未排序的 processed 数据会静默地得出错误结果。
  - 建议：默认抽样检查，或者在 `load_head_data()` 时做一次轻量检查。

## 3. 性能

- [ ] **PERF-01 [P1][读码] UMAP 在 Zarr/H5 离子图上的读取量是平方级的**
  - 位置：[umap_kmeans.py:459-510](src/massflow/segmentation/umap_kmeans.py#L459-L510)、[:514-588](src/massflow/segmentation/umap_kmeans.py#L514-L588)
  - 问题：transform 阶段每处理一个像素 batch，都会把全部 ion block（也就是全部数据）读一遍，再从中切出需要的列。10 万像素、batch 为 1024 时，大约要把全量数据读 100 次；H5 每读一次还要解压整个数据集。
  - 建议：改成外层按 ion block 循环，一次性为所有像素累积；或者先把数据转置成 pixel-major 的临时文件。

- [ ] **PERF-02 [P1][读码] 峰提取对每条谱做两遍完整的检测**
  - 位置：[peak_pick_numba.py:346-365](src/massflow/preprocess/numba/peak_pick_numba.py#L346-L365)、[:373-404](src/massflow/preprocess/numba/peak_pick_numba.py#L373-L404)
  - 问题：第一遍用来计数，第二遍用来写出，局部极大值检测、prominence 过滤和 SNR 过滤都完整做了两次。
  - 建议：第一遍把峰的索引存进按上界预分配的缓冲区，第二遍只做拷贝。

- [ ] **PERF-03 [P2][读码] 统计聚合的内存访问是跨步的**
  - 位置：[stat_filter_numba.py:24-34](src/massflow/preprocess/numba/stat_filter_numba.py#L24-L34)、[:49-57](src/massflow/preprocess/numba/stat_filter_numba.py#L49-L57)、[:76-81](src/massflow/preprocess/numba/stat_filter_numba.py#L76-L81)
  - 问题：外层按峰循环（`prange`）、内层按谱循环，每次访问的步长都是一整条谱的长度，对缓存很不友好。
  - 建议：外层按谱分块，每个线程使用局部累加器，最后再归约。

- [ ] **PERF-04 [P2][读码] Zarr reader 存在重复读取**
  - 位置：[msi_zarr/reader.py:173-174](src/massflow/msi_zarr/reader.py#L173-L174)、[:217-219](src/massflow/msi_zarr/reader.py#L217-L219)、[:359-366](src/massflow/msi_zarr/reader.py#L359-L366)
  - 问题：
    - `read_spectrum` 每读一条谱都会把共享 m/z 轴整个读一遍，而 `SpectrumZarr` 懒加载时每条谱都会触发它；
    - `read_flat_batch` 每个 batch 都重读一次 m/z 轴；
    - 索引不连续时逐条谱读取 zarr，同一个 chunk 会被反复解压。
  - 建议：缓存共享 m/z 轴；对不连续的索引按 chunk 合并读取。

- [ ] **PERF-05 [P2][读码] 给了 reference 但没给 tolerance 时，参考轴会被白白重建一遍**
  - 位置：[api.py:280-288](src/massflow/preprocess/api.py#L280-L288)、[peak_align_helper.py:451-456](src/massflow/preprocess/helper/peak_align_helper.py#L451-L456)
  - 问题：`_reference_computer` 没有收到用户传入的 reference，于是完整地重新分箱一遍（两次全量读取），最后却只用到了 tolerance。
  - 建议：把 reference 传进去，只估计 tolerance。

- [ ] **PERF-06 [P2][读码] 对齐内核在内循环里分配数组**
  - 位置：[peak_align_numba.py:77-78](src/massflow/preprocess/numba/peak_align_numba.py#L77-L78)、[:89-90](src/massflow/preprocess/numba/peak_align_numba.py#L89-L90)
  - 问题：在上采样分支中，每处理一个参考点，都要新建一个 `np.array([ref])` 和一个 `search_nearest_jit` 的结果数组。
  - 建议：改用标量版本的最近邻查找。

- [ ] **PERF-07 [P2][读码] Zarr 写出器的 compaction 与切分效率低**
  - 位置：[msi_zarr/writer.py:387-405](src/massflow/msi_zarr/writer.py#L387-L405)、[:529-551](src/massflow/msi_zarr/writer.py#L529-L551)、[:988-1000](src/massflow/msi_zarr/writer.py#L988-L1000)、[:439-443](src/massflow/msi_zarr/writer.py#L439-L443)
  - 问题：
    - 坐标去重之后的 compaction 逐行读写 zarr，同一个 chunk 会被反复解压、重写；
    - `_split_source_indices` 对每个 batch 都扫描一遍整个 `source_idx`，复杂度是 O(batch 数 × 像素数)，而 `source_idx` 本来就有序，可以用 `searchsorted`；
    - offsets 的 chunk 大小取第一个 batch 的长度加 1，大数据集会因此生成大量小文件。

- [ ] **PERF-08 [P2][实测] 重复调用 `load_head_data()` 的耗时是 O(n²)，计数也会翻倍**
  - 位置：[mass_spectrum_set.py:128-139](src/massflow/module/mass_spectrum_set.py#L128-L139)、[ms_data_manager_imzml.py:218-219](src/massflow/data_manager/ms_data_manager_imzml.py#L218-L219)
  - 问题：坐标已经存在时，程序用 `_queue.index(existing)` 线性查找旧谱（每次比较都要调用 Python 层的 `__eq__`）。实测 3000 条谱，第一次加载 0.09 秒，第二次 0.38 秒；按平方增长推算，10 万条谱的第二次加载需要数分钟。另外 `current_spectrum_num` 从 3000 变成了 6000。
  - 建议：维护一个坐标到下标的字典；`load_head_data()` 开始时先清空 `ms` 和计数。

- [ ] **PERF-09 [P2][读码] ion-image finalize 的内存峰值随 `ion_block × 图像面积` 线性增长**
  - 位置：[ion_image_data_manager_h5.py:223-259](src/massflow/data_manager/ion_image_data_manager_h5.py#L223-L259)、[:407-429](src/massflow/data_manager/ion_image_data_manager_h5.py#L407-L429)；[ion_image_data_manager_binary.py:527-556](src/massflow/data_manager/ion_image_data_manager_binary.py#L527-L556)；[msi_zarr/writer.py:599-639](src/massflow/msi_zarr/writer.py#L599-L639)
  - 问题：H5 按稠密的 `H × W` 分配数据块，而且带预取，同一时间有两块数据在内存里。以 2000×2000 的图像、`ion_block=256`、float32 为例，一块约 4 GB，峰值接近 8 GB。Binary 的内存占用随 `merge_workers` 并发数成倍增加，Zarr 同样有双缓冲。
  - 建议：根据内存预算自动推算块大小，并在文档中说明。

- [ ] **PERF-10 [P3][读码] 其他较小的性能问题**
  - `Spectrum.__str__` 用 Python 内置的 `min`/`max` 遍历 numpy 数组（[spectrum.py:132-133](src/massflow/module/spectrum.py#L132-L133)），每次 `load_head_data()` 中的 `inspect_data` 都会用到；
  - UMAP 对 imzML 采样时，每条谱都重新打开一次 `.ibd`（[umap_kmeans.py:338-343](src/massflow/segmentation/umap_kmeans.py#L338-L343)）；
  - Binary 的 `get_ion_image` 每次调用都重新计算散射索引（[ion_image_data_manager_binary.py:603-614](src/massflow/data_manager/ion_image_data_manager_binary.py#L603-L614)）；
  - Zarr 的 `read_ion_image` 每次调用都读取完整的坐标数组（[msi_zarr/reader.py:236](src/massflow/msi_zarr/reader.py#L236)）；
  - `MSDataManagerZarr.storage_mode` 每次访问都重新校验并覆盖 dtype（[ms_data_manager_zarr.py:76-85](src/massflow/data_manager/ms_data_manager_zarr.py#L76-L85)），而 `matrix_generator` 每个 batch 都会访问它；
  - `_local_poly_interp_numba` 对每个查询点都线性扫描一遍节点（[est_noise_numba.py:528-532](src/massflow/preprocess/numba/est_noise_numba.py#L528-L532)）；
  - H5 数据集的 `__getitem__` 多做了一次 `.clone()`（[ion_image_h5_dataset.py:57](src/massflow/morch/ion_image_h5_dataset.py#L57)）；
  - `IonMajorTempStore.append_batch` 会把调用方已经复制过的数组再复制一次。

## 4. 文档与代码不一致

- [ ] **DOC-01 [P1] 文档引用了已经不存在的方法和类名**
  - `load_full_data_from_file()` 在 `ms-data-structures`、`metadata`、`contribution`、`naming-conventions` 的中英文版，以及 [imzml_monkey_patch.py:14](src/massflow/tools/imzml_monkey_patch.py#L14) 中共出现 27 处，而实际的方法名是 `load_head_data()`。
  - [naming-conventions.md:22-29](docs/en/naming-conventions.md#L22-L29) 中列出的 `MSIDataManager`、`MSIDataManagerZYS`、`load_meta_from_file`、`allocate_data_from_meta`、`add_msi_slice`、`get_msi_by_mz`、`plot_msi`、`write2local` 在代码中都不存在。

- [ ] **DOC-02 [P1] peak_align 的默认值与代码不符**
  - [peak_alignment.md:36-40](docs/en/peak_alignment.md#L36-L40)（中文版相同）写的是 `binfun="median"`、`min_frequency=None`，而代码中（[api.py:258-261](src/massflow/preprocess/api.py#L258-L261)）分别是 `"min"` 和 `0.01`。文档还说 `"median"` 是推荐的默认值。

- [ ] **DOC-03 [P1] peak_pick 的文档有三处与代码不符（中英文都一样）**
  - `method` 的默认值：文档写 `"quantile"`，代码是 `"mad"`（[peak_pick.md:37](docs/en/peak_pick.md#L37)）。
  - `width`：文档称它是"半宽"（[:52](docs/en/peak_pick.md#L52)），但代码中 `radius = width // 2`，也就是说 width 实际是整个窗口的宽度。
  - `relheight`：文档说它只在 `return_type="area"` 时用来确定积分边界（[:59-61](docs/en/peak_pick.md#L59-L61)）。代码里它其实是一个过滤条件：峰高与最高峰之比小于 `relheight` 的峰被丢弃。它与 `return_type` 无关，计算面积边界时也不用它。

- [ ] **DOC-04 [P2] noise_reduction 的签名和正文自相矛盾**
  - [noise_reduction.md:33](docs/en/noise_reduction.md#L33) 的签名写的是 `method="ma_numba"`，但第 43 行的正文和代码都是 `gaussian_numba`。

- [ ] **DOC-05 [P2] getting-started 对元数据输出的描述不准确**
  - [getting-started.md:144](docs/en/getting-started.md#L144) 和 [:149](docs/en/getting-started.md#L149) 说 centroid 或 processed 数据会显示为 `False`。实际上这些字段的值是 None，而 `__str__` 会跳过 None，所以对应的行根本不会显示。

- [ ] **DOC-06 [P2] 没有说明任务会被按固定顺序重排**
  - planner 按 [planner.py:20-30](src/massflow/preprocess/pipeline/planner.py#L20-L30) 中的固定顺序执行任务，用户注册的顺序会被忽略。按这个顺序，adduct（80）在 stats（90）之前执行，而 [getting-started.md:334-335](docs/en/getting-started.md#L334-L335) 表格里的顺序正好相反。

- [ ] **DOC-07 [P2] 归一化文档没有说明默认参考 m/z 会随 batch 变化**
  - [normalization.md:126](docs/en/normalization.md#L126) 说不设 `ref` 时取"当前提供的 m/z 数据的中间位置"，但没有说明：对 processed 数据，每个 batch 取到的值都不一样（见 BUG-09）；找不到参考峰的谱会被清零。

- [ ] **DOC-08 [P3] 部分开发者文档不在侧边栏中**
  - `contribution`、`naming-conventions`、`ms-data-structures`、`metadata`、`testing_guide`、`collaboration_guide`、`zh/zarr_schema` 这 7 个页面在侧边栏中没有入口。DOC-01 所说的过时内容正好集中在这几页。另外 `zarr_schema` 只有中文版。

- [ ] **DOC-09 [P3] R 环境文档对 locale 的描述不准确**
  - [r_environment_setup.md](docs/en/r_environment_setup.md) 说"把 R 的 locale 设为英文"，实际修改的是整个 Python 进程的 `LANGUAGE` 和 `LC_ALL`（见 LOGIC-09）。

- [ ] **DOC-10 [P3] 示例代码和 README 过时**
  - [main.py:1](main.py#L1) 的 docstring 写着 "imzML → binary ion image"，实际输出的是 `ion_image_h5`，而且 import 了没有用到的 `IonImageDataManagerBinary`。
  - README 的项目结构只列出了两个 DataManager，还提到了没有纳入版本控制的 `.vscode/`。
  - README 和 README_CN 的安装步骤不一致：中文版多了 `uv sync --group dev`、`pip install -e .` 等步骤。

- [ ] **DOC-11 [P3] 过时的 docstring 和日志文案**
  - [tools/\_\_init\_\_.py:3-5](src/massflow/tools/__init__.py#L3-L5) 提到了一个不存在的 `__init___.py` 文件。
  - [logger.py:218](src/massflow/tools/logger.py#L218) 提示调用 `bootstrap.init_logging()`，这个函数并不存在。
  - [infer_spectrum_type.py:31](src/massflow/tools/infer_spectrum_type.py#L31) 的日志说 "Setting metadata accordingly"，实际上并没有写回元数据。

## 5. 测试、CI 与工程化

- [ ] **TEST-01 [P1] CI 不跑测试，推送到 main 就直接发布到 PyPI**
  - 仓库里只有 [deploy.yml](.github/workflows/deploy.yml)（部署文档）和 [publish-to-pypi.yml](.github/workflows/publish-to-pypi.yml)（发布）两个 workflow，没有任何运行 pytest 或 pylint 的步骤。发布前只检查版本号是否已经存在，即使本地测试失败也照样发布。
  - 建议：增加一个在 PR 上运行 pytest 和 pylint 的 workflow，并让发布 job 依赖它。

- [ ] **TEST-02 [P1] 新 clone 下来的仓库，核心功能几乎没有测试覆盖**
  - 35 个失败/错误中有 34 个、以及全部 64 个跳过，都是因为测试数据（`data/example.imzML`、`baseline_test.imzML`、`example_file.imzML`、`example.iib` 等）不在仓库中；9 个异步流水线测试全部被跳过。
  - [test_ion_image_zarr.py:25](tests/test_ion_image_zarr.py#L25) 硬编码了另一台机器上的路径 `/Users/dre/Desktop/...`。
  - 缺少数据时，有的测试会跳过，有的直接报 `FileNotFoundError`，处理方式不统一。
  - 建议：用合成数据的 fixture 现场生成小的 imzML/Zarr 文件（几十条谱，毫秒级即可生成），用来覆盖本报告中的回归用例；使用真实数据的基准测试单独打上 marker。

- [ ] **TEST-03 [P2] 唯一真实失败的测试已经过时，转换脚本的相应参数也因此失效**
  - `test_mean_spectrum_is_absent_by_default` 断言默认情况下不写 mean_spectrum。但 `IonImageDataManagerZarr` 默认 `include_spectra=True`，writer 会强制写入 mean_spectrum（[msi_zarr/writer.py:77-79](src/massflow/msi_zarr/writer.py#L77-L79)）。这与 [zarr_schema.md:301](docs/zh/zarr_schema.md#L301) 中"combined 文件必须包含 tic 与 mean_spectrum"的规定一致，所以需要修改的是测试。
  - 同样的原因，[scripts/imzml_to_zarr.py:78-82](scripts/imzml_to_zarr.py#L78-L82) 的 `--no-mean-spectrum` 对 `ion_image_zarr` 不起作用。
  - 该脚本的 `summarize()` 会调用 `reader.n_ions`（[:148](scripts/imzml_to_zarr.py#L148)），输入是 processed 数据时，转换虽然成功，却会在这一步报 ValueError。此外，它打印的格式版本号是硬编码的字符串 "1.1"（[:146](scripts/imzml_to_zarr.py#L146)）。

## 6. 代码质量与小问题

- [ ] **MINOR-01 [P3]** 日志 group 名拼错成了 `"preprocesss"`（[est_noise_helper.py:13](src/massflow/preprocess/helper/est_noise_helper.py#L13)），会多生成一个 `logs/preprocesss.log`。另外，`ColorFormatter` 用 `str.replace` 替换级别名，消息正文里出现的 `INFO`、`ERROR` 等字样也会被染色并补齐空格（[logger.py:39](src/massflow/tools/logger.py#L39)）。
- [ ] **MINOR-02 [P3]** `Preprocessor` 用 `batch_size > 9056` 作为上限，这是一个魔数（[preprocessor.py:62](src/massflow/preprocess/preprocessor.py#L62)）。
- [ ] **MINOR-03 [P3]** `filter_data_by_mask` 把 `processed` 设为 False（[stat_filter_helper.py:201](src/massflow/preprocess/helper/stat_filter_helper.py#L201)），而其他地方的约定是设为 None。
- [ ] **MINOR-04 [P3]** 元数据的小问题：`version` 为字符串时，校验会抛 TypeError；把 `mask` 设为 None 会被静默忽略，因此无法清除 mask；描述符上的 `MetaField.__int__` 意义不明（[ms_meta_data.py:69-71](src/massflow/module/ms_meta_data.py#L69-L71)、[:137-163](src/massflow/module/ms_meta_data.py#L137-L163)）。
- [ ] **MINOR-05 [P3]** `risky_imzml_loader` 通过 monkey patch 修改 pyimzml 的模块级函数，在多线程下不安全（[imzml_monkey_patch.py:23-38](src/massflow/tools/imzml_monkey_patch.py#L23-L38)）。
- [ ] **MINOR-06 [P3]** `read_contiguous_block` 不检查实际读到的字节数，文件被截断时报出的是难以理解的 reshape 错误；`is_disk_contiguous` 只比较首尾两个偏移（[read_matrix.py:118-152](src/massflow/tools/read_matrix.py#L118-L152)）。
- [ ] **MINOR-07 [P3]** 死代码：`tools/dm_process.flat_speed_process`、`SpectrumImzML.swap_out2disk`、`MetaDataBase.update_meta` 都没有调用方；`runner._del_last_dm` 中的 `del del_dm` 不起任何作用。
- [ ] **MINOR-08 [P3]** [tools/plot.py](src/massflow/tools/plot.py) 文件开头带有 UTF-8 BOM。
- [ ] **MINOR-09 [P3]** `flat_executor`、`export_imzml`、`umap_kmeans` 等多处直接调用 `MSDataManagerZarr._target_indices()`、`_require_reader()` 等私有方法（并用 pylint disable 压掉警告），建议把它们提升为公开 API。
- [ ] **MINOR-10 [P3][实测]** TIC 归一化的求和在 float32 上累加（[normalization_numba.py:12-14](src/massflow/preprocess/numba/normalization_numba.py#L12-L14)）。200 万个点时，相对误差为 2.9e-4；`_get_rms` 使用 float64 累加，没有这个问题。建议改成与 `_get_rms` 一样的显式 float64 循环。

## 7. 实测确认没有问题的部分

以下各处在审阅时存在疑点，但经脚本实测，结果都是正确的，不需要修改：

- **ion-image Zarr 导出为 imzML 的往返**：分块转置时的偏移计算正确，读回的每个像素、每个离子的强度都与源数据一致。
- **ion-major Zarr 写入时的坐标去重**：重复坐标保留最后一条；去重后，`/ion_image`、`/spectra` 和 `stats/tic` 三者保持对齐。
- **加合物过滤中流式 Pearson 相关系数的合并**：batch 大小不等（包括只有一条谱的 batch）时，Chan 合并公式的结果与 `np.corrcoef` 一致。
- **子区域读取**：imzML 的 `target_locs` 子区域读取，以及 spectra Zarr 在非连续索引下的子区域读取，返回的 m/z、强度和坐标都一一对应。
- **流水线中某个阶段抛出异常时**：异常约 0.1 秒就传回调用方（以 `RuntimeError` 的形式，并附带原始异常），不会等到停滞超时，也不会留下残留的输出文件。
- **峰提取后包含空谱的数据**：写出为 imzML 时正常（只有写出为 Zarr 时会出错，见 BUG-07）。
- **`diff` 噪声估计**：公式 `mean(|x − mean(diff(x))|)` 与 matter 原版的 `estnoise_diff` 一致。
- **R_HOME 的查找优先级**：函数参数 → 环境变量 → `set_default_r_home` → `R RHOME`，与文档描述一致。
- **文档中函数签名的默认值**：用脚本与 `inspect.signature` 自动比对，除 DOC-02、DOC-03、DOC-04 列出的 4 处外全部一致；umap、segmentation、离子图、谱图可视化文档中的参数默认值也与代码一致。
- **`IonImageBinaryDataset` 的 pickle**：访问过之后仍然可以正常 pickle，只有 H5 版本有问题（见 LOGIC-05）。
- **Binary 离子图的坐标处理**：坐标中 y 含 0、x 不含 0 时，Binary 输出正确；H5 在同一组坐标下会崩溃（见 BUG-04）。
- **"最佳实践"文档的分步流程**：每一步都显式传入了不同的 `output_path`，这种写法本身规避了 BUG-03（我没有端到端运行这个流程）。
- **现有测试**：能够运行的 97 个测试全部通过。

## 8. 附录：关键问题复现代码

以下代码在仓库根目录下、使用项目自带的 `.venv` 运行（macOS，commit `a70b5b3`）。这些代码会在当前目录写入文件，并且 import massflow 时会创建 `logs/`（见 BUG-06），建议在临时目录中执行。

**B1. auto 模式写出器写错 m/z（BUG-01）**

```python
import numpy as np
from pyimzml.ImzMLParser import ImzMLParser
from massflow.tools.stream_imzml_writer import ImzMLWriter

w = ImzMLWriter("t", mode="auto")
w.add_spectrum(np.array([100.0, 150.2, 200.0]), np.array([1, 2, 3], np.float32), (1, 1, 1))
w.add_spectrum(np.array([100.0, 170.5, 200.0]), np.array([4, 5, 6], np.float32), (2, 1, 1))
w.close()
print(ImzMLParser("t.imzML").getspectrum(1)[0])  # [100. 150.2 200.]，应为 [100. 170.5 200.]
```

**B2. fastmath 让 NaN 过滤失效（BUG-02）**

```python
import numpy as np
from massflow.preprocess.numba.baseline_correction_numba import _median_of_finite

print(_median_of_finite(np.array([1.0, np.nan, 3.0, 5.0])))  # 4.0，应为 3.0
```

**B3. 连续调用 Preprocessor 读坏自己的输入（BUG-03），以及已注册任务被改写（BUG-10）**

```python
import numpy as np
from massflow.tools.stream_imzml_writer import ImzMLWriter
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import Preprocessor

rng = np.random.default_rng(0)
w = ImzMLWriter("proc", mode="processed", spec_type="centroid")
for i in range(64):
    mz = np.sort(rng.choice(np.linspace(100, 300, 3000), 40, replace=False))
    w.add_spectrum(mz, rng.random(40).astype(np.float32) + 0.1, (i % 8 + 1, i // 8 + 1, 1))
w.close()
dm = MSDataManagerImzML(filepath="proc.imzML")
dm.load_head_data()

# BUG-10：第一次运行得到的参考轴被写回已注册的任务
p = Preprocessor(dm, temp_dir="T0").peak_align()
p.start().close()
task = next(t for t in p._tasks if t.name == "peak_align")
print(task.kwargs["reference"] is None)  # False

# BUG-03：两步都不传 output_path
step1 = Preprocessor(dm, temp_dir="T", batch_size=16).filter_peaks_by_isotopes().start()
step2 = Preprocessor(step1, temp_dir="T", batch_size=16).peak_align().start()
# RuntimeError: Async preprocessing pipeline failed.
#   原始异常：ValueError: could not broadcast input array from shape (0,) into shape (40,)
```

**B4. H5 离子图 finalize 崩溃（BUG-04）**

```python
import numpy as np
from massflow.data_manager import IonImageDataManagerH5

h5 = IonImageDataManagerH5("x.h5", mode="w", width=3, height=2)
coords = np.array([[1, 0], [2, 0], [1, 1], [2, 1]])  # y 从 0 开始，x 不含 0
h5.write_flat_batch(np.arange(8, dtype=np.float32), np.full(4, 2), coords, np.array([100.0, 200.0]))
h5.finalize()  # ValueError: shape mismatch；同一组坐标在 IonImageDataManagerBinary 下输出正确
```

**B5. 临时文件描述符耗尽（BUG-05）**

```python
import resource
import numpy as np
from massflow.data_manager import IonImageDataManagerBinary

_, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (256, hard))  # macOS launchd 的默认软上限

mz = np.linspace(100, 200, 8)
m = IonImageDataManagerBinary("big.iib", mode="w", width=400, height=1)
for b in range(300):
    m.write_flat_batch(np.random.rand(8).astype(np.float32), np.array([8]), np.array([[b + 1, 1]]), mz)
m.finalize()  # OSError: [Errno 24] Too many open files
```

**B6. import 时改动全局日志配置（BUG-06）**

```python
import logging, os

root = logging.getLogger()
h1, h2 = logging.StreamHandler(), logging.StreamHandler()
root.addHandler(h1)
root.addHandler(h2)

import massflow.tools.logger  # noqa: E402

print(os.path.isdir("logs"))                        # True：在当前目录创建了 logs/
print([h for h in root.handlers if h in (h1, h2)])  # 只剩 1 个：删掉了用户的 handler，而且还漏删了一个
```

**B7. 峰提取后存在空谱时，Zarr 输出崩溃（BUG-07）**

```python
import numpy as np
from massflow.tools.stream_imzml_writer import ImzMLWriter
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import Preprocessor

rng = np.random.default_rng(0)
mz = np.linspace(100, 300, 400)
w = ImzMLWriter("cont", mode="continuous", spec_type="profile")
for i in range(32):
    background = i % 4 == 0  # 模拟组织外的全零背景像素
    y = np.zeros(400, np.float32) if background else (
        rng.random(400) + 50 * np.exp(-((mz - 200) / 2) ** 2)
    ).astype(np.float32)
    w.add_spectrum(mz, y, (i % 8 + 1, i // 8 + 1, 1))
w.close()
dm = MSDataManagerImzML(filepath="cont.imzML")
dm.load_head_data()

Preprocessor(dm, temp_dir="T", batch_size=8, output_format="spectra_zarr", output_path="out").peak_pick().start()
# RuntimeError: Async preprocessing pipeline failed.
#   原始异常：ValueError: Writer does not accept empty spectra.
# 同样的调用改为 output_format="imzml" 时正常
```

**B8. SNIP 结果随 batch 大小变化（BUG-08）**

```python
import numpy as np
from massflow.preprocess.numba.baseline_correction_numba import baseline_snip_numba

rng = np.random.default_rng(0)
spec = np.sin(np.linspace(0, 20, 300)) + 2 + rng.random(300)
one = baseline_snip_numba(spec, lengths=np.array([300]))
batch = baseline_snip_numba(np.tile(spec, 64), lengths=np.full(64, 300))
print(np.abs(one - batch[:300]).max())  # 0.17，应为 0
```

**B9. H5 数据集无法 pickle（LOGIC-05）**

```python
import pickle
from massflow.data_manager import IonImageDataManagerH5
from massflow.morch import IonImageH5Dataset

ds = IonImageH5Dataset(IonImageDataManagerH5("ion_images.h5", mode="r"))  # 任意一个已 finalize 的 H5 离子图文件
_ = ds[0]
pickle.dumps(ds)  # TypeError: h5py objects cannot be pickled
```
