# Preprocessor pipeline parity and executor decoupling plan

## 结论

`src/massflow/preprocess/preprocessor_new.py` 是把旧版 `src/massflow/preprocess/preprocessor.py` 的内联逻辑拆到 `pipeline/` 包后的薄入口。按当前工作区源码静态对齐看，核心处理路径设计上基本等价：任务注册仍来自同一个 `PreprocessorAPI`，profile/centroid 任务过滤规则一致，连续 batch task 合并规则一致，python `peak_align` 仍是先全数据集计算 reference/tolerance、再走 flat batch pipeline，Cardinal `peak_align/peak_pick` 仍走 dataset callable，最终 `imzml` 与 `ion_image_*` 输出逻辑也基本对应。

但当前不能认为二者已经“处理效果完全一致”，原因有三点：

1. `preprocessor_new.py` 现在不可直接导入。`src/massflow/preprocess/pipeline/__init__.py` 仍从 `massflow.preprocess.pipeline.preprocessor` 导入，但该文件在当前工作区已删除，实际代码被移动到了 `src/massflow/preprocess/preprocessor_new.py`。
2. 新入口没有接管公共 API。`src/massflow/preprocess/__init__.py` 仍导出旧版 `Preprocessor`，所以外部 `from massflow.preprocess import Preprocessor` 还不会走新 pipeline。
3. 有少量非核心行为差异需要显式确认：新版本多校验 `output_max_workers > 0`，报错信息和 logger 名称不同；默认 ion image 输出路径由字符串拼接改为 `Path(base_dir) / ...`，Windows 下字符串分隔符可能不同。

因此建议判断为：核心算法链路目标一致，但当前工程状态还没达到可替换旧版的程度，需要先修导入/导出，再做 parity tests。

## dataset_executor 与 flat_executor 的耦合判断

`DatasetTaskExecutor` 现在持有 `FlatPipelineExecutor`，主要是为了处理 python backend 的 `peak_align`：这个任务被注册为 `scope="dataset"`，因为它需要扫描全数据集计算 reference/tolerance；但真正的数据变换仍是 `FlatPreprocess.peak_align_flat`，适合走 flat batch pipeline。

这确实是耦合，但它是“过渡期可理解、长期不理想”的耦合：

- 好处：复用了 flat pipeline，避免在 dataset executor 里重复 batch 读写、队列、sink、错误处理逻辑。
- 问题：dataset executor 同时知道 `FlatPreprocess.peak_align_flat`、sink 创建和 flat executor 调度，职责越界；后续增加类似“先 dataset 准备、再 batch 执行”的任务时会继续堆特殊判断。
- 设计上不应要求 dataset executor 和 flat executor “完全一致”。它们应该是两个不同执行引擎，但要共享一致的执行契约：输入/输出类型、metadata 更新、final output format 处理、异常关闭资源、任务分段语义。

## 推荐方案

### 阶段 1：修正新入口可用性

1. 明确新入口位置：要么恢复 `src/massflow/preprocess/pipeline/preprocessor.py` 作为兼容 shim，要么修改 `pipeline/__init__.py` 指向 `src/massflow/preprocess/preprocessor_new.py`。
2. 保留或补上 `Preprocessor = PipelinePreprocessor` alias，避免外部调用方迁移成本。
3. 决定是否让 `src/massflow/preprocess/__init__.py` 导出新 `Preprocessor`。在 parity tests 通过前，不建议直接替换公共入口。

### 阶段 2：建立旧/新一致性测试

1. 用同一份小型 imzML fixture，对旧版 `Preprocessor` 和新版 `PipelinePreprocessor` 跑相同任务链：
   - `noise_reduction`
   - `baseline_correction`
   - `normalization`
   - `peak_pick` python backend
   - `peak_align` python backend
   - Cardinal backend 的 `peak_pick/peak_align`，可按环境条件跳过
2. 校验内容不要只看类型，应比较：
   - spectrum 数量
   - 每个 spectrum 的 mz/intensity shape
   - 关键 mz/intensity 数值误差
   - `profile_spectrum/centroid_spectrum` metadata
   - `imzml` 与 `ion_image_binary/h5` 输出的 header、mz_axis、n_pixels、n_ions
3. 单独覆盖边界：
   - centroid 输入时 profile-only task 被 drop
   - 所有任务被 drop 后的行为
   - `output_path=None` 的默认输出路径
   - `output_max_workers <= 0` 是否保持旧行为或采用新校验

### 阶段 3：把“dataset 准备 + flat 执行”从 DatasetTaskExecutor 中移走

推荐引入一个轻量的 task lowering/preparation 层，让 executor 只负责执行自己的 engine。

目标结构：

1. `TaskPlanner` 仍负责静态分段：连续 batch task 合并、dataset task 单独成段、final output format 标记。
2. 新增 `TaskLowerer` 或 `TaskPreparationExecutor`，运行时拿到 `data_manager` 后处理需要 dataset 信息的任务。
3. 对 python `peak_align`，lowering 层计算 reference/tolerance，并把原 dataset task 转成一个 batch task。
4. `PipelineRunner` 持有 `FlatPipelineExecutor` 和 `DatasetTaskExecutor`，由 runner 根据 lowering 结果选择执行器。
5. `DatasetTaskExecutor` 只执行真正的 dataset callable，例如 Cardinal `peak_align/peak_pick`，不再 import 或持有 `FlatPipelineExecutor`。

这样可以保留现有处理效果，同时把依赖方向统一为：

`PipelineRunner -> TaskLowerer -> FlatPipelineExecutor / DatasetTaskExecutor`

而不是：

`DatasetTaskExecutor -> FlatPipelineExecutor`

## 任务模型改进建议

当前代码用 `task.name == "peak_align"` 和 `task.apply_fn is FlatPreprocess.peak_align_flat` 判断执行策略，这比较脆。建议后续把 `PreprocessTask` 扩展为显式策略字段，例如：

- `scope`: `batch | dataset`
- `engine`: `flat | dataset`
- `prepare_fn`: 可选，用于 reference/tolerance 等运行时准备
- `output_spectrum_effect`: 例如 `preserve | centroid`

这样 planner/runner 不需要猜 task name，也不需要直接依赖具体函数对象。

## 最小落地顺序

1. 先修导入/导出，让 `PipelinePreprocessor` 可导入但不替换默认 `Preprocessor`。
2. 加旧/新 parity tests，确认实际输出一致。
3. 把 `prepare_flat_peak_align_task` 移到 lowering 层。
4. 移除 `DatasetTaskExecutor` 对 `FlatPipelineExecutor` 和 `OutputSinkFactory` 的依赖。
5. parity tests 通过后，再考虑让公共 `Preprocessor` 指向新 pipeline。
