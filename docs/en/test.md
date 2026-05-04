# Preprocessor 重构方案:修复 spectrum-type 元数据 bug + 厘清 dataset/batch 任务边界

## Context

`src/massflow/preprocess/preprocessor.py` 在 commit `8c392a4`(Feature pytorch #76)合并时,把语义完全不同的两个参数错误地替换:

| 参数 | 类型 | 职责 |
|------|------|------|
| `output_spectrum_type` | `Literal["profile", "centroid"]` | 写入 `set_spectrum_type_metadata`,标记结果是 profile 还是 centroid |
| `final_output_format` | `str \| None`(`"imzml"` / `"ion_image_binary"` / `"ion_image_h5"`) | 选择最终输出**文件格式**,控制 `StageCSink` 的 kind |

合并后的代码:
- `_run_flat_batch_task`([preprocessor.py:317-393](src/massflow/preprocess/preprocessor.py#L317-L393)):`output_spectrum_type` 形参被改为 `final_output_format`(职责完全不同),并且 `set_spectrum_type_metadata(...)` 调用整体被删除 → batch 路径丢失元数据(回归 bug)。
- `_run_dataset_task`([preprocessor.py:438-463](src/massflow/preprocess/preprocessor.py#L438-L463)):形参更名为 `final_output_format`,但 line 460 仍引用旧名 `output_spectrum_type` → 运行时 `NameError`。
- 测试 `tests/test_async_pipeline.py::TestAsyncPipelineProfile::test_pipeline_resolve_spectrum_type` 会失败:运行 `noise_reduction → peak_pick → peak_align` 后断言 `centroid_spectrum is True`,但当前代码不再写入元数据。

同时,当前 `start()` 方法逻辑不够清晰:
1. 嵌套 while 循环把连续 batch 任务"边走边收集",而不是先分组。
2. `fmt = self.output_format if (_ion_out and ...) else None` 在两处重复。
3. `_run_dataset_task` 承担两种职责:(a) 真正的 cardinal 后端 dataset 操作;(b) 给 `peak_align(flat backend)` 计算 reference 后**反过来委托** `_run_flat_batch_task`。
4. `_export_ion_image` 作为兜底分支挂在 `start()` 末尾,与主循环里的 ion_image 输出路径并存。

## 目标

1. 还原并修正 spectrum-type 元数据写入逻辑 — 与 `final_output_format` 解耦,作为独立参数 `output_spectrum_type` 重新加入。
2. 让 `_run_dataset_task` 只处理**真正的 dataset 任务**(cardinal 后端的 `peak_pick` / `peak_align`)。
3. 把 `peak_align(python backend)` 的 reference 预计算从 `_run_dataset_task` 抽到 `start()`,在那里转换为 batch 任务再交给 `_run_flat_batch_task`。
4. 用显式 `_Segment` 列表 + 单层 for 循环替换嵌套 while,消除 `fmt` 重复计算。
5. spectrum-type 元数据**只**在 imzml 返回路径设置(IonImage 输出仅接受 continuous 数据,不需要 profile/centroid 标记)。

## 设计要点

### 1. 参数职责分离

`_run_flat_batch_task` 与 `_run_dataset_task` 的签名都接收**两个**正交参数:

```python
def _run_flat_batch_task(
    self, *,
    data_manager: MSDataManager,
    tasks: list[PreprocessTask],
    output_spectrum_type: SpectrumType,           # 新增:写入元数据的 profile/centroid
    final_output_format: str | None = None,       # 已有:选择 sink 文件格式
) -> MSDataManagerImzML | IonImageDataManagerBinary | IonImageDataManagerH5: ...

def _run_dataset_task(
    self, *,
    data_manager: MSDataManager,
    task: PreprocessTask,
    output_spectrum_type: SpectrumType,           # 新增
) -> MSDataManagerImzML: ...
```

注意:`_run_dataset_task` 重构后**不再**接收 `final_output_format` — 因为 cardinal 后端只会返回 `MSDataManagerImzML`,后续若需要 ion_image 转换,统一由 `start()` 在所有段执行完后通过 `_export_ion_image` 处理。

### 2. 元数据写入位置(单一规则)

只在两处写入 `set_spectrum_type_metadata`:

- `_run_flat_batch_task`:仅当 sink kind 是 `"imzml"` 时,在 `sink.finalize()` 返回的 `MSDataManagerImzML` 上调用一次。
  - 当 sink kind 是 `"ion_image"` 时,跳过(IonImage 不持有 `.ms.meta`)。
- `_run_dataset_task`:总是在结果(必为 `MSDataManagerImzML`)的 `.ms.meta` 上调用一次。

### 3. Segment 分组(显式)

新增模块级数据类:

```python
@dataclass(slots=True)
class _Segment:
    kind: Literal["batch", "dataset"]
    tasks: list[PreprocessTask]   # dataset 段长度恒为 1
```

新增方法 `_build_segments`:**只做分类**,不解析 peak_align 的 reference(reference 必须用循环里**当时的** `current_dm`,不能提前)。

```python
def _build_segments(self, ordered_tasks: list[PreprocessTask]) -> list[_Segment]:
    segments: list[_Segment] = []
    batch_buffer: list[PreprocessTask] = []
    for task in ordered_tasks:
        if task.scope == "batch":
            batch_buffer.append(task)
            continue
        # dataset 任务来了 — 先把缓冲的 batch 任务封段
        if batch_buffer:
            segments.append(_Segment(kind="batch", tasks=batch_buffer))
            batch_buffer = []
        segments.append(_Segment(kind="dataset", tasks=[task]))
    if batch_buffer:
        segments.append(_Segment(kind="batch", tasks=batch_buffer))
    return segments
```

### 4. peak_align(flat backend) 的处理

约定:在 `start()` 段循环里,如果当前段是 `dataset` 且 task 是 `peak_align + flat backend`,就**就地转换**:用 `current_dm` 调用 `_prepare_peak_align_task` 得到 prepared_task,然后把它当作一个**单任务的 batch 段**走 `_run_flat_batch_task`。`_run_dataset_task` 不再处理这种情况。

这保持现有行为(reference 用前序处理后的中间数据计算),也让 `_run_dataset_task` 只剩 cardinal 后端单一职责。

### 5. start() 新控制流(伪码)

```python
def start(self) -> Any:
    current_dm: MSDataManager = self.data_manager
    input_st = self._resolve_input_spectrum_type()
    ordered = self._arrange_tasks(input_st)
    final_st = self._resolve_output_spectrum_type(input_spectrum_type=input_st, tasks=ordered)
    ion_out = self.output_format in self._ION_IMAGE

    if not ordered:
        # 空管线 + ion_image 输出 → 直接走转换
        if ion_out:
            return self._export_ion_image(current_dm)
        return current_dm

    segments = self._build_segments(ordered)
    last_idx = len(segments) - 1

    for i, segment in enumerate(segments):
        is_last = (i == last_idx)
        # 只有最后一段且要 ion_image 输出时,把 fmt 传给 batch executor
        fmt = self.output_format if (ion_out and is_last and segment.kind == "batch") else None
        # 段执行后的 spectrum 类型 = 截至此段累计 tasks 推断
        cumulative_tasks = [t for s in segments[: i + 1] for t in s.tasks]
        seg_st = self._resolve_output_spectrum_type(input_spectrum_type=input_st, tasks=cumulative_tasks)

        if segment.kind == "batch":
            next_dm = self._run_flat_batch_task(
                data_manager=current_dm,
                tasks=segment.tasks,
                output_spectrum_type=seg_st,
                final_output_format=fmt,
            )
        else:  # dataset
            task = segment.tasks[0]
            # peak_align(flat) 就地转换为 batch
            if task.name == "peak_align" and task.apply_fn is FlatPreprocess.peak_align_flat:
                prepared = self._prepare_peak_align_task(data_manager=current_dm, task=task)
                next_dm = self._run_flat_batch_task(
                    data_manager=current_dm,
                    tasks=[prepared],
                    output_spectrum_type=seg_st,
                    final_output_format=fmt,
                )
            else:
                # cardinal 后端 dataset 任务
                next_dm = self._run_dataset_task(
                    data_manager=current_dm,
                    task=task,
                    output_spectrum_type=seg_st,
                )

        self._maybe_close(current_dm)
        if fmt is not None:
            return next_dm   # batch 段已直出 ion_image
        current_dm = cast(MSDataManager, next_dm)

    # 所有段跑完后,若 ion_image 输出但最后一段是 dataset(cardinal),才需要二次转换
    if ion_out and not isinstance(current_dm, (IonImageDataManagerBinary, IonImageDataManagerH5)):
        out = self._export_ion_image(current_dm)
        self._maybe_close(current_dm)
        return out

    return cast(MSDataManagerImzML, current_dm)
```

## 关键文件与修改点

- [src/massflow/preprocess/preprocessor.py](src/massflow/preprocess/preprocessor.py)
  - 新增模块级 `_Segment` dataclass(放在 `PreprocessTask` 之后)
  - `_run_flat_batch_task`(317-393):增加 `output_spectrum_type` 参数;`sink.finalize()` 后,若 sink.kind == "imzml",调用 `set_spectrum_type_metadata(finalized.ms.meta, output_spectrum_type)` 再返回
  - `_run_dataset_task`(438-463):删除 peak_align(flat) 的委托分支(交给 start());移除 `final_output_format` 形参;在 `result` 的 meta 上调用 `set_spectrum_type_metadata(result.ms.meta, output_spectrum_type)`(取代 line 460 的 `NameError`);改为 `if task.name in ("peak_align", "peak_pick"): ... else raise NotImplementedError`
  - 新增 `_build_segments`
  - 重写 `start()`(517-560):按上述伪码,使用 `_build_segments` + 单层 for 循环

## 复用既有工具(避免新写)

- `_resolve_input_spectrum_type` / `_resolve_output_spectrum_type`(167-174)— 直接调用,不动
- `set_spectrum_type_metadata` / `SpectrumType`(`src/massflow/tools/infer_spectrum_type.py`)
- `_prepare_peak_align_task`(395-436)— 不动,仅调用方从 `_run_dataset_task` 改到 `start()`
- `_make_ion_image_dm` / `_export_ion_image` / `_default_ion_image_output_path`(465-511)— 不动
- `StageCSink`(45-78)— 不动
- `_arrange_tasks`(146-161)— 不动

## 风险与缓解

| 风险 | 处理 |
|------|------|
| `_run_flat_batch_task` 中,sink kind 在 `__init__` 后内部已知,需要在 finalize 后做 `if sink.kind == "imzml"` 判断设置元数据。`StageCSink.kind` 字段已存在 → 直接读取即可 | OK |
| ion_image 输出 + 最后一段是 batch:批管线已直接产出 IonImage,跳过元数据写入(IonImage 无 `.ms.meta`)| 在 `_run_flat_batch_task` 内根据 `sink.kind` 判断,自然短路 |
| ion_image 输出 + 最后一段是 cardinal dataset:dataset 段返回 imzml,主循环结束后再走 `_export_ion_image` 转换 | 主循环末尾的 `if ion_out and not isinstance(...)` 兜底分支 |
| peak_align 是**第一个**任务(没有前序 batch)的情形:`current_dm = self.data_manager`,`_prepare_peak_align_task(self.data_manager, task)` — 与现有行为一致 | OK,无回归 |
| `_resolve_output_spectrum_type` 在每段都重新计算 O(n²) | 任务数 n 极小(通常 ≤5),代价可忽略,胜在显式可读 |
| 公开 API(链式接口、`start()` 返回类型)| 不变 |

## 验证方案

1. **静态:** 重读修改后的 `preprocessor.py`,确认:
   - `output_spectrum_type` 不再有 `NameError` 风险(grep 全文,只在两个 executor 内使用,作为参数传入)
   - `_run_dataset_task` 体内不再有 `_run_flat_batch_task` 调用
   - `start()` 中 `fmt =` 只出现一次

2. **单元测试:** 跑 `tests/test_async_pipeline.py`(尤其 `test_pipeline_resolve_spectrum_type`):
   ```bash
   uv run pytest tests/test_async_pipeline.py -v
   ```
   断言通过:`centroid_spectrum is True` / `profile_spectrum is None`(此前因 batch 路径丢失元数据写入而失败)。

3. **集成:** 跑完整 `tests/test_async_pipeline.py` 套件确认 benchmark 与任务排序行为不变。

4. **回归:** 跑相关测试:
   ```bash
   uv run pytest tests/test_swap.py tests/test_main_plot.py -v
   ```

5. **手动 smoke(可选):** 用 `data/example_file.imzML` 跑一个 `noise_reduction → peak_pick → peak_align(python)` 链,导出 imzml 与 ion_image_binary 各一份,确认:
   - imzml 输出的 meta 标记为 centroid
   - ion_image_binary 文件能被读取,数据形状正确

## 实施顺序

1. 调整 `_run_flat_batch_task` 签名,增加 `output_spectrum_type`,在 finalize 后按 sink.kind 写元数据
2. 调整 `_run_dataset_task`:删除 peak_align flat 委托分支,移除 `final_output_format` 形参,修复 line 460
3. 新增 `_Segment` dataclass + `_build_segments` 方法
4. 重写 `start()` 主控制流
5. 跑测试套件验证
