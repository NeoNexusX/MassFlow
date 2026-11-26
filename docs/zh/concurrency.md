# 并发控制

本文介绍如何通过 `massflow.concurrency` 设置 MassFlow 的并发预算，以及预算在数据读取、预处理、写出和 UMAP/KMeans 计算之间的分配方式。

作者：董乐强，董京伟

**问题反馈：** 如果在使用 MassFlow 或阅读文档时发现问题，可通过 [GitHub Issues](https://github.com/NeoNexusX/MassFlow/issues) 提交反馈，或联系邮箱：707244007@qq.com。

## 概述

MassFlow 在处理数据时，会使用读取线程、Numba 计算线程和写入线程。并发预算用于协调这些工作，避免每个线程池都按整台机器的核心数独立运行。

默认情况下，并发控制已启用，初始预算为检测到的可用 CPU 容量。通常可以直接运行预处理或分割流程；需要为其他程序保留计算资源时，可在运行前通过 `set_budget` 调整预算。

- **预算基数**：当前 Python 进程中用来独立计算各池额度的线程数。启动时的自动配置以检测到的 CPU 容量为基数。
- **处理阶段**：例如批次预处理、整数据集任务或 UMAP。不同阶段使用不同的分配比例。
- **资源池**：读取、计算、写入和 OpenMP 对应的额度类别。它们不一定各自拥有独立的线程池。

预算限制的是受管理工作的并发额度，不能直接等同于系统监视器中显示的进程线程总数。流水线协调线程、事件循环和空闲的底层线程仍可能存在。

## 使用示例

在预处理前设置预算，业务调用方式如下：

```python
import massflow.concurrency as mc
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor

mc.set_budget(8)
print(mc.get_status())

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()

normalized = (
    Preprocessor(dm, output_format="spectra_zarr", output_path="normalized.zarr")
    .normalization(method="tic_numba")
    .start()
)

normalized.close()
dm.close()
```

此处 `8` 是配置值；如果进程只能使用 4 个逻辑核，实际预算为 4。后续调用 UMAP 或 KMeans 时也使用同一份全局配置。

也可以在启动 Python 前设置初始预算：

```bash
export MASSFLOW_BUDGET=8
```

`MASSFLOW_BUDGET` 在并发模块首次导入时读取，应为正整数；无效值会被忽略并回退到自动检测。导入后需要调整预算时，使用 `set_budget`。

## 功能说明

### set_budget（设置预算基数与资源池上限）

```python
mc.set_budget(8, pools={"compute": 4})
```

- 功能：为后续处理阶段设置预算基数和可选的资源池上限。
- 参数：
  - `budget`：正整数或 `None`。省略时保留现有预算；传入 `None` 时关闭预算基数并清除已有的池级设置。
  - `pools`：字典，可使用 `compute`、`io_read`、`io_write` 和 `openmp`。正整数设置该池上限，`None` 移除该池的设置；未列出的池保持原配置。
- 输出：与 `get_status()` 相同的状态字典。

常见配置方式：

```python
mc.set_budget(8)                          # 设置预算基数，保留已有池级设置
mc.set_budget(pools={"compute": 4})       # 只调整计算池上限
mc.set_budget(pools={"compute": None})    # 移除计算池上限
mc.set_budget(None)                       # 关闭预算基数并清除所有池级设置
mc.set_budget(mc.get_logical_cpu_count()) # 重新按当前可用 CPU 容量启用
```

有预算基数时，池级设置只能收紧该阶段的分配额，不能增加额度，空出的额度也不会重新分配。`set_budget()` 不带参数会保留当前配置，不会恢复默认值。

如果在同一次调用中指定 `budget=None` 和 `pools`，会先清除旧配置，再应用新的池级设置。例如 `mc.set_budget(None, pools={"compute": 4})` 仍会限制计算池，但不再限制预算基数。要完全关闭控制，使用不带 `pools` 的 `mc.set_budget(None)`，并且不要为后续 `phase` 单独指定池级上限。

配置应在任务之间调整。任一阶段执行期间调用 `set_budget` 会抛出 `RuntimeError`；无效池名、非整数或非正数也会报错，且不会部分修改配置。

### get_status（查看配置与当前阶段）

```python
status = mc.get_status()
print(status)
```

| 字段 | 含义 |
|------|------|
| `budget` | 配置的预算基数；`None` 表示未设置预算基数 |
| `pool_overrides` | 当前全局池级设置 |
| `phase` | 调用上下文所在的阶段名；阶段外为 `None` |
| `allocations` | 当前阶段的资源池分配；阶段外为空字典 |
| `effective_budget` | 经可用 CPU 容量限制后的实际预算；未设置预算基数时为 `None` |
| `sequential` | 当前阶段是否采用顺序执行模式；阶段外为 `False` |

该函数返回调用上下文的状态，不是所有工作线程的实时监控。若要查看阶段分配，应在该阶段内部调用。

## 预算怎样分配

### 资源池与处理阶段

| 资源池 | 主要用途 |
|--------|----------|
| `compute` | Numba 并行计算，以及 UMAP 的计算线程额度 |
| `io_read` | imzML 分片读取等 IO 工作的读取份额 |
| `io_write` | 临时写出、转置和合并等 IO 工作的写入份额 |
| `openmp` | KMeans 等 OpenMP 计算的线程额度 |

| 阶段 | 使用位置 | 资源池倍率 |
|------|----------|----------|
| `preprocess_pipeline` | `Preprocessor.start()` 的批次处理段 | 读取 10% / 计算 100% / 写入 10% |
| `dataset_task` | 峰对齐、统计过滤等整数据集任务 | 计算 80% / 读取 20% |
| `io_transfer` | 预处理中将 imzML 中间结果转为离子图存储 | 读取 20% / 写入 80% |
| `io_read` | 预处理输入校验，或显式指定的读取流程 | 读取 100% |
| `compute` | `flat_speed_process` 等计算入口 | 计算 100% |
| `umap` | UMAP 降维 | 计算 100% |
| `kmeans` | KMeans 聚类 | OpenMP 100% |

`run_umap_kmeans` 依次使用 `umap` 和 `kmeans` 阶段；`plot_umap_image` 只使用 `umap` 阶段。`io_write` 是资源池名称，不是可传给 `phase` 的阶段名。

各池额度等于实际预算乘以本池倍率。倍率最高按 100% 计算，结果向下取整，每个池至少分到 1 个工作线程，再应用池级上限。各池独立计算，额度之和可以超过预算基数。设置 `set_budget(10)` 且至少有 10 个可用逻辑核时，批次预处理的默认分配如下。

```text
io_read = 1, compute = 10, io_write = 1
```

### 读取、写入与 Zarr 共享 IO 额度

启用预算基数，或阶段包含 `io_read` / `io_write` 分配时，imzML 分片读取、临时写出、binary 合并、H5/Zarr 后台预取，以及 Zarr 的本地 IO 和编解码工作，共用一个阶段 IO 执行器。

- 允许读写重叠时，共享 IO 容量为 `io_read + io_write`；只有一个 IO 池时使用该池额度。
- 顺序模式下，容量取读写额度的最大值。
- 有预算基数但阶段未分配 IO 份额时，保留 1 个辅助 IO worker。

上述预算为 10 的例子分配 10 个计算线程和 2 个共享 IO 工作线程。Zarr、其他读写工作和预取共用同一个 IO 执行器。

使用共享 IO 时，Zarr 的异步并发数与共享 IO 容量一致，事件循环使用同一个执行器。阶段退出时恢复原有 Zarr 设置和执行器，并关闭本阶段创建的执行器；借用者不会提前关闭共享资源。

### 小预算下的顺序执行

实际预算小于阶段中的资源池数量时，会标记为顺序模式。例如批次预处理预算为 1–2 时，切换为“读取一批 → 计算一批 → 写出一批”，完成后再读下一批。

对批次预处理而言，此时计算池先分到全部实际预算，再应用池级上限；读写各保留 1 个 worker 并分时执行，不能把显示的各池额度相加理解为同时工作的线程数。该路径仍按批次处理数据，不会因预算减小而加载整个数据集。实际预算为 1 时，H5/Zarr 最终写出关闭后台预取。

### CPU 容量与底层计算库

实际预算为配置值与可用 CPU 容量的较小值。容量检测考虑逻辑核数、平台支持的进程 CPU affinity（允许使用的 CPU 集合），以及 Linux cgroup v1/v2 配额。默认预算基数按可用逻辑核计算，不会自动缩减为性能核数。

Numba 通过调用线程的线程掩码限制并行计算，生效值还受 Numba 固定线程池大小限制。UMAP 的 `n_jobs` 使用计算额度；设置 `random_state` 时，当前 UMAP 会将 `n_jobs` 改为 1，因此提高预算不一定提高其拟合速度。

有预算基数时，BLAS 限制为 1 个线程；OpenMP 在 `kmeans` 阶段使用其分配额，在其他阶段限制为 1。关闭预算基数但指定了 `openmp` 额度时，也会施加相应的 OpenMP 限制和单线程 BLAS 限制。

## 使用范围与注意事项

- **同一进程中的多个任务**：受控的顶层阶段依次执行，各自使用阶段预算；串行范围是阶段，不是整个 `Preprocessor.start()` 调用。
- **独立读取或写入**：部分接口会查询全局额度，但直接调用数据管理器不会自动建立阶段。需要完整的阶段限制时，使用 `phase` 包住实际工作。
- **运行时设置恢复**：受控阶段正常或异常退出时，都会恢复其修改的 Numba 掩码、原生库限制和 Zarr 设置。原生库与 Zarr 设置具有进程级影响，独立进程、任意外部库调用及 R/Cardinal 内部线程不由阶段准入管理。
- **自定义线程与嵌套阶段**：内部线程需通过 `contextvars.copy_context()` 继承上下文；执行 Numba 内核的线程还需应用计算掩码，并在阶段退出前结束工作。只有进入阶段的线程可以嵌套阶段，嵌套前应先等待外层后台任务完成；继承上下文的工作线程不能再次进入阶段。
- **异步调用**：如果 asyncio 任务在 `await` 期间仍持有受控阶段，同一线程上的其他任务必须等该阶段退出后，再进入自己的阶段。
- **停滞与退出**：流水线可以检测停滞，但不能强制中断已经运行的原生计算或阻塞 IO；退出仍需等待工作线程结束。

预算、批次大小和存储布局分别影响并发、内存和 IO 行为。增大预算不保证提速，比较性能时应保持数据、批次大小、输出格式和缓存条件一致。

## 扩展接口

以下接口主要用于接入自定义执行流程。常规预处理和分割通常只需 `set_budget` 与 `get_status`。

| 接口 | 用途 |
|------|------|
| `get_effective_budget()` | 返回当前实际预算基数；未设置预算基数时为 `None` |
| `runs_sequentially()` | 查询当前阶段是否采用顺序模式 |
| `get_allocation(pool, default=None)` | 查询池额度，不创建阶段；有预算基数的阶段内，未分配的辅助池返回 1 |
| `resolve_workers(pool, requested, default=...)` | 将请求值限制在池额度内；没有额度时使用请求值或默认值 |
| `apply_compute_mask()` | 将计算额度应用到调用线程，返回生效值 |
| `apply_mask(target)` / `restore_mask(previous)` | 设置并返回原 Numba 掩码 / 恢复原掩码 |
| `resolve_compute_workers(target)` | 将计算线程数限制在可用 CPU 与 Numba 池上限内 |
| `get_logical_cpu_count()` | 返回进程可用的逻辑 CPU 容量 |
| `detect_performance_core_workers()` | 按性能核数量估算计算线程数，并受可用容量与 Numba 池大小限制；不设置 CPU 绑定 |

`apply_compute_mask()` 只有在没有可用计算额度时才回退到性能核检测。直接调用掩码接口不会自动建立阶段或恢复掩码，自定义代码需管理其作用范围。

## 参考

- [`concurrency/controller.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/controller.py)：预算配置与阶段分配
- [`concurrency/runtime.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/runtime.py)：共享 IO 执行器与运行时限制
- [`concurrency/numba_bridge.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/numba_bridge.py)、[`concurrency/capacity.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/concurrency/capacity.py)：Numba 掩码与 CPU 容量检测
