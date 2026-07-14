# MassFlow 发布前代码审核报告

- 审核日期：2026-07-14
- 审核分支：`feat/filter_by_adduct`
- 审核提交：`b910a5842f3657572cdd2c1326e9125f3fc5aae9`
- 对比基线：`origin/main`（当前分支领先 61 个提交，未落后）
- 当前版本：`0.1.0a3`
- 当前结论：**No-Go，不应从当前状态发布**

## 1. 审核范围与限制

本报告针对当前工作树、当前分支相对默认分支的发布状态、Python 包配置、GitHub Actions 发布流程、当前未提交改动涉及的预处理 API，以及相关依赖、测试和文档入口进行静态审核。

按照本次要求，本轮没有运行 `pytest`、完整测试套件、构建命令、安装命令、文档站点构建或任何会生成发布产物的命令。结论不能替代运行时、真实数据和目标平台验证，也不是对所有算法数值正确性的证明。

严重度定义：

- **P1 / 阻断**：发布前必须解决，否则可能发布错误内容、覆盖已有版本、暴露失效 API，或导致发布/文档流程直接失败。
- **P2 / 高优先级**：建议在本次发行前解决；若延期，必须明确支持边界并记录风险。
- **P3 / 改进项**：不一定阻断功能，但影响发行质量、可维护性或用户体验。

## 2. 总体结论

当前代码库不满足安全发行条件，主要原因如下：

1. 当前版本号仍为已经发布并打过标签的 `0.1.0a3`，发布工作流会把同名产物直接复制到自建索引，存在覆盖已有发行物的风险。
2. 工作树包含 6 个已修改文件和 3 个未跟踪文档；如果现在给 `HEAD` 打标签，这些内容不会进入标签和发行包。
3. VitePress 配置文件末尾存在多余的 `)`，文档配置有明确语法错误。
4. dataset 级流水线任务会无条件改写任务自身的 `output_path`，与公开 API 和新增文档描述冲突。
5. `massflow.morch` 的公开数据集接口运行时使用 PyTorch，但 `torch` 仅在开发依赖组中，安装正式包的用户得不到该依赖或可安装 extra。
6. 发布工作流的验证仅导入顶层包并枚举模块名称，没有导入子模块、校验标签与包版本的一致性，也没有阻止同名产物被覆盖。

在以下 P1 项全部关闭之前，不建议创建新标签或手动触发发布工作流。

## 3. 审核发现

### R-01 [P1] 工作树不是可发行快照

**证据**

当前状态包含以下未提交内容：

- 已修改：
  - `src/massflow/preprocess/api.py`
  - `src/massflow/preprocess/helper/adduct_filter_helper.py`
  - `src/massflow/preprocess/helper/stat_filter_helper.py`
  - `src/massflow/preprocess/pipeline/planner.py`
  - `src/massflow/r_preprocess/adapter.py`
  - `src/massflow/tools/plot_adduct_filter.py`
- 未跟踪：
  - `docs/zh/adduct_filter.md`
  - `docs/zh/isotope_filter.md`
  - `docs/zh/stat_filter.md`

`git describe` 为 `v0.1.0a3-1-gb910a58-dirty`。其中 `Preprocessor.filter_peaks_by_adducts()` 只存在于当前工作树的 `src/massflow/preprocess/api.py:333`，不在已提交的 `HEAD` 中。

**影响**

Git 标签只指向提交。若现在在 `HEAD` 上打标签，新增 API、修正和三份中文文档不会进入发行快照，最终发布内容会与当前本地审核对象不同。

**建议修复**

逐项确认这 9 个文件哪些属于本次发行，完成必要修复后统一提交；不属于本次发行的内容应明确排除，而不是依赖未提交状态保留。

**验收条件**

- 标签前 `git status --short` 无输出。
- 待发布提交包含所有预期代码和文档。
- 审核提交、标签提交和工作流 checkout 的提交完全一致。

### R-02 [P1] 版本号复用会产生同名、不同内容的发行物

**证据**

- `pyproject.toml:7` 仍声明 `version = "0.1.0a3"`。
- 标签 `v0.1.0a3` 已存在，指向提交 `f339cbd04a276f843ccf384487971c297f145328`。
- 当前 `HEAD` 为该标签之后的 `b910a584...`，且还有未提交改动。
- `.github/workflows/publish.yml:49` 使用 `cp dist/* /tmp/pypi-repo/packages/`，同名文件会被直接覆盖。

**影响**

再次发布 `0.1.0a3` 会生成与历史版本同名但内容不同的 wheel/sdist。自建索引中的既有文件和哈希可能被替换，破坏发行物不可变性、缓存一致性和问题追溯。

**建议修复**

- 将下一版本提升为未使用的新版本，例如 `0.1.0a4`，并同步锁文件。
- 标签必须与包版本严格对应，例如 `v0.1.0a4` 对应 `0.1.0a4`。
- 发布工作流增加“同名文件已存在则比较 SHA-256；不同则失败”的保护，不允许静默覆盖。
- 在构建前增加标签名与 `pyproject.toml` 版本一致性校验。

**验收条件**

- 新版本在仓库标签和自建索引中均未使用。
- 标签版本、项目元数据版本、wheel/sdist 文件名一致。
- 同名不同哈希的文件会使工作流失败。

### R-03 [P1] dataset 任务的 `output_path` 被流水线无条件覆盖

**证据**

- `src/massflow/preprocess/api.py:300-331` 和 `:333-358` 分别向 dataset 任务注册方法级 `output_path`。
- `src/massflow/preprocess/pipeline/planner.py:133-138` 在最后一个任务为 dataset 任务且输出为 imZML 时，无条件执行：

  ```python
  task_segments[-1].tasks[0].kwargs["output_path"] = input_meta.result_path
  ```

- `docs/zh/adduct_filter.md:71` 将该参数描述为“该过滤阶段的 imZML 输出路径”。
- `docs/zh/stat_filter.md:61` 同时描述了阶段路径和 `Preprocessor(..., output_path=...)` 的最终路径。

**影响**

当 dataset 过滤任务位于流水线末尾时，调用者显式传给过滤方法的 `output_path` 不一定生效。代码、API 参数和文档目前没有统一的路径优先级合同，可能写入意外位置或返回与预期不一致的结果管理器。

**建议修复**

先确定唯一合同，然后同步实现和文档。建议的最小合同是：

1. 方法级 `output_path` 非空时优先，表示该 dataset 阶段的显式持久化路径。
2. 方法级 `output_path` 为 `None` 时，才由 planner 注入 `Preprocessor` 的最终 `result_path`。
3. 返回的数据管理器和 `output_meta.file_path` 必须指向实际写出的同一路径。

不能只把赋值改成 `setdefault()`，因为任务注册时键已存在且值可能为 `None`；实现应显式判断当前值。

**验收条件**

- 分别覆盖“仅方法路径”“仅 Preprocessor 路径”“两者同时提供”“两者都不提供”四种情况。
- 实际写入路径、返回管理器路径和文档描述一致。
- 显式输出文件属于持久化结果，不被临时文件清理逻辑删除。

### R-04 [P1] 文档配置存在明确语法错误

**证据**

`docs/.vitepress/config.mts:132-133` 当前结尾为：

```ts
})
)
```

最后一个 `)` 没有匹配的左括号。该问题已存在于已提交的 `HEAD`，不是当前未提交改动引入。

**影响**

VitePress 配置无法被正常解析，文档构建/部署流程预期会失败。由于 `.github/workflows/deploy.yml` 只在 `main` 上构建站点，若不在合并前处理，问题会进入默认分支部署链路。

**建议修复**

删除多余右括号，并在修复后由发行人员执行文档配置解析和站点构建验证。

**验收条件**

- 配置文件可被 TypeScript/VitePress 正常解析。
- 文档构建命令退出码为 0。
- 生成站点中的导航和页面链接可访问。

### R-05 [P1] 发布工作流尚未进入默认分支

**证据**

- 当前分支相对 `origin/main` 为 `0 behind / 61 ahead`。
- `origin/main` 中不存在 `.github/workflows/publish.yml`。
- 工作流声明了 `workflow_dispatch`（`.github/workflows/publish.yml:7`）。

GitHub 官方规则指出，`workflow_dispatch` 只有在工作流文件存在于默认分支时才会接收触发事件：<https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#workflow_dispatch>。

**影响**

当前无法把手动发布作为可靠的发行入口。直接在功能分支上推标签虽然可能触发 tag workflow，但会绕过“先进入默认分支、再从可追溯主线发布”的常规门禁，并放大 R-01 和 R-02 的风险。

**建议修复**

先将审核通过的发布工作流和代码合并到默认分支，确认默认分支上的工作流可见，再从默认分支的已审核提交创建标签。

**验收条件**

- 默认分支包含发布工作流。
- 标签提交可从默认分支到达。
- Actions 页面可见手动触发入口；正式发行仍优先使用版本标签触发。

### R-06 [P1] `morch` 运行时依赖 PyTorch，但正式包没有声明可安装依赖

**证据**

- `src/massflow/morch/ion_image_binary_dataset.py:95-139` 在 `__getitem__()` 和 `get_batch()` 中导入并返回 `torch.Tensor`。
- H5 和 Zarr 数据集实现同样在运行时导入 `torch`。
- `pyproject.toml:27-30` 仅定义 `r` extra。
- `torch>=2.11.0` 仅存在于 `pyproject.toml:32-40` 的 `dependency-groups.dev`。

PyPA 的 Dependency Groups 规范明确说明依赖组用于开发环境，不会作为项目构建元数据发布；面向包用户的可选功能应使用 `[project.optional-dependencies]`：

- <https://packaging.python.org/en/latest/specifications/dependency-groups/>
- <https://packaging.python.org/en/latest/guides/writing-pyproject-toml/>

**影响**

普通安装不会提供 PyTorch，也没有 `massflow[torch]` 或类似 extra 可供用户安装。公开导出的 `IonImageBinaryDataset`、`IonImageH5Dataset` 和 `IonImageZarrDataset` 会在实际取样时失败。

**建议修复**

在 `[project.optional-dependencies]` 中增加明确命名的 extra（例如 `torch` 或 `morch`），并在文档中使用 `pip install "massflow[torch]"`。若这些接口被视为核心功能，则应把 PyTorch 放入核心依赖；考虑 PyTorch 的体积和平台差异，独立 extra 更合理。

**验收条件**

- 包元数据中存在用户可安装的 PyTorch extra，或明确将 `morch` 排除在本次正式支持范围外。
- 三种数据集接口的安装说明与实际依赖一致。

### R-07 [P2] 发布验证不能发现子模块导入失败

**证据**

`.github/workflows/publish.yml:30-38` 安装 wheel 后只导入顶层 `massflow`，随后使用 `pkgutil.walk_packages()` 收集模块名称。枚举名称不会导入这些子模块。

**影响**

子模块缺少依赖、导入时语法错误或导出名称错误时，发布验证仍可能显示成功。R-06 和 R-09 这类问题不会被当前步骤可靠发现。

**建议修复**

- 在隔离环境中导入承诺支持的公开模块和公共符号。
- 将可选模块分为“基础安装应成功”和“安装对应 extra 后应成功”两组。
- 增加包版本、标签版本和产物文件名的一致性检查。
- 不要把“枚举到模块名”当作“模块可导入”的证据。

### R-08 [P2] 直接使用了未声明的一层传递依赖

**证据**

`src/massflow/tools/stream_imzml_writer.py:18` 直接导入 `wheezy.template`，但 `pyproject.toml` 没有声明 `wheezy-template`。当前可用性依赖 `pyimzml` 间接安装该包。

**影响**

上游依赖关系变化、裁剪安装或不同解析结果都可能使该工具导入失败。直接使用的第三方库不应依赖另一依赖的实现细节来间接获得。

**建议修复**

把 `wheezy-template` 声明为直接依赖；如果该 writer 是可选工具，则放入对应 extra，并使文档明确安装方式。

### R-09 [P2] 公共 docstring 示例使用不存在的类名

**证据**

- `src/massflow/morch/ion_image_binary_dataset.py:38,42` 使用 `IonImageDatasetBinary`。
- 实际类名和 `src/massflow/morch/__init__.py:4,10` 的公开导出是 `IonImageBinaryDataset`。

**影响**

用户复制示例后会得到 `ImportError`，且示例不能作为公开 API 的有效使用说明。

**建议修复**

将两处示例统一为 `IonImageBinaryDataset`，并把公开示例导入纳入静态或文档示例检查。

### R-10 [P2] 新增中文文档没有接入站点导航

**证据**

- 三份新增文档当前均为未跟踪文件。
- `docs/.vitepress/config.mts:100` 声称侧边栏列出英文和中文页面，但 `:106-120` 只有 English 分组。
- 新增的 adduct、isotope 和 stat filter 页面均未出现在侧边栏。

**影响**

即使文件被提交并修复语法错误，普通用户也难以从文档站点导航发现这些页面；配置注释与实际内容也不一致。

**建议修复**

将确认要发行的中文页面加入中文导航/侧边栏，或建立明确的中文索引页；同时修正配置注释。`socialLinks` 当前指向 `NeoNexusX/MassFlow`，而本地 `origin` 为 `dil666555/MassFlow`，应确认哪个仓库是正式入口后统一链接。

### R-11 [P2] 当前测试文件没有覆盖新增流水线接线合同

**证据**

- `tests/test_adduct_filter.py` 主要覆盖 helper、报告和绘图函数。
- 静态搜索没有发现通过 `Preprocessor.filter_peaks_by_adducts()` 验证 planner 注册、路径优先级和最终返回管理器的用例。
- 当前 R adapter 新增的 baseline/noise/normalization 暴露也未发现相应的流水线级用例。

**影响**

helper 层正确并不能证明公开流水线 API 正确；R-03 正是接线层合同缺口。由于本次按要求没有运行测试，本报告只能确认“测试源中缺少该合同覆盖”，不能确认现有测试的运行结果。

**建议修复**

至少补齐 planner/API 的小型合同用例，重点验证任务注册参数、顺序、路径优先级和返回路径。是否运行 `pytest` 由项目发行政策另行决定，本次审核未执行。

### R-12 [P3] 包元数据和格式门禁不完整

**证据**

- `[project]` 当前没有 SPDX `license`、`license-files`、authors、classifiers 或 `project.urls`。
- 仓库有 Black 开发依赖，但发布/CI 工作流没有格式检查门禁。
- 当前改动的部分文件不符合仓库锁定的 Black 版本；全仓也存在历史格式欠账。

**影响**

该问题通常不直接破坏运行，但会降低包索引信息完整性和发布一致性。对全仓直接格式化会产生大量无关改动，不符合本次最小范围原则。

**建议修复**

- 补充最基本的 license 和项目链接元数据，并与仓库实际许可证保持一致。
- 明确格式策略：只处理本次触及文件，或单独建立格式基线；不要把全仓机械格式化混入功能发行提交。

## 4. 已确认的正向结果

在本次静态审核范围内，以下项目未发现阻断问题：

- 当前分支包含 `origin/main` 的全部提交，没有分叉或落后。
- 当前已跟踪差异通过 `git diff --check`，未发现空白错误。
- `pyproject.toml` 与当前 `uv.lock` 的项目版本均为 `0.1.0a3`；问题是版本需要提升，而不是两者当前互相不一致。
- 当前工作树没有发现明显的凭据、令牌或大型二进制文件进入待发布改动。
- 对当前安装的 Cardinal 3.12.0 接口签名进行静态核对时，新增 adapter 使用的 baseline、smooth 和 normalize 方法名/主要参数未发现明显不匹配。
- 新增中文 Markdown 中的仓库内相对链接静态检查未发现缺失目标。

这些结果不抵消第 3 节的发布阻断项。

## 5. 修复后的发行门禁

以下条件应全部满足后，才能把结论从 No-Go 改为 Go：

- [ ] R-01：工作树干净，标签包含全部预期改动。
- [ ] R-02：使用新版本号，标签/元数据/产物一致，索引禁止覆盖同名不同哈希文件。
- [ ] R-03：明确并验证 dataset 任务的输出路径优先级。
- [ ] R-04：VitePress 配置语法修复并通过文档构建验证。
- [ ] R-05：发布工作流进入默认分支，标签提交可从默认分支到达。
- [ ] R-06：PyTorch 支持边界和可安装依赖声明一致。
- [ ] R-07：发布验证可以发现承诺支持的子模块导入失败。
- [ ] R-08 至 R-11：解决，或由维护者书面确认延期原因、影响范围和后续版本计划。
- [ ] 重新审核最终提交相对本报告提交的差异，避免修复过程引入新问题。

## 6. 问题解决后的发行步骤

以下是修复完成后由发行人员执行的步骤，本次审核没有执行这些命令，也没有生成任何产物。

### 6.1 固化发行快照

1. 确认本次发行范围，只保留已经审核的代码和文档。
2. 将版本从 `0.1.0a3` 提升到未使用的新版本，例如 `0.1.0a4`，同步 `uv.lock`。
3. 提交所有预期文件，确认工作树干净。
4. 把功能分支通过正常评审流程合并到 `main`，不要直接从脏工作树或仅存在于功能分支的提交发布。

建议的静态确认项：

```bash
git status --short
git diff --check v0.1.0a3..HEAD
uv lock --check --offline
```

第一条必须无输出，其余命令必须成功。

### 6.2 验证版本和发布工作流

1. 确认 `main` 包含 `.github/workflows/publish.yml`。
2. 确认 `pyproject.toml`、`uv.lock`、目标标签和产物命名使用同一版本。
3. 确认新标签尚不存在：

   ```bash
   git tag --list v0.1.0a4
   git ls-remote --tags origin refs/tags/v0.1.0a4
   ```

   两条命令在创建标签前都不应返回已有标签。
4. 确认发布工作流已经加入：
   - 标签与项目版本一致校验；
   - 同名产物哈希冲突保护；
   - 基础安装的公共模块导入检查；
   - 安装对应 extra 后的可选模块导入检查。

### 6.3 由发行人员生成并检查发行物

在隔离环境中构建 wheel 和 sdist，随后检查：

- wheel 与 sdist 文件名版本正确；
- 包元数据、许可证和依赖声明正确；
- wheel 中包含所有预期 Python 模块；
- 基础安装不依赖开发依赖组；
- `morch` 按已声明的 extra 安装后可导入并使用；
- sdist 能独立重建 wheel；
- 文档站点可构建，新增中文页面可从导航进入。

本报告不把 `pytest` 列为已执行项目；如项目自己的发行政策要求测试，应在独立验证记录中执行并留存结果。

### 6.4 创建标签并发布

1. 在已经合并并验证的 `main` 提交上创建与版本一致的标签，例如 `v0.1.0a4`。
2. 先推送 `main`，再推送该标签。
3. 观察 `Publish Pre-release Package` 工作流，任何验证失败都应终止发布，不要手工跳过。
4. 确认自建 PyPI 索引新增的是新文件，而不是替换旧版本文件。
5. 确认 GitHub Release 标记为 prerelease，附件版本和哈希与工作流产物一致。

### 6.5 发布后核验

1. 从自建索引查询新版本，确认旧版本仍可访问且哈希未变化。
2. 在干净环境中按基础安装和各可选 extra 分别安装新版本，检查公开导入路径。
3. 核对文档站点、GitHub Release、项目版本说明和安装命令均指向同一版本。
4. 保存标签提交、工作流运行链接、产物 SHA-256 和发布结论，作为可追溯发行记录。

## 7. 最终发行建议

当前建议是：**停止发布，先完成 R-01 至 R-07 的修复和验收，再评估 R-08 至 R-11 是否全部纳入本次发行。**

最短安全路径为：修复文档配置和 planner 路径合同，补齐运行时依赖声明与发布工作流保护，提交全部预期改动，提升到新版本，合并到 `main`，再由发行人员按照第 6 节生成、验证并发布新产物。
