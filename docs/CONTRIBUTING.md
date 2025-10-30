# 贡献指南（Contributing Guidelines）

欢迎为 MassFlow 做出贡献！本文档为多人协作提供统一的流程与质量规范，确保分工有序、代码稳定、发布可控。

## 快速开始
- Fork 仓库或从 `dev` 切出你的工作分支（推荐分支命名见下文）。
- 在本地安装推荐的 Cursor/Trae/VSCode 插件（见 README 与 `.vscode/extensions.json`）。Pycharm 可以略过
- 遵循命名规范与代码风格（`NAMING_CONVENTIONS.md` / `NAMING_CONVENTIONS_EN.md`、`.pylintrc`）。
- 完成任务后提交 PR，按照评审清单自检并请求评审。

## 分支策略
- `dev`：主要开发分支，用于日常开发工作。
- `main`：受保护分支，仅通过评审的 PR 合入；保持稳定可发布。
  - 功能：`feature/<topic>`（例：`feature/imzml-import`）
- 合并策略：优先使用 squash merge 保持提交历史简洁；必要时 rebase 清理提交记录。

## 提交流程（Issue → 分支 → PR）
1. 创建 Issue：描述背景、目标与验收标准（DoD）
2. 认领与拆解：确认标签（如 `type/feature`, `area/data-manager`, `priority/high`）。
3. 开发与自检：遵循代码规范，运行本地检查与必要测试。
4. 提交 PR：
   - 标题精炼（建议遵循 Conventional Commits）：例如 `feat(msi): add HDF5 group writer`
   - 描述完整：关键改动、影响范围、测试要点（`Closes #<id>`）
   - 请求评审：指定模块负责人或代码所有者
5. 评审与合并：修复评论项 → CI 通过 → 满足评审清单 → 合并至 `dev`。

## 提交信息规范（Commit Messages）
- 推荐使用 Conventional Commits：
  - `feat: <描述>` 新功能
  - `fix: <描述>` Bug 修复
  - `docs: <描述>` 文档更新
  - `refactor: <描述>` 重构但不影响功能
  - `test: <描述>` 测试相关
- 示例：`feat(data-manager): support split/merge write modes`

## 代码规范与命名
- 遵循 `NAMING_CONVENTIONS.md` / `NAMING_CONVENTIONS_EN.md`：
  - 类：`PascalCase`，领域缩写保持大写（如 `MSIDataManager`）
  - 函数/变量/文件：`snake_case`（如 `load_full_data_from_file`）
  - 元数据：私有 `_meta_*` 与公开 `meta_*` 属性映射（通过 `@property`）
  - HDF5：组名 `mz_{:.5f}`；数据集 `mz`、`msroi`；元数据数据集 `meta_*`
- 统一风格与静态检查：
  - 运行 `ruff`、`pylint`（遵循 `.pylintrc`）、`black`、`isort`
  - 断言与错误信息明确：`assert condition, "message"`

## 编辑器与工具
- Cursor/Trae/VSCode：
  - 安装推荐扩展（Python、Pylance、Pylint、Ruff、Black、isort、Jupyter、Markdownlint、GitLens、H5Web）
  - 工作区 `.vscode/extensions.json` 将自动提示安装
- 可选：启用保存自动格式化（Format on Save）并配置 Black/isort 一致性。

## 测试与质量保障
- 测试范围：
  - 单元测试：模块方法与关键路径（加载、写入、过滤、可视化）
  - 集成测试：从数据文件到处理输出的端到端路径（轻量数据）
- 本地检查（建议作为提交前步骤）：
  - 无报错

## PR 评审检查清单（Checklist）
- 命名风格与模块边界一致，接口与数据约定符合规范
- 无明显性能与内存问题；避免不必要拷贝/分配
- 断言与错误处理到位，边界情况覆盖
- 文档更新（README/中文文档/示例）到位
- 通过本地与 CI 检查；新增/更新测试覆盖关键逻辑
- 变更说明清晰

## Issue 模板与自动化（建议）
- 模板：
  - `ISSUE_TEMPLATE/feature.md`（背景、方案、验收标准）
  - `ISSUE_TEMPLATE/bug.md`（复现步骤、预期/实际、环境信息）
  -  可以直接选用对应模板使用
- 自动评审指派：通过 `CODEOWNERS` 指定模块负责人
- 本地钩子：配置 `pre-commit` 在提交前自动进行格式化与 lint
- CI：在 PR 上运行 `ruff/pylint/black/isort` 与测试

## 发布与版本管理
- 语义化版本：`MAJOR.MINOR.PATCH`；对外接口与数据格式变更提升 `MAJOR`
- 变更日志：维护 `CHANGELOG.md`，记录功能、修复、破坏性变更
- Release 流程：合并到 `dev` → 打 Tag → 产出 Release Notes → 同步文档

## 联系方式
- 如有问题与支持需求，请在仓库提交 Issue；或在 PR 中 @负责人 进行讨论。