# Contributing Guidelines

Thank you for contributing to MassFlow! This document defines the shared workflow and quality standards for collaborating, so work stays organized, code remains stable, and releases are controlled.

## Quick Start
- Fork the repository or branch off `main` for your work branch (recommended naming below).
- Install the recommended VSCode extensions locally (see README and `.vscode/extensions.json`).
- Follow naming conventions and code style (`NAMING_CONVENTIONS.md` / `NAMING_CONVENTIONS_EN.md`, `.pylintrc`).
- When your task is done, open a PR, self-check with the review checklist, and request review.

## Branch Strategy
- `main`: Protected branch. Only PRs that pass review are merged; keep it stable and releasable.
- Working branches (short-lived, purpose-specific):
  - Feature: `feature/<topic>` (e.g., `feature/imzml-import`)
- Merge policy: Prefer squash merge to keep history clean; use rebase as needed to tidy commits.

## Workflow (Issue → Branch → PR)
1. Create an Issue: describe background, goals, and acceptance criteria (DoD).
2. Claim and break down: confirm assignee and labels (e.g., `type/feature`, `area/data-manager`, `priority/high`).
3. Develop and self-check: follow code standards; run local checks and required tests.
4. Open a PR:
   - Concise title (Conventional Commits recommended), e.g., `feat(msi): add HDF5 group writer`.
   - Complete description: key changes, impact scope, testing points (`Closes #<id>`).
   - Request reviewers: tag module owners or code owners.
5. Review and merge: address comments → CI green → meets checklist → merge into `main`.

## Commit Message Guidelines
- Conventional Commits are recommended:
  - `feat: <description>` new feature
  - `fix: <description>` bug fix
  - `docs: <description>` documentation update
  - `refactor: <description>` refactor without behavior change
  - `test: <description>` tests
- Example: `feat(data-manager): support split/merge write modes`

## Code Style and Naming
- Follow `NAMING_CONVENTIONS.md` / `NAMING_CONVENTIONS_EN.md`:
  - Classes: `PascalCase`; keep domain acronyms uppercase (e.g., `MSIDataManager`).
  - Functions/variables/files: `snake_case` (e.g., `load_full_data_from_file`).
  - Metadata: private `_meta_*` mapped to public `meta_*` via `@property`.
  - HDF5: group names `mz_{:.5f}`; datasets `mz`, `msroi`; metadata datasets `meta_*`.
- Consistent style and static checks:
  - Run `ruff`, `pylint` (per `.pylintrc`), `black`, `isort`.
  - Clear assertions and error messages: `assert condition, "message"`.

## Editor and Tools
- VSCode:
  - Install recommended extensions (Python, Pylance, Pylint, Ruff, Black, isort, Jupyter, Markdownlint, GitLens, H5Web).
  - The workspace `.vscode/extensions.json` will prompt you to install them automatically.
- Optional: enable Format on Save and configure Black/isort consistency.

## Testing and Quality Assurance
- Test scope:
  - Unit tests: module methods and critical paths (load, write, filter, visualize).
  - Integration tests: end-to-end from data files to processed outputs (use lightweight data).
- Local checks (recommended before committing):
  - Ensure no errors.

## PR Review Checklist
- Naming style and module boundaries are consistent; interfaces and data contracts comply with conventions.
- No obvious performance or memory issues; avoid unnecessary copies/allocations.
- Assertions and error handling are adequate; edge cases covered.
- Documentation updates (README/Chinese docs/examples) are complete.
- Passes local and CI checks; new/updated tests cover critical logic.
- Clear change description.

## Issue Templates and Automation (Recommended)
- Templates:
  - `ISSUE_TEMPLATE/feature.md` (background, approach, acceptance criteria).
  - `ISSUE_TEMPLATE/bug.md` (repro steps, expected/actual, environment info).
- Auto reviewer assignment: use `CODEOWNERS` to designate module owners.
- Local hooks: configure `pre-commit` to auto-format and lint before commits.
- CI: run `ruff/pylint/black/isort` and tests on PRs.

## Release and Versioning
- Semantic Versioning: `MAJOR.MINOR.PATCH`; bump `MAJOR` for public API or data format changes.
- Changelog: maintain `CHANGELOG.md` with features, fixes, and breaking changes.
- Release flow: merge to `main` → tag → generate Release Notes → update docs.

## Contact
- For questions or support, open an Issue in the repository; or @ module owners in PR discussions.