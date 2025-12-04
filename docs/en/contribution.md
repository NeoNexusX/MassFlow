# Contributing Guidelines

Welcome to contribute to MassFlow! This document provides a unified workflow and quality standards for multi-person collaboration, ensuring organized division of labor, stable code, and controllable releases.

## Quick Start
- After forking the repository, first switch to the remote `dev` branch and create your working branch (see branch naming rules below).
- Install the recommended Cursor/Trae/VSCode extensions locally (see README and `.vscode/extensions.json`). PyCharm users can skip this.
- Install dependencies and development environment:
  ```bash
  # Ensure uv is installed
  uv sync && uv pip install -e .  
  ```
- Follow naming conventions and code style (`NAMING_CONVENTIONS.md` / `NAMING_CONVENTIONS_EN.md`, `.pylintrc`).
- After completing the task, submit a PR, self-check against the review checklist, and request a review.

## Branch Strategy
- `dev`: Main development branch, used for daily development work, merged after review by designated personnel.
  - Feature: `feature/<topic>` (e.g., `feature/imzml-import`)

- `main`: Protected branch, only merged via admin-reviewed PRs; kept stable and releasable.
- Merge policy: Prefer squash merge to keep commit history clean; use rebase when necessary to clean up commit records.

## Submission Workflow (Issue → Branch → PR)
1. Claim Issue: Before starting each new issue, delete the previous repository [including local and online content].
3. Develop and Self-check: Follow code standards, run local checks, and write necessary tests [you can submit a PR first and write test code as requested by the reviewer].
4. Submit PR:
   - Concise title (Conventional Commits recommended): e.g., `feat(msi): add HDF5 group writer`
   - Complete description: key changes, impact scope, testing points (`Closes #<id>`). AI generation is recommended.
   - Request review: designate module owner or code owner.
5. Review and Merge: AI preliminary review → Fix comments → Owner review → Meet checklist → Merge to `dev`.

## Commit Message Guidelines
- Conventional Commits are recommended:
  - `feat: <description>` New feature
  - `fix: <description>` Bug fix
  - `docs: <description>` Documentation update
  - `refactor: <description>` Refactor without affecting function
  - `test: <description>` Test related
- Example: `feat(data-manager): support split/merge write modes`

## Code Standards and Naming
- Follow `NAMING_CONVENTIONS.md` / `NAMING_CONVENTIONS_EN.md`:
  - Classes: `PascalCase`, keep domain acronyms uppercase (e.g., `MSIDataManager`)
  - Functions/variables/files: `snake_case` (e.g., `load_full_data_from_file`)
  - Attributes: Private `_meta_*` mapped to public `meta_*` attributes (via `@property`)
- Unified style and static checks:
  - Run `pylint` (following `.pylintrc`), `black`, `isort`

## Editor and Tools
- Cursor/Trae/VSCode:
  - Install recommended extensions (Python, Pylance, Pylint, Ruff, Black, isort, Jupyter, Markdownlint, GitLens, H5Web)
  - Workspace `.vscode/extensions.json` will automatically prompt for installation
- Optional: Enable Format on Save and configure Black/isort consistency.

## Testing and Quality Assurance
- Test scope:
  - Unit tests: Module methods and critical paths (loading, writing, filtering, visualization)
  - Integration tests: End-to-end path from data file to processing output (lightweight data)
- Run tests:
  ```bash
  uv run pytest
  ```
- Local checks (recommended as a pre-commit step):
  - No errors, 10 rounds of running passed

## PR Review Checklist
- Naming style and module boundaries are consistent; interfaces and data contracts comply with conventions
- No obvious performance or memory issues; avoid unnecessary copies/allocations
- Assertions and error handling are in place; edge cases are covered
- Documentation updates (README/Chinese docs/examples) are in place
- Passed local and CI checks; new/updated tests cover critical logic
- Clear change description

## Contact
- For questions and support needs, please submit an Issue in the repository; or @ the person in charge in the PR for discussion.
