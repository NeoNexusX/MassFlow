"""Contracts between the agent skills in .claude/skills and mzorbit's tool catalog.

The skills name mzorbit tools, their parameters and parameter values in prose.
When a tool or parameter is renamed, nothing else fails: the skill just starts
telling the model to call something that does not exist. These tests keep the
two in step, and check that the skills, their references and the MCP server
instructions point at things that exist.
"""

import re
from pathlib import Path

import pytest

catalog = pytest.importorskip("mzorbit.catalog")

ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = ROOT / ".claude" / "skills"
SKILL_DIRS = sorted(p.parent for p in SKILLS_DIR.glob("*/SKILL.md"))
TOOLS = catalog.TOOL_CATALOG_BY_NAME

# Prefixes shared by mzorbit tool names; a backticked identifier with one of
# them is read as a tool reference.
TOOL_PREFIXES = (
    "load_", "set_", "get_", "profile_", "render_", "plan_", "try_", "compare_",
    "preview_", "finalize_", "run_", "apply_", "export_", "cluster_", "score_",
    "search_", "resolve_", "annotate_",
)
# Result fields that share a tool prefix (e.g. details.finalize_step).
RESULT_FIELDS = {"finalize_step"}

SPAN = re.compile(r"`([^`]+)`", re.S)
CALL = re.compile(r"([a-z]+(?:_[a-z0-9]+)+)\((.*)\)\s*$", re.S)
KEYWORD = re.compile(r"(?<![=!<>*\w])([a-z_][a-z0-9_]*)=(?!=)(\"[^\"]*\"|[^,)]*)")
SKILL_NAME = re.compile(r"\bmsi-[a-z]+(?:-[a-z]+)*(\.md)?")

pytestmark = pytest.mark.skipif(not SKILL_DIRS, reason="no agent skills in .claude/skills")


def _skill_files():
    return sorted(SKILLS_DIR.rglob("*.md"))


def _spans(path):
    """Backticked spans of a markdown file, whitespace-normalised."""
    return [" ".join(m.group(1).split()) for m in SPAN.finditer(path.read_text(encoding="utf-8"))]


def _frontmatter(path):
    match = re.match(r"---\n(.*?)\n---\n", path.read_text(encoding="utf-8"), re.S)
    assert match, f"{path} has no YAML frontmatter"
    fields = {}
    for line in match.group(1).splitlines():
        key, _, value = line.partition(":")
        fields[key.strip()] = value.strip().strip('"')
    return fields


def _tool_references(path):
    """(tool name, keyword arguments) for every tool the file names."""
    refs = []
    for span in _spans(path):
        call = CALL.match(span)
        if call:
            refs.append((call.group(1), KEYWORD.findall(call.group(2))))
        elif re.fullmatch(r"[a-z]+(?:_[a-z0-9]+)+", span) and span.startswith(TOOL_PREFIXES):
            if span not in RESULT_FIELDS:
                refs.append((span, []))
    return refs


@pytest.mark.parametrize("skill_dir", SKILL_DIRS, ids=lambda p: p.name)
def test_frontmatter_names_the_skill_and_fits_the_description_limit(skill_dir):
    fields = _frontmatter(skill_dir / "SKILL.md")
    assert fields.get("name") == skill_dir.name
    # opencode skips skills whose name breaks this pattern.
    assert re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", skill_dir.name) and len(skill_dir.name) <= 64
    description = fields.get("description", "")
    assert description, "description is the trigger: it must not be empty"
    assert len(description) <= 1024, f"description has {len(description)} characters (limit 1024)"


@pytest.mark.parametrize("path", _skill_files(), ids=lambda p: str(p.relative_to(SKILLS_DIR)))
def test_skills_call_existing_tools_with_existing_parameters(path):
    problems = []
    for name, keywords in _tool_references(path):
        spec = TOOLS.get(name)
        if spec is None:
            problems.append(f"unknown tool {name!r}")
            continue
        params = {param.name: param for param in spec.params}
        for keyword, value in keywords:
            param = params.get(keyword)
            if param is None:
                problems.append(f"{name}: unknown parameter {keyword!r}")
            elif param.choices and value.startswith('"'):
                if value.strip('"') not in param.choices:
                    problems.append(f"{name}: {keyword}={value} not in {param.choices}")
    assert not problems, "\n".join(problems)


def test_every_tool_is_covered_by_a_skill():
    mentioned = {name for path in _skill_files() for name, _ in _tool_references(path)}
    uncovered = sorted(set(TOOLS) - mentioned)
    assert not uncovered, f"tools no skill tells the model how to use: {uncovered}"


@pytest.mark.parametrize("path", _skill_files(), ids=lambda p: str(p.relative_to(SKILLS_DIR)))
def test_skill_and_reference_links_resolve(path):
    skill_dir = SKILLS_DIR / path.relative_to(SKILLS_DIR).parts[0]
    text = path.read_text(encoding="utf-8")
    problems = []
    for match in SKILL_NAME.finditer(text):
        if not match.group(1) and not (SKILLS_DIR / match.group(0)).is_dir():
            problems.append(f"unknown skill {match.group(0)!r}")
    for span in _spans(path):
        if span.endswith(".md") and " " not in span:
            bases = (path.parent, skill_dir, ROOT)
            if not any((base / span).resolve().is_file() for base in bases):
                problems.append(f"missing file {span!r}")
    assert not problems, "\n".join(problems)


def test_mcp_instructions_name_existing_tools_and_every_skill():
    pytest.importorskip("mcp.server.fastmcp")
    from mzorbit.interfaces.mcp_server import build_server

    instructions = build_server().instructions
    words = set(re.findall(r"\b[a-z]+(?:_[a-z0-9]+)+\b", instructions))
    unknown = sorted(w for w in words if w.startswith(TOOL_PREFIXES) and w not in TOOLS)
    assert not unknown, f"instructions name unknown tools: {unknown}"
    named = {m.group(0) for m in SKILL_NAME.finditer(instructions) if not m.group(1)}
    assert named == {p.name for p in SKILL_DIRS}
