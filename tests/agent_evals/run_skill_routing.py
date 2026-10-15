"""Skill-routing eval: which MSI skill does Claude load for a request?

Each query in skill_routing.json runs through ``claude -p`` in the repo, with
the real skills in .claude/skills and the mzorbit MCP server, and the first
decisive action is recorded: a Skill call (which skill), another tool call, or
an answer without tools. The three skills compete with each other and with the
mzorbit tools themselves, so routing, not isolated triggering, is what matters.

Safety: read-only lookups (Read, Glob, Grep, ToolSearch, TodoWrite) may run;
the process is killed as soon as any other tool call *starts*, before its input
is complete, so no tool with side effects executes. Writing tools are also
disallowed outright.

    python tests/agent_evals/run_skill_routing.py --claude <path-to-claude> \
        --model claude-opus-5-5 --runs 3
"""
from __future__ import annotations

import argparse
import json
import os
import select
import shutil
import subprocess
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
READ_ONLY = {"Read", "Glob", "Grep", "ToolSearch", "TodoWrite"}
DISALLOWED = [
    "Edit", "Write", "NotebookEdit",
    "mcp__mzorbit__finalize_preprocessing", "mcp__mzorbit__apply_tissue_mask",
    "mcp__mzorbit__export_dataset", "mcp__mzorbit__set_biometa", "mcp__mzorbit__set_acquisition",
]
MAX_READ_ONLY = 4  # a skill is chosen up front; long read-only exploration means none was


def _short(tool: str) -> str:
    return tool.replace("mcp__mzorbit__", "mzorbit:")


def run_query(claude: str, query: str, *, model: str | None, mcp_config: str, timeout: int) -> dict:
    """Run one query; return its outcome and the tool trail that led to it."""
    cmd = [
        claude, "-p", query,
        "--output-format", "stream-json", "--verbose", "--include-partial-messages",
        "--no-session-persistence", "--mcp-config", mcp_config, "--strict-mcp-config",
    ]
    if model:
        cmd += ["--model", model]
    cmd += ["--disallowedTools", *DISALLOWED]
    env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, cwd=ROOT, env=env)
    trail: list[str] = []
    current: list[str] | None = None  # [tool name, partial input json]
    buffer, start = b"", time.time()

    def done(outcome: str) -> dict:
        return {"outcome": outcome, "trail": trail, "seconds": round(time.time() - start, 1)}

    try:
        while time.time() - start < timeout:
            ready, _, _ = select.select([proc.stdout], [], [], 1.0)
            if not ready:
                if proc.poll() is not None:
                    break
                continue
            chunk = os.read(proc.stdout.fileno(), 65536)
            if not chunk:
                break
            buffer += chunk
            while b"\n" in buffer:
                line, buffer = buffer.split(b"\n", 1)
                try:
                    event = json.loads(line)
                except ValueError:
                    continue
                if event.get("type") == "result":
                    return done("answered")
                if event.get("type") != "stream_event":
                    continue
                inner = event.get("event", {})
                kind = inner.get("type")
                block = inner.get("content_block", {})
                if kind == "content_block_start" and block.get("type") == "tool_use":
                    name = block.get("name", "")
                    trail.append(_short(name))
                    if name != "Skill" and name not in READ_ONLY:
                        return done(f"tool:{_short(name)}")  # killed before its input completes
                    current = [name, ""]
                elif kind == "content_block_delta" and current:
                    delta = inner.get("delta", {})
                    if delta.get("type") == "input_json_delta":
                        current[1] += delta.get("partial_json", "")
                elif kind == "content_block_stop" and current:
                    name, raw = current
                    current = None
                    try:
                        args = json.loads(raw) if raw else {}
                    except ValueError:
                        args = {}
                    if name == "Skill":
                        return done(f"skill:{args.get('skill')}")
                    path = str(args.get("file_path", ""))
                    if name == "Read" and path.endswith("SKILL.md"):
                        return done(f"skill:{Path(path).parent.name}")
                    trail[-1] = f"{name}({path or args.get('pattern') or args.get('query') or ''})"
                    if sum(t.split('(')[0] in READ_ONLY for t in trail) >= MAX_READ_ONLY:
                        return done("read-only exploration")
        return done("timeout")
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()


def passed(expected: str, outcome: str) -> bool:
    if expected == "none":
        return not outcome.startswith("skill:msi-")
    return outcome == f"skill:{expected}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--eval-set", default=str(Path(__file__).with_name("skill_routing.json")))
    parser.add_argument("--claude", default=shutil.which("claude") or "claude")
    parser.add_argument("--model", default=None)
    parser.add_argument("--mcp-config", default=str(ROOT / ".mcp.json"))
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=150)
    parser.add_argument("--only", default=None, help="Run only queries whose expected route equals this")
    parser.add_argument("--out", default=None, help="Write full results as JSON here")
    args = parser.parse_args()

    items = json.loads(Path(args.eval_set).read_text(encoding="utf-8"))
    if args.only:
        items = [item for item in items if item["expected"] == args.only]
    jobs = [(i, item) for i, item in enumerate(items) for _ in range(args.runs)]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        runs = list(pool.map(
            lambda job: (job[0], run_query(args.claude, job[1]["query"], model=args.model,
                                           mcp_config=args.mcp_config, timeout=args.timeout)),
            jobs,
        ))

    results = []
    for i, item in enumerate(items):
        outcomes = [run for j, run in runs if j == i]
        rate = sum(passed(item["expected"], run["outcome"]) for run in outcomes) / len(outcomes)
        results.append({**item, "pass_rate": rate, "pass": rate >= 0.5, "runs": outcomes})

    for expected in ("msi-dataset-evaluation", "msi-preprocessing", "msi-spatial-analysis", "none"):
        group = [r for r in results if r["expected"] == expected]
        if not group:
            continue
        print(f"\n== {expected}: {sum(r['pass'] for r in group)}/{len(group)} queries pass")
        for r in group:
            seen = Counter(run["outcome"] for run in r["runs"])
            mark = "PASS" if r["pass"] else "FAIL"
            print(f"  [{mark}] {dict(seen)}  {r['query'][:70]}")
    total = sum(r["pass"] for r in results)
    print(f"\nTOTAL {total}/{len(results)} queries pass "
          f"({sum(r['pass_rate'] for r in results) / len(results):.0%} of runs routed correctly)")
    if args.out:
        Path(args.out).write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
