from __future__ import annotations

from dataclasses import asdict
from importlib.util import find_spec
from pathlib import Path
from typing import Any
import json
import os
import platform
import sys
import tempfile

from massflow.runtime_config import load_runtime_config

from ._utils import make_json_safe
from .schemas import ToolResult


def run_self_check(
    input_path: str | None = None,
    *,
    output_dir: str | None = None,
    config_path: str | None = None,
    require_mcp: bool = False,
) -> ToolResult:
    """Run lightweight environment and file checks without touching core workflows."""
    config = load_runtime_config(config_path)
    effective_output_dir = output_dir if output_dir is not None else config.agent_tools.default_output_dir

    warnings: list[str] = []
    errors: list[str] = []

    python_ok = sys.version_info >= (3, 12)
    if not python_ok:
        errors.append("Python 3.12 or newer is required.")

    core_deps = {
        "numpy": _module_available("numpy"),
        "scipy": _module_available("scipy"),
        "pyimzml": _module_available("pyimzml"),
    }
    if not all(core_deps.values()):
        missing = [name for name, available in core_deps.items() if not available]
        errors.append(f"Missing required runtime dependencies: {', '.join(missing)}")

    mcp_available = _module_available("mcp")
    if require_mcp and not mcp_available:
        errors.append("Missing optional dependency 'mcp'. Install it with `pip install -e .[mcp]`.")
    elif not mcp_available:
        warnings.append("Optional dependency 'mcp' is not installed; Claude Code MCP startup will be unavailable.")

    writable_probe = _check_writable_directory(effective_output_dir)
    if not writable_probe["ok"]:
        errors.append(writable_probe["message"])

    input_probe = _check_input_pair(input_path) if input_path else None
    if input_probe is not None:
        warnings.extend(input_probe["warnings"])
        errors.extend(input_probe["errors"])

    result = ToolResult(
        ok=not errors,
        tool_name="run_self_check",
        message="MassFlow self-check completed." if not errors else "MassFlow self-check found issues.",
        input_path=input_path,
        metrics={
            "python_version": platform.python_version(),
            "core_dependency_count": len(core_deps),
            "core_dependency_ready_count": sum(1 for ready in core_deps.values() if ready),
            "mcp_available": mcp_available,
        },
        details={
            "runtime": {
                "python_executable": sys.executable,
                "python_version_info": list(sys.version_info[:3]),
                "platform": platform.platform(),
                "cwd": os.getcwd(),
            },
            "config": {
                "config_path": config.config_path,
                "agent_tools": make_json_safe(asdict(config.agent_tools)),
                "mcp_server": make_json_safe(asdict(config.mcp_server)),
            },
            "dependencies": {
                "core": core_deps,
                "optional": {
                    "mcp": mcp_available,
                },
            },
            "output_probe": writable_probe,
            "input_probe": input_probe,
            "usage_modes": {
                "claude_stdio": {
                    "description": "Use Claude Code with MCP config; the client auto-starts the server.",
                    "command": [sys.executable, "-m", "massflow.mcp_server.server"],
                },
                "manual_debug_http": {
                    "description": "Use local HTTP transport for manual debugging only.",
                    "command": [
                        sys.executable,
                        "-m",
                        "massflow.mcp_server.server",
                        "--transport",
                        "streamable-http",
                        "--host",
                        config.mcp_server.host,
                        "--port",
                        config.mcp_server.port,
                    ],
                },
            },
        },
        warnings=warnings,
        errors=errors,
    )
    return result


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="massflow-doctor",
        description="Run local environment checks for MassFlow agent tools and MCP startup.",
    )
    parser.add_argument("--input", help="Optional .imzML input path used to verify the imzML/ibd pair.")
    parser.add_argument("--output-dir", help="Optional writable output directory to probe.")
    parser.add_argument("--config", help="Optional path to massflow.toml.")
    parser.add_argument("--require-mcp", action="store_true", help="Fail when the optional MCP dependency is unavailable.")
    args = parser.parse_args()

    result = run_self_check(
        args.input,
        output_dir=args.output_dir,
        config_path=args.config,
        require_mcp=args.require_mcp,
    )
    sys.stdout.write(json.dumps(make_json_safe(result.to_dict()), indent=2, ensure_ascii=False))
    sys.stdout.write("\n")
    return 0 if result.ok else 1


def _module_available(name: str) -> bool:
    return find_spec(name) is not None


def _check_writable_directory(output_dir: str | None) -> dict[str, Any]:
    target = Path(output_dir) if output_dir else Path(tempfile.gettempdir()) / "massflow_agent_probe"
    try:
        target.mkdir(parents=True, exist_ok=True)
        probe = target / ".massflow_write_test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink(missing_ok=True)
        return {"ok": True, "path": str(target), "message": "Writable output directory verified."}
    except OSError as exc:
        return {"ok": False, "path": str(target), "message": f"Output directory is not writable: {exc}"}


def _check_input_pair(input_path: str) -> dict[str, Any]:
    path = Path(input_path)
    warnings: list[str] = []
    errors: list[str] = []

    if not path.exists():
        errors.append(f"Input imzML file does not exist: {path}")
        return {
            "exists": False,
            "imzml_path": str(path),
            "ibd_path": str(path.with_suffix('.ibd')),
            "warnings": warnings,
            "errors": errors,
        }

    if path.suffix.lower() != ".imzml":
        warnings.append("Input path does not end with .imzML; ensure this is the correct MSI metadata file.")

    ibd_path = path.with_suffix(".ibd")
    if not ibd_path.exists():
        errors.append(f"Missing paired .ibd file: {ibd_path}")

    return {
        "exists": True,
        "imzml_path": str(path),
        "ibd_path": str(ibd_path),
        "warnings": warnings,
        "errors": errors,
    }
