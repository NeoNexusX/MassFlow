from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from massflow.runtime_config import find_runtime_config_path, load_runtime_config

from .registry import register_massflow_tools


def build_server(**kwargs: Any) -> Any:
    """Create the FastMCP server and register the stable MassFlow tools."""
    try:
        from mcp.server.fastmcp import FastMCP
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Missing optional dependency 'mcp'. Install it with "
            "`pip install -e .[mcp]` before starting the MassFlow MCP server."
        ) from exc

    server = FastMCP(
        name="MassFlow",
        instructions=(
            "Use MassFlow tools to inspect imzML datasets, run lightweight QC, "
            "apply stable preprocessing, and orchestrate the current non-msiPL workflow."
        ),
        **kwargs,
    )
    register_massflow_tools(server)
    return server


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="massflow-mcp-server",
        description="Run the MassFlow MCP server.",
    )
    parser.add_argument(
        "--config",
        help="Optional path to massflow.toml. When omitted, the current directory is searched.",
    )
    parser.add_argument(
        "--transport",
        choices=["stdio", "sse", "streamable-http"],
        help="Transport for MCP server. Use stdio for MCP clients, streamable-http for local debugging.",
    )
    parser.add_argument(
        "--host",
        help="Host for HTTP transports.",
    )
    parser.add_argument(
        "--port",
        type=int,
        help="Port for HTTP transports.",
    )
    parser.add_argument(
        "--mount-path",
        help="Optional mount path for SSE transport.",
    )
    parser.add_argument(
        "--print-claude-config",
        action="store_true",
        help="Print a ready-to-paste Claude Code MCP config snippet and exit.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    config = load_runtime_config(args.config)
    transport = args.transport or config.mcp_server.transport
    host = args.host or config.mcp_server.host
    port = args.port or config.mcp_server.port
    mount_path = args.mount_path if args.mount_path is not None else config.mcp_server.mount_path

    if args.print_claude_config:
        sys.stdout.write(
            json.dumps(
                build_claude_code_config(config_path=find_runtime_config_path(args.config)),
                indent=2,
                ensure_ascii=False,
            )
        )
        sys.stdout.write("\n")
        return 0

    try:
        server = build_server(host=host, port=port)
    except ModuleNotFoundError as exc:
        sys.stderr.write(f"{exc}\n")
        return 2

    server.run(transport=transport, mount_path=mount_path)
    return 0


def build_claude_code_config(*, config_path: str | None | object = None) -> dict[str, Any]:
    args = ["-m", "massflow.mcp_server.server"]
    if config_path:
        args.extend(["--config", str(config_path)])
    return {
        "mcpServers": {
            "massflow": {
                "command": sys.executable,
                "args": args,
            }
        }
    }


if __name__ == "__main__":
    raise SystemExit(main())
