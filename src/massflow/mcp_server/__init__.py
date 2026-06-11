"""MCP server entrypoints for exposing stable MassFlow agent tools."""

__all__ = ["build_server", "main"]


def __getattr__(name: str):
    if name in {"build_server", "main"}:
        from .server import build_server, main

        return {
            "build_server": build_server,
            "main": main,
        }[name]
    raise AttributeError(name)
