from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import os
import tomllib


DEFAULT_CONFIG_CANDIDATES = ("massflow.toml", ".massflow.toml")


@dataclass(slots=True)
class AgentToolDefaults:
    default_output_dir: str | None = None
    max_threads: int = 8
    qc_sample_size: int = 16
    preprocess_qc_sample_size: int = 12
    noise_method: str = "sd"
    batch_size: int = 256


@dataclass(slots=True)
class MCPServerDefaults:
    transport: str = "stdio"
    host: str = "127.0.0.1"
    port: int = 8000
    mount_path: str | None = None


@dataclass(slots=True)
class MassFlowRuntimeConfig:
    config_path: str | None = None
    agent_tools: AgentToolDefaults = field(default_factory=AgentToolDefaults)
    mcp_server: MCPServerDefaults = field(default_factory=MCPServerDefaults)


def load_runtime_config(
    config_path: str | None = None,
    *,
    search_start: str | Path | None = None,
) -> MassFlowRuntimeConfig:
    """Load optional runtime config from TOML while preserving current defaults."""
    resolved_path = find_runtime_config_path(config_path, search_start=search_start)
    if resolved_path is None:
        return MassFlowRuntimeConfig()

    raw = tomllib.loads(resolved_path.read_text(encoding="utf-8"))
    config = MassFlowRuntimeConfig(config_path=str(resolved_path))

    _apply_mapping(config.agent_tools, raw.get("agent_tools"), allowed_keys={
        "default_output_dir",
        "max_threads",
        "qc_sample_size",
        "preprocess_qc_sample_size",
        "noise_method",
        "batch_size",
    })
    _apply_mapping(config.mcp_server, raw.get("mcp_server"), allowed_keys={
        "transport",
        "host",
        "port",
        "mount_path",
    })
    return config


def find_runtime_config_path(
    config_path: str | None = None,
    *,
    search_start: str | Path | None = None,
) -> Path | None:
    """Resolve the runtime config path from explicit path, env var, or cwd."""
    if config_path:
        path = Path(config_path).expanduser().resolve()
        return path if path.exists() else None

    env_path = os.environ.get("MASSFLOW_CONFIG")
    if env_path:
        path = Path(env_path).expanduser().resolve()
        return path if path.exists() else None

    base = Path(search_start) if search_start is not None else Path.cwd()
    for candidate in DEFAULT_CONFIG_CANDIDATES:
        path = (base / candidate).resolve()
        if path.exists():
            return path
    return None


def _apply_mapping(target: Any, raw_section: Any, *, allowed_keys: set[str]) -> None:
    if not isinstance(raw_section, dict):
        return

    for key, value in raw_section.items():
        if key in allowed_keys and hasattr(target, key):
            setattr(target, key, value)
