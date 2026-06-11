from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any
import uuid


@dataclass(slots=True)
class ToolResult:
    """Standard result envelope for agent-facing MassFlow tools."""

    ok: bool
    tool_name: str
    message: str
    input_path: str | None = None
    job_id: str = field(default_factory=lambda: uuid.uuid4().hex)
    artifacts: dict[str, str] = field(default_factory=dict)
    metrics: dict[str, Any] = field(default_factory=dict)
    details: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @classmethod
    def success(
        cls,
        tool_name: str,
        message: str,
        input_path: str | None = None,
        *,
        artifacts: dict[str, str] | None = None,
        metrics: dict[str, Any] | None = None,
        details: dict[str, Any] | None = None,
        warnings: list[str] | None = None,
    ) -> "ToolResult":
        return cls(
            ok=True,
            tool_name=tool_name,
            message=message,
            input_path=input_path,
            artifacts=artifacts or {},
            metrics=metrics or {},
            details=details or {},
            warnings=warnings or [],
        )

    @classmethod
    def failure(
        cls,
        tool_name: str,
        message: str,
        input_path: str | None = None,
        *,
        errors: list[str] | None = None,
        warnings: list[str] | None = None,
    ) -> "ToolResult":
        return cls(
            ok=False,
            tool_name=tool_name,
            message=message,
            input_path=input_path,
            warnings=warnings or [],
            errors=errors or [],
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
