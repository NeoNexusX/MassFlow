from __future__ import annotations

from typing import Any, Iterable

from ._utils import write_result_json
from .chemical_annotation_engine import annotate_mz_values
from .chemical_db_loader import load_chemical_library_directory, summarize_chemical_library
from .chemical_schema import build_chemical_annotation_context
from .schemas import ToolResult


def annotate_mz_with_local_library(
    mz_values: Iterable[float],
    *,
    library_dir: str,
    source_name: str = "auto",
    ion_mode: str = "positive",
    ppm_tolerance: float = 5.0,
    max_candidates_per_mz: int = 5,
    output_dir: str | None = None,
) -> ToolResult:
    """Annotate a list of observed m/z values using local chemical library files."""

    tool_name = "annotate_mz_with_local_library"
    mz_list = [float(value) for value in mz_values]
    try:
        library_entries = load_chemical_library_directory(library_dir, source_name=source_name)
        annotations = annotate_mz_values(
            mz_list,
            library_entries=library_entries,
            ion_mode=ion_mode,
            ppm_tolerance=ppm_tolerance,
            max_candidates_per_mz=max_candidates_per_mz,
        )
    except Exception as exc:  # pylint: disable=broad-except
        return ToolResult.failure(
            tool_name,
            "Failed to annotate m/z values with the local chemical library.",
            errors=[str(exc)],
        )

    details: dict[str, Any] = {
        "library_summary": summarize_chemical_library(library_entries),
        "chemical_context": build_chemical_annotation_context(annotations),
        "annotations": [annotation.to_dict() for annotation in annotations],
    }
    result = ToolResult.success(
        tool_name,
        "Annotated m/z values with the local chemical library.",
        artifacts={},
        metrics={
            "mz_query_count": len(mz_list),
            "annotation_count": len(annotations),
            "library_entry_count": len(library_entries),
        },
        details=details,
    )
    artifact_path = write_result_json(result, output_dir, f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result
