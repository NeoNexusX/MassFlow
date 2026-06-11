from __future__ import annotations

from pathlib import Path
from typing import Any

from ._utils import ensure_directory, write_result_json
from .chemical_tools import annotate_mz_with_local_library
from .conclusion_engine import build_integrated_conclusion
from .evidence_tools import build_unified_evidence
from .preprocess_tools import run_preprocessing
from .preprocessing_recommender import recommend_preprocessing_plan
from .qc_tools import dataset_qc_summary, load_imzml
from .schemas import ToolResult


def analyze_imzml_end_to_end(
    imzml_path: str,
    *,
    run_preprocessing_first: bool = True,
    preprocessing_operations: list[dict[str, Any]] | None = None,
    biological_metadata: dict[str, Any] | None = None,
    biological_metadata_json: str | None = None,
    biological_metadata_path: str | None = None,
    chemical_library_dir: str | None = None,
    chemical_source_name: str = "auto",
    chemical_ion_mode: str | None = None,
    chemical_ppm_tolerance: float = 5.0,
    chemical_max_mz_values: int = 12,
    chemical_max_candidates_per_mz: int = 5,
    output_dir: str | None = None,
) -> ToolResult:
    """Run the current end-to-end bootstrap workflow for one imzML dataset."""
    tool_name = "analyze_imzml_end_to_end"
    warnings: list[str] = []
    output_root = ensure_directory(output_dir, temp_prefix="massflow_workflow")

    load_result = load_imzml(
        imzml_path,
        biological_metadata=biological_metadata,
        biological_metadata_json=biological_metadata_json,
        biological_metadata_path=biological_metadata_path,
        output_dir=str(output_root / "load"),
    )
    if not load_result.ok:
        return ToolResult.failure(
            tool_name,
            "End-to-end workflow failed during load stage.",
            imzml_path,
            errors=load_result.errors,
            warnings=load_result.warnings,
        )

    qc_result = dataset_qc_summary(
        imzml_path,
        biological_metadata=biological_metadata,
        biological_metadata_json=biological_metadata_json,
        biological_metadata_path=biological_metadata_path,
        output_dir=str(output_root / "qc"),
    )
    if not qc_result.ok:
        return ToolResult.failure(
            tool_name,
            "End-to-end workflow failed during QC stage.",
            imzml_path,
            errors=qc_result.errors,
            warnings=qc_result.warnings,
        )

    warnings.extend(load_result.warnings)
    warnings.extend(qc_result.warnings)
    biological_context = load_result.details.get("biological_context") or load_result.details.get(
        "dataset", {}
    ).get("biological_context")
    statistical_evidence = qc_result.details.get("statistical_evidence")
    resolved_chemical_library_dir = _resolve_chemical_library_dir(chemical_library_dir)
    resolved_chemical_ion_mode = chemical_ion_mode or _infer_ion_mode(biological_context)

    analysis_input_path = imzml_path
    preprocess_result: ToolResult | None = None
    preprocessing_recommendation_result = recommend_preprocessing_plan(
        imzml_path,
        biological_metadata=biological_metadata,
        biological_metadata_json=biological_metadata_json,
        biological_metadata_path=biological_metadata_path,
        output_dir=str(output_root / "preprocess_recommendation"),
    )
    warnings.extend(preprocessing_recommendation_result.warnings)
    recommended_plan = preprocessing_recommendation_result.details.get("plan", {})
    recommended_operations = recommended_plan.get("recommended_operations")
    if run_preprocessing_first:
        preprocess_result = run_preprocessing(
            imzml_path,
            operations=preprocessing_operations or recommended_operations,
            output_dir=str(output_root / "preprocess"),
        )
        warnings.extend(preprocess_result.warnings)
        if not preprocess_result.ok:
            return ToolResult.failure(
                tool_name,
                "End-to-end workflow failed during preprocessing stage.",
                imzml_path,
                errors=preprocess_result.errors,
                warnings=warnings,
            )
        analysis_input_path = preprocess_result.artifacts.get("processed_imzml", imzml_path)

    chemical_result: ToolResult | None = None
    chemical_context: dict[str, Any] | None = None
    target_mz_values = _collect_candidate_mz_values(
        qc_result.details.get("sampled_spectra", []),
        max_mz_values=chemical_max_mz_values,
    )
    if resolved_chemical_library_dir and target_mz_values:
        chemical_result = annotate_mz_with_local_library(
            target_mz_values,
            library_dir=resolved_chemical_library_dir,
            source_name=chemical_source_name,
            ion_mode=resolved_chemical_ion_mode,
            ppm_tolerance=chemical_ppm_tolerance,
            max_candidates_per_mz=chemical_max_candidates_per_mz,
            output_dir=str(output_root / "chemical"),
        )
        warnings.extend(chemical_result.warnings)
        if chemical_result.ok:
            chemical_context = chemical_result.details.get("chemical_context")
        else:
            warnings.extend(chemical_result.errors)
    elif resolved_chemical_library_dir is None:
        warnings.append("Chemical annotation skipped because no local chemical library directory was found.")
    elif not target_mz_values:
        warnings.append("Chemical annotation skipped because no representative m/z values were extracted from QC samples.")

    artifacts: dict[str, str] = {}
    artifacts.update(_collect_stage_artifacts("load", load_result))
    artifacts.update(_collect_stage_artifacts("qc", qc_result))
    if preprocess_result is not None:
        artifacts.update(_collect_stage_artifacts("preprocess", preprocess_result))
    if preprocessing_recommendation_result.ok:
        artifacts.update(_collect_stage_artifacts("preprocess_recommendation", preprocessing_recommendation_result))
    if chemical_result is not None:
        artifacts.update(_collect_stage_artifacts("chemical", chemical_result))

    unified_evidence = build_unified_evidence(
        biological_context=biological_context,
        chemical_context=chemical_context,
        statistical_evidence=statistical_evidence,
    )
    conclusion = build_integrated_conclusion(
        biological_context=biological_context,
        chemical_context=chemical_context,
        statistical_evidence=statistical_evidence,
    )
    metrics = {
        "load_spectra_count": load_result.metrics.get("spectra_count"),
        "qc_sampled_spectra_count": qc_result.metrics.get("sampled_spectra_count"),
        "biological_metadata_provided": load_result.metrics.get("biological_metadata_provided", False),
        "statistical_evidence_provided": bool(statistical_evidence and statistical_evidence.get("provided")),
        "chemical_annotation_enabled": chemical_result is not None,
        "chemical_context_provided": bool(chemical_context and chemical_context.get("provided")),
        "conclusion_ready": conclusion.get("ready_for_integrated_conclusion", False),
        "conclusion_confidence": conclusion.get("confidence"),
        "preprocessing_plan_available": preprocessing_recommendation_result.ok,
        "preprocessing_enabled": run_preprocessing_first,
    }

    result = ToolResult.success(
        tool_name,
        "End-to-end bootstrap workflow completed.",
        imzml_path,
        artifacts=artifacts,
        metrics=metrics,
        details={
            "load": load_result.to_dict(),
            "qc": qc_result.to_dict(),
            "preprocessing_recommendation": preprocessing_recommendation_result.to_dict(),
            "preprocess": None if preprocess_result is None else preprocess_result.to_dict(),
            "chemical": None if chemical_result is None else chemical_result.to_dict(),
            "final_analysis_input": analysis_input_path,
            "biological_context": biological_context,
            "statistical_evidence": statistical_evidence,
            "chemical_context": chemical_context,
            "chemical_query_mz_values": target_mz_values,
            "evidence": unified_evidence,
            "conclusion": conclusion,
        },
        warnings=_deduplicate(warnings),
    )

    artifact_path = write_result_json(result, str(output_root), f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def _collect_stage_artifacts(prefix: str, result: ToolResult) -> dict[str, str]:
    return {f"{prefix}_{key}": value for key, value in result.artifacts.items()}


def _deduplicate(items: list[str]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for item in items:
        if item not in seen:
            ordered.append(item)
            seen.add(item)
    return ordered


def _resolve_chemical_library_dir(explicit_dir: str | None) -> str | None:
    if explicit_dir:
        candidate = Path(explicit_dir)
        return str(candidate) if candidate.exists() and candidate.is_dir() else None

    default_dir = Path.cwd() / "data" / "chemical_db"
    if default_dir.exists() and default_dir.is_dir():
        return str(default_dir)
    return None


def _infer_ion_mode(biological_context: dict[str, Any] | None) -> str:
    if biological_context:
        experiment = biological_context.get("experiment") or {}
        ion_mode = experiment.get("ion_mode")
        if isinstance(ion_mode, str) and ion_mode.strip().lower() in {"positive", "negative"}:
            return ion_mode.strip().lower()
    return "positive"


def _collect_candidate_mz_values(
    sampled_spectra: list[dict[str, Any]],
    *,
    max_mz_values: int,
) -> list[float]:
    values: list[float] = []
    seen: set[float] = set()
    for row in sampled_spectra:
        for peak in row.get("top_peaks", []):
            mz_value = round(float(peak["mz"]), 6)
            if mz_value in seen:
                continue
            seen.add(mz_value)
            values.append(mz_value)
            if len(values) >= max_mz_values:
                return values
    return values
