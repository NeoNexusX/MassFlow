from __future__ import annotations

from typing import Any

from ._utils import write_result_json
from .schemas import ToolResult
from .workflow_tools import analyze_imzml_end_to_end


def validate_integrated_conclusion(
    imzml_path: str,
    *,
    biological_metadata: dict[str, Any] | None = None,
    biological_metadata_json: str | None = None,
    biological_metadata_path: str | None = None,
    chemical_library_dir: str | None = None,
    chemical_source_name: str = "LIPID MAPS",
    chemical_ion_mode: str | None = None,
    chemical_ppm_tolerance: float = 5.0,
    chemical_max_mz_values: int = 12,
    chemical_max_candidates_per_mz: int = 5,
    run_preprocessing_first: bool = False,
    preprocessing_operations: list[dict[str, Any]] | None = None,
    output_dir: str | None = None,
) -> ToolResult:
    """Run the integrated workflow and compare it with a legacy statistical-only baseline."""

    tool_name = "validate_integrated_conclusion"
    workflow_result = analyze_imzml_end_to_end(
        imzml_path,
        run_preprocessing_first=run_preprocessing_first,
        preprocessing_operations=preprocessing_operations,
        biological_metadata=biological_metadata,
        biological_metadata_json=biological_metadata_json,
        biological_metadata_path=biological_metadata_path,
        chemical_library_dir=chemical_library_dir,
        chemical_source_name=chemical_source_name,
        chemical_ion_mode=chemical_ion_mode,
        chemical_ppm_tolerance=chemical_ppm_tolerance,
        chemical_max_mz_values=chemical_max_mz_values,
        chemical_max_candidates_per_mz=chemical_max_candidates_per_mz,
        output_dir=output_dir,
    )
    if not workflow_result.ok:
        return ToolResult.failure(
            tool_name,
            "Failed to validate integrated conclusion because the workflow did not complete.",
            imzml_path,
            errors=workflow_result.errors,
            warnings=workflow_result.warnings,
        )

    integrated_conclusion = workflow_result.details.get("conclusion") or {}
    statistical_evidence = workflow_result.details.get("statistical_evidence") or {}
    evidence = workflow_result.details.get("evidence") or {}
    legacy_statistical_baseline = build_legacy_statistical_baseline(statistical_evidence)
    validation_checks = build_validation_checks(workflow_result.details)
    comparison_summary = build_comparison_summary(
        integrated_conclusion=integrated_conclusion,
        legacy_statistical_baseline=legacy_statistical_baseline,
        evidence=evidence,
    )

    result = ToolResult.success(
        tool_name,
        "Validated integrated conclusion against a legacy statistical-only baseline.",
        imzml_path,
        artifacts=dict(workflow_result.artifacts),
        metrics={
            "conclusion_ready": workflow_result.metrics.get("conclusion_ready"),
            "integrated_confidence": integrated_conclusion.get("confidence"),
            "legacy_confidence": legacy_statistical_baseline.get("confidence"),
            "check_pass_count": sum(1 for item in validation_checks.values() if item is True),
            "check_count": len(validation_checks),
        },
        details={
            "workflow": workflow_result.to_dict(),
            "integrated_conclusion": integrated_conclusion,
            "legacy_statistical_baseline": legacy_statistical_baseline,
            "validation_checks": validation_checks,
            "comparison_summary": comparison_summary,
        },
        warnings=list(workflow_result.warnings),
    )
    artifact_path = write_result_json(result, output_dir, f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["validation_json"] = artifact_path
    return result


def build_legacy_statistical_baseline(statistical_evidence: dict[str, Any] | None) -> dict[str, Any]:
    """Approximate the older statistics-first behavior for comparison."""

    evidence = statistical_evidence or {}
    if not evidence.get("provided"):
        return {
            "provided": False,
            "summary": "旧统计基线不可用，因为没有统计辅助证据。",
            "confidence": "very_low",
            "main_flags": [],
            "assessment": "insufficient_statistics",
        }

    qc_summary = evidence.get("qc_summary") or {}
    risk_flags = list(evidence.get("risk_flags") or [])
    occupancy_ratio = ((evidence.get("dataset_summary") or {}).get("occupancy") or {}).get("occupancy_ratio")
    snr_median = (qc_summary.get("snr_summary") or {}).get("median")
    peak_count_mean = (qc_summary.get("peak_count_summary") or {}).get("mean")

    if risk_flags:
        assessment = "warning_detected"
        confidence = "low" if len(risk_flags) >= 2 else "medium"
        summary = "旧统计基线认为当前数据存在质量风险，结论需谨慎解释。"
    elif _safe_ge(occupancy_ratio, 0.2) and _safe_ge(snr_median, 3.0) and _safe_ge(peak_count_mean, 50.0):
        assessment = "stable_statistics"
        confidence = "medium"
        summary = "旧统计基线认为当前数据统计表现稳定，可支持初步数据级判断。"
    else:
        assessment = "limited_statistics"
        confidence = "low"
        summary = "旧统计基线认为统计信号有限，只适合输出弱结论。"

    return {
        "provided": True,
        "summary": summary,
        "confidence": confidence,
        "main_flags": risk_flags[:3],
        "assessment": assessment,
        "occupancy_ratio": occupancy_ratio,
        "snr_median": snr_median,
        "peak_count_mean": peak_count_mean,
    }


def build_validation_checks(workflow_details: dict[str, Any]) -> dict[str, bool]:
    biological_context = workflow_details.get("biological_context") or {}
    chemical_context = workflow_details.get("chemical_context") or {}
    statistical_evidence = workflow_details.get("statistical_evidence") or {}
    conclusion = workflow_details.get("conclusion") or {}
    return {
        "statistical_evidence_present": bool(statistical_evidence.get("provided")),
        "chemical_context_present": bool(chemical_context.get("provided")),
        "biological_context_present": bool(biological_context.get("provided")),
        "conclusion_present": bool(conclusion.get("provided")),
        "conclusion_has_basis": bool(conclusion.get("primary_basis")),
        "conclusion_has_actions": bool(conclusion.get("recommended_next_actions")),
    }


def build_comparison_summary(
    *,
    integrated_conclusion: dict[str, Any],
    legacy_statistical_baseline: dict[str, Any],
    evidence: dict[str, Any],
) -> dict[str, Any]:
    messages: list[str] = []
    primary_sources = list(evidence.get("primary_sources") or [])
    auxiliary_sources = list(evidence.get("auxiliary_sources") or [])

    if primary_sources:
        messages.append(f"新流程的主证据来源为: {', '.join(primary_sources)}。")
    if auxiliary_sources:
        messages.append(f"辅助证据来源为: {', '.join(auxiliary_sources)}。")

    integrated_confidence = integrated_conclusion.get("confidence")
    legacy_confidence = legacy_statistical_baseline.get("confidence")
    if integrated_confidence != legacy_confidence:
        messages.append(
            f"新流程置信度为 {integrated_confidence}，旧统计基线为 {legacy_confidence}。"
        )
    if integrated_conclusion.get("conclusion_summary"):
        messages.append(integrated_conclusion["conclusion_summary"])
    if legacy_statistical_baseline.get("summary"):
        messages.append(legacy_statistical_baseline["summary"])

    return {
        "messages": messages,
        "recommended_validation_scenarios": _build_recommended_validation_scenarios(evidence),
    }


def _build_recommended_validation_scenarios(evidence: dict[str, Any]) -> list[str]:
    scenarios = [
        "对同一个 imzML 分别运行：仅统计基线、化学注释开启、化学注释+真实生物信息三种模式。",
        "检查 conclusion_summary、confidence、limitations 是否随生物信息补全而收敛。",
    ]
    if not (evidence.get("biological_context") or {}).get("provided"):
        scenarios.append("补充真实 biological metadata 后再次运行，验证主证据是否从 chemical-only 升级为 biological+chemical。")
    if not (evidence.get("chemical_context") or {}).get("provided"):
        scenarios.append("调整本地化学库目录、ppm 容差和离子模式后重试，检查 chemical_context 是否成功生成。")
    return scenarios


def _safe_ge(value: Any, threshold: float) -> bool:
    try:
        return float(value) >= threshold
    except (TypeError, ValueError):
        return False
