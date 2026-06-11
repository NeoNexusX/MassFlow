from __future__ import annotations

from typing import Any


def build_integrated_conclusion(
    *,
    biological_context: dict[str, Any] | None = None,
    chemical_context: dict[str, Any] | None = None,
    statistical_evidence: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a conservative, evidence-first conclusion object."""

    biological = biological_context or {"provided": False}
    chemical = chemical_context or {"provided": False}
    statistical = statistical_evidence or {"provided": False, "risk_flags": []}

    valid_biological = _has_meaningful_biological_context(biological)
    valid_chemical_hits = _collect_valid_chemical_hits(chemical)
    risk_flags = list(statistical.get("risk_flags") or [])

    primary_basis: list[str] = []
    supporting_evidence: list[str] = []
    limitations: list[str] = []
    recommended_next_actions: list[str] = []

    if valid_biological:
        primary_basis.append(_summarize_biological_basis(biological))
    else:
        limitations.append("缺少可用的结构化生物信息，无法建立明确的生物学上下文。")
        recommended_next_actions.append("补充样本物种、组织、区域、条件和离子模式等结构化生物信息。")

    if valid_chemical_hits:
        primary_basis.append(_summarize_chemical_basis(valid_chemical_hits))
    else:
        limitations.append("当前本地化学注释未提供可靠候选，化学证据仍不足。")
        recommended_next_actions.append("扩充本地化学数据库，或调整 ppm 容差和离子加合物规则后重新注释。")

    if statistical.get("provided"):
        supporting_evidence.append(_summarize_statistical_support(statistical))
    else:
        limitations.append("未生成统计辅助证据，无法对信号稳定性和质量风险做补强判断。")

    supporting_evidence.extend(risk_flags[:3])
    confidence = _estimate_confidence(
        valid_biological=valid_biological,
        valid_chemical_hit_count=len(valid_chemical_hits),
        risk_flag_count=len(risk_flags),
    )

    if valid_biological and valid_chemical_hits:
        conclusion_summary = (
            "已形成以生物信息和化学注释为主、统计信息为辅的综合证据，可用于生成解释型结论。"
        )
    elif valid_chemical_hits:
        conclusion_summary = "当前结论主要由化学注释驱动，生物学解释仍需补充。"
    elif valid_biological:
        conclusion_summary = "当前结论主要由生物信息驱动，仍需化学注释补充分子层证据。"
    else:
        conclusion_summary = "当前证据不足以形成可靠综合结论，仅能输出数据级提示。"

    if risk_flags:
        recommended_next_actions.append("优先处理统计风险标记，再评估结论是否需要降级或人工复核。")
    if valid_chemical_hits and not valid_biological:
        recommended_next_actions.append("将主要注释到的分子与样本组织来源进行人工交叉核对。")
    if valid_biological and not valid_chemical_hits:
        recommended_next_actions.append("增加针对目标组织/代谢物类别的专用数据库，提高命中率。")

    return {
        "provided": True,
        "conclusion_summary": conclusion_summary,
        "confidence": confidence,
        "primary_basis": _deduplicate(primary_basis),
        "supporting_evidence": _deduplicate(supporting_evidence),
        "limitations": _deduplicate(limitations),
        "recommended_next_actions": _deduplicate(recommended_next_actions),
        "ready_for_integrated_conclusion": bool(valid_biological or valid_chemical_hits),
    }


def _has_meaningful_biological_context(biological_context: dict[str, Any]) -> bool:
    if not biological_context.get("provided"):
        return False
    sample = biological_context.get("sample") or {}
    meaningful_fields = ("species", "organ", "tissue_region", "condition", "sample_type")
    return any(_is_meaningful_text(sample.get(field)) for field in meaningful_fields)


def _collect_valid_chemical_hits(chemical_context: dict[str, Any]) -> list[dict[str, Any]]:
    if not chemical_context.get("provided"):
        return []
    top_hits = chemical_context.get("top_hits") or []
    valid_hits: list[dict[str, Any]] = []
    for hit in top_hits:
        formula = str(hit.get("formula") or "").strip().lower()
        if not formula or formula == "unannotated":
            continue
        valid_hits.append(hit)
    return valid_hits


def _summarize_biological_basis(biological_context: dict[str, Any]) -> str:
    sample = biological_context.get("sample") or {}
    experiment = biological_context.get("experiment") or {}
    species = _fallback_text(sample.get("species"), "unknown species")
    organ = _fallback_text(sample.get("organ"), "unknown organ")
    condition = _fallback_text(sample.get("condition"), "unknown condition")
    ion_mode = _fallback_text(experiment.get("ion_mode"), "unknown ion mode")
    return f"生物信息显示样本来自 {species} / {organ}，条件为 {condition}，离子模式为 {ion_mode}。"


def _summarize_chemical_basis(valid_chemical_hits: list[dict[str, Any]]) -> str:
    top_hit = valid_chemical_hits[0]
    formula = _fallback_text(top_hit.get("formula"), "unknown formula")
    compound_name = _fallback_text(top_hit.get("compound_name"), "未命名候选")
    source_name = _fallback_text(top_hit.get("source_name"), "unknown source")
    confidence = _fallback_text(top_hit.get("confidence_level"), "unknown confidence")
    return (
        f"化学注释已命中 {len(valid_chemical_hits)} 个代表性 m/z，首要候选为 {compound_name} "
        f"({formula})，来源 {source_name}，置信度 {confidence}。"
    )


def _summarize_statistical_support(statistical_evidence: dict[str, Any]) -> str:
    qc_summary = statistical_evidence.get("qc_summary") or {}
    sampled_count = qc_summary.get("sampled_spectra_count")
    snr_summary = qc_summary.get("snr_summary") or {}
    snr_median = snr_summary.get("median")
    occupancy = (statistical_evidence.get("dataset_summary") or {}).get("occupancy") or {}
    occupancy_ratio = occupancy.get("occupancy_ratio")
    return (
        "统计辅助证据覆盖 "
        f"{sampled_count} 条抽样谱，occupancy_ratio={_format_float(occupancy_ratio)}，"
        f"SNR median={_format_float(snr_median)}。"
    )


def _estimate_confidence(
    *,
    valid_biological: bool,
    valid_chemical_hit_count: int,
    risk_flag_count: int,
) -> str:
    score = 0
    if valid_biological:
        score += 2
    if valid_chemical_hit_count >= 5:
        score += 2
    elif valid_chemical_hit_count > 0:
        score += 1
    if risk_flag_count == 0:
        score += 1
    elif risk_flag_count >= 3:
        score -= 1

    if score >= 4:
        return "high"
    if score >= 2:
        return "medium"
    if score >= 1:
        return "low"
    return "very_low"


def _is_meaningful_text(value: Any) -> bool:
    text = str(value).strip().lower() if value is not None else ""
    return bool(text and text not in {"unknown", "pending", "n/a", "none", "null"})


def _fallback_text(value: Any, default: str) -> str:
    return str(value).strip() if _is_meaningful_text(value) else default


def _format_float(value: Any) -> str:
    if value is None:
        return "n/a"
    try:
        return f"{float(value):.4f}"
    except (TypeError, ValueError):
        return "n/a"


def _deduplicate(items: list[str]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for item in items:
        if item not in seen:
            ordered.append(item)
            seen.add(item)
    return ordered
