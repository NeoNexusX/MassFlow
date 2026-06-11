from __future__ import annotations

from typing import Any

from ._utils import write_result_json
from .preprocess_tools import DEFAULT_PREPROCESS_OPERATIONS
from .qc_tools import dataset_qc_summary
from .schemas import ToolResult


def recommend_preprocessing_plan(
    imzml_path: str,
    *,
    biological_metadata: dict[str, Any] | None = None,
    biological_metadata_json: str | None = None,
    biological_metadata_path: str | None = None,
    sample_size: int = 16,
    noise_method: str = "sd",
    max_threads: int = 8,
    output_dir: str | None = None,
) -> ToolResult:
    """Recommend preprocessing operations from QC evidence."""

    tool_name = "recommend_preprocessing_plan"
    qc_result = dataset_qc_summary(
        imzml_path,
        biological_metadata=biological_metadata,
        biological_metadata_json=biological_metadata_json,
        biological_metadata_path=biological_metadata_path,
        sample_size=sample_size,
        noise_method=noise_method,
        max_threads=max_threads,
        output_dir=output_dir,
    )
    if not qc_result.ok:
        return ToolResult.failure(
            tool_name,
            "Failed to recommend preprocessing plan because QC did not complete.",
            imzml_path,
            errors=qc_result.errors,
            warnings=qc_result.warnings,
        )

    dataset_summary = qc_result.details.get("dataset", {})
    qc_metrics = qc_result.metrics
    plan = _build_recommendation_plan(dataset_summary=dataset_summary, qc_metrics=qc_metrics)

    result = ToolResult.success(
        tool_name,
        "Recommended a preprocessing plan from dataset QC evidence.",
        imzml_path,
        metrics={
            "recommended_operation_count": len(plan["recommended_operations"]),
            "needs_preprocessing": bool(plan["recommended_operations"]),
            "risk_level": plan["risk_level"],
        },
        details={
            "qc": qc_result.to_dict(),
            "plan": plan,
        },
        warnings=list(qc_result.warnings),
    )
    artifact_path = write_result_json(result, output_dir, f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def _build_recommendation_plan(
    *,
    dataset_summary: dict[str, Any],
    qc_metrics: dict[str, Any],
) -> dict[str, Any]:
    recommended_operations: list[dict[str, Any]] = []
    reasons: list[str] = []
    expected_effects: list[str] = []
    risks: list[str] = []

    occupancy_ratio = (dataset_summary.get("occupancy") or {}).get("occupancy_ratio")
    tic_summary = qc_metrics.get("tic_summary") or {}
    snr_summary = qc_metrics.get("snr_summary") or {}
    peak_count_summary = qc_metrics.get("peak_count_summary") or {}

    snr_median = snr_summary.get("median")
    peak_count_mean = peak_count_summary.get("mean")
    tic_cv = tic_summary.get("cv")

    if _safe_lt(snr_median, 3.5):
        recommended_operations.append(
            {"name": "noise_reduction", "params": {"method": "gaussian_numba"}}
        )
        reasons.append("抽样谱 SNR 中位数偏低，建议先做去噪。")
        expected_effects.append("提升弱峰可分辨性，减少噪声对后续化学注释和结论的干扰。")

    if _safe_lt(snr_median, 4.0) or _safe_gt(peak_count_mean, 1200.0):
        recommended_operations.append(
            {"name": "baseline_correction", "params": {"method": "locmin_numba"}}
        )
        reasons.append("谱图存在基线漂移风险，建议做基线校正。")
        expected_effects.append("压低背景抬升，改善峰强度比较和后续归一化稳定性。")

    if _safe_gt(tic_cv, 0.25):
        recommended_operations.append(
            {"name": "normalization", "params": {"method": "tic_numba"}}
        )
        reasons.append("TIC 波动较大，建议做 TIC 归一化。")
        expected_effects.append("降低整体信号量级差异，提升不同像素/谱图之间的可比性。")

    if _safe_gt(peak_count_mean, 1500.0):
        recommended_operations.append({"name": "peak_pick", "params": {}})
        reasons.append("平均峰数较高，建议峰提取以压缩冗余峰并突出主要信号。")
        expected_effects.append("减少化学注释查询量，提高候选排序效率。")

    if occupancy_ratio is not None and occupancy_ratio < 0.1:
        risks.append("像素占用率较低，预处理前应先确认 ROI 或稀疏采样是否符合预期。")
    if not recommended_operations:
        recommended_operations = [dict(item) for item in DEFAULT_PREPROCESS_OPERATIONS]
        reasons.append("当前 QC 未显示极端异常，采用默认稳健预处理链。")
        expected_effects.append("获得基础一致的预处理输出，便于进入化学注释和综合结论阶段。")

    recommended_operations = _deduplicate_operations(recommended_operations)
    risk_level = _estimate_risk_level(
        snr_median=snr_median,
        tic_cv=tic_cv,
        occupancy_ratio=occupancy_ratio,
    )
    return {
        "recommended_operations": recommended_operations,
        "reasons": reasons,
        "expected_effects": expected_effects,
        "risks": risks,
        "risk_level": risk_level,
        "should_run_before_conclusion": True,
    }


def _deduplicate_operations(operations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    deduplicated: list[dict[str, Any]] = []
    seen: set[str] = set()
    for operation in operations:
        name = operation["name"]
        if name in seen:
            continue
        deduplicated.append(operation)
        seen.add(name)
    return deduplicated


def _estimate_risk_level(
    *,
    snr_median: Any,
    tic_cv: Any,
    occupancy_ratio: Any,
) -> str:
    score = 0
    if _safe_lt(snr_median, 3.0):
        score += 2
    elif _safe_lt(snr_median, 4.0):
        score += 1
    if _safe_gt(tic_cv, 0.35):
        score += 2
    elif _safe_gt(tic_cv, 0.2):
        score += 1
    if occupancy_ratio is not None and occupancy_ratio < 0.1:
        score += 2

    if score >= 4:
        return "high"
    if score >= 2:
        return "medium"
    return "low"


def _safe_lt(value: Any, threshold: float) -> bool:
    try:
        return float(value) < threshold
    except (TypeError, ValueError):
        return False


def _safe_gt(value: Any, threshold: float) -> bool:
    try:
        return float(value) > threshold
    except (TypeError, ValueError):
        return False
