from massflow.agent_tools.validation_tools import (
    build_comparison_summary,
    build_legacy_statistical_baseline,
    build_validation_checks,
)


def test_build_legacy_statistical_baseline_detects_stable_statistics() -> None:
    baseline = build_legacy_statistical_baseline(
        {
            "provided": True,
            "dataset_summary": {"occupancy": {"occupancy_ratio": 0.82}},
            "qc_summary": {
                "snr_summary": {"median": 4.2},
                "peak_count_summary": {"mean": 120.0},
            },
            "risk_flags": [],
        }
    )

    assert baseline["assessment"] == "stable_statistics"
    assert baseline["confidence"] == "medium"


def test_build_validation_checks_reports_presence_flags() -> None:
    checks = build_validation_checks(
        {
            "biological_context": {"provided": True},
            "chemical_context": {"provided": False},
            "statistical_evidence": {"provided": True},
            "conclusion": {"provided": True, "primary_basis": ["x"], "recommended_next_actions": ["y"]},
        }
    )

    assert checks["biological_context_present"] is True
    assert checks["chemical_context_present"] is False
    assert checks["conclusion_has_basis"] is True


def test_build_comparison_summary_mentions_primary_sources() -> None:
    summary = build_comparison_summary(
        integrated_conclusion={"confidence": "medium", "conclusion_summary": "新流程结论"},
        legacy_statistical_baseline={"confidence": "low", "summary": "旧流程摘要"},
        evidence={"primary_sources": ["chemical"], "auxiliary_sources": ["statistical"], "biological_context": {"provided": False}, "chemical_context": {"provided": True}},
    )

    assert any("chemical" in message for message in summary["messages"])
    assert any("旧流程摘要" in message for message in summary["messages"])
