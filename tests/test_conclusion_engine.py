from massflow.agent_tools.conclusion_engine import build_integrated_conclusion


def test_build_integrated_conclusion_with_biological_and_chemical_evidence() -> None:
    conclusion = build_integrated_conclusion(
        biological_context={
            "provided": True,
            "sample": {
                "species": "Homo sapiens",
                "organ": "liver",
                "condition": "tumor",
                "tissue_region": "segment-4",
                "sample_type": "tissue",
            },
            "experiment": {"ion_mode": "positive"},
        },
        chemical_context={
            "provided": True,
            "top_hits": [
                {
                    "mz": 760.58,
                    "formula": "C42H82NO8P",
                    "compound_name": "PC(34:1)",
                    "confidence_level": "high",
                    "source_name": "LIPID MAPS",
                },
                {
                    "mz": 782.56,
                    "formula": "C44H84NO8P",
                    "compound_name": "PC(36:2)",
                    "confidence_level": "medium",
                    "source_name": "LIPID MAPS",
                },
            ],
        },
        statistical_evidence={
            "provided": True,
            "dataset_summary": {"occupancy": {"occupancy_ratio": 0.82}},
            "qc_summary": {"sampled_spectra_count": 8, "snr_summary": {"median": 4.8}},
            "risk_flags": [],
        },
    )

    assert conclusion["ready_for_integrated_conclusion"] is True
    assert conclusion["confidence"] in {"medium", "high"}
    assert any("生物信息显示样本来自" in item for item in conclusion["primary_basis"])
    assert any("化学注释已命中" in item for item in conclusion["primary_basis"])


def test_build_integrated_conclusion_reports_missing_evidence() -> None:
    conclusion = build_integrated_conclusion(
        biological_context={"provided": False},
        chemical_context={"provided": True, "top_hits": [{"formula": "unannotated"}]},
        statistical_evidence={"provided": False, "risk_flags": []},
    )

    assert conclusion["ready_for_integrated_conclusion"] is False
    assert conclusion["confidence"] in {"very_low", "low"}
    assert any("缺少可用的结构化生物信息" in item for item in conclusion["limitations"])
    assert any("当前本地化学注释未提供可靠候选" in item for item in conclusion["limitations"])
