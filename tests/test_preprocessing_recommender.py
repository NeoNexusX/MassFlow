from massflow.agent_tools.preprocessing_recommender import _build_recommendation_plan


def test_build_recommendation_plan_recommends_noise_and_baseline() -> None:
    plan = _build_recommendation_plan(
        dataset_summary={"occupancy": {"occupancy_ratio": 0.8}},
        qc_metrics={
            "tic_summary": {"cv": 0.3},
            "snr_summary": {"median": 2.8},
            "peak_count_summary": {"mean": 1400.0},
        },
    )

    operation_names = [item["name"] for item in plan["recommended_operations"]]
    assert "noise_reduction" in operation_names
    assert "baseline_correction" in operation_names
    assert "normalization" in operation_names
    assert plan["risk_level"] in {"medium", "high"}


def test_build_recommendation_plan_falls_back_to_default_chain() -> None:
    plan = _build_recommendation_plan(
        dataset_summary={"occupancy": {"occupancy_ratio": 0.8}},
        qc_metrics={
            "tic_summary": {"cv": 0.05},
            "snr_summary": {"median": 5.0},
            "peak_count_summary": {"mean": 200.0},
        },
    )

    operation_names = [item["name"] for item in plan["recommended_operations"]]
    assert operation_names == ["baseline_correction", "noise_reduction", "normalization"]
    assert plan["risk_level"] == "low"
