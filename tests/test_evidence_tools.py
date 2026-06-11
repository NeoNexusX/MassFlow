from massflow.agent_tools.evidence_tools import build_statistical_evidence, build_unified_evidence


def test_build_statistical_evidence_marks_auxiliary_role() -> None:
    evidence = build_statistical_evidence(
        dataset_summary={
            "coordinates": {"spectra_count": 12},
            "occupancy": {"occupancy_ratio": 0.8},
            "shared_mz_length": 1200,
        },
        qc_metrics={
            "sampled_spectra_count": 4,
            "tic_summary": {"median": 100.0},
            "base_peak_summary": {"median": 25.0},
            "peak_count_summary": {"median": 80.0},
            "snr_summary": {"median": 5.0},
        },
        qc_flags=["Sampled median SNR is low; denoising or baseline correction may be needed."],
        sampled_spectra=[{"index": 1}, {"index": 2}, {"index": 3}, {"index": 4}, {"index": 5}, {"index": 6}],
    )

    assert evidence["provided"] is True
    assert evidence["role"] == "auxiliary"
    assert len(evidence["sample_preview"]) == 5


def test_build_unified_evidence_tracks_primary_and_auxiliary_sources() -> None:
    evidence = build_unified_evidence(
        biological_context={"provided": True, "sample": {"organ": "liver"}},
        chemical_context={"provided": False, "top_hits": []},
        statistical_evidence={"provided": True, "role": "auxiliary", "risk_flags": []},
    )

    assert evidence["primary_sources"] == ["biological"]
    assert evidence["auxiliary_sources"] == ["statistical"]
    assert evidence["ready_for_integrated_conclusion"] is True
