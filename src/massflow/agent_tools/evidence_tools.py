from __future__ import annotations

from typing import Any


def build_statistical_evidence(
    *,
    dataset_summary: dict[str, Any],
    qc_metrics: dict[str, Any],
    qc_flags: list[str],
    sampled_spectra: list[dict[str, Any]],
) -> dict[str, Any]:
    """Normalize QC statistics into an auxiliary evidence object."""

    return {
        "provided": True,
        "role": "auxiliary",
        "dataset_summary": {
            "coordinates": dataset_summary.get("coordinates"),
            "occupancy": dataset_summary.get("occupancy"),
            "shared_mz_length": dataset_summary.get("shared_mz_length"),
        },
        "qc_summary": {
            "sampled_spectra_count": qc_metrics.get("sampled_spectra_count"),
            "tic_summary": qc_metrics.get("tic_summary"),
            "base_peak_summary": qc_metrics.get("base_peak_summary"),
            "peak_count_summary": qc_metrics.get("peak_count_summary"),
            "snr_summary": qc_metrics.get("snr_summary"),
        },
        "risk_flags": list(qc_flags),
        "sample_preview": sampled_spectra[:5],
    }


def build_unified_evidence(
    *,
    biological_context: dict[str, Any] | None = None,
    chemical_context: dict[str, Any] | None = None,
    statistical_evidence: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one evidence object for downstream agent reasoning."""

    biological = biological_context or _empty_biological_context()
    chemical = chemical_context or _empty_chemical_context()
    statistical = statistical_evidence or _empty_statistical_evidence()
    primary_sources = []
    if biological.get("provided"):
        primary_sources.append("biological")
    if chemical.get("provided"):
        primary_sources.append("chemical")
    auxiliary_sources = []
    if statistical.get("provided"):
        auxiliary_sources.append("statistical")

    return {
        "primary_sources": primary_sources,
        "auxiliary_sources": auxiliary_sources,
        "biological_context": biological,
        "chemical_context": chemical,
        "statistical_evidence": statistical,
        "ready_for_integrated_conclusion": bool(primary_sources),
    }


def _empty_biological_context() -> dict[str, Any]:
    return {
        "provided": False,
        "source": None,
        "schema_version": None,
        "sample": None,
        "experiment": None,
        "binding": None,
        "additional_fields": {},
    }


def _empty_chemical_context() -> dict[str, Any]:
    return {
        "provided": False,
        "schema_version": None,
        "annotation_count": 0,
        "candidate_count": 0,
        "top_hits": [],
    }


def _empty_statistical_evidence() -> dict[str, Any]:
    return {
        "provided": False,
        "role": "auxiliary",
        "dataset_summary": None,
        "qc_summary": None,
        "risk_flags": [],
        "sample_preview": [],
    }
