from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .chemical_db_loader import ChemicalLibraryEntry
from .chemical_schema import ChemicalAnnotationRecord


PROTON_MASS = 1.007276466812
SODIUM_MASS = 22.989218
POTASSIUM_MASS = 38.963158
AMMONIUM_MASS = 18.033823
WATER_MASS = 18.010565


@dataclass(frozen=True, slots=True)
class AdductRule:
    name: str
    ion_mode: str
    mass_delta: float
    charge: int


COMMON_ADDUCT_RULES: tuple[AdductRule, ...] = (
    AdductRule(name="[M+H]+", ion_mode="positive", mass_delta=PROTON_MASS, charge=1),
    AdductRule(name="[M+Na]+", ion_mode="positive", mass_delta=SODIUM_MASS, charge=1),
    AdductRule(name="[M+K]+", ion_mode="positive", mass_delta=POTASSIUM_MASS, charge=1),
    AdductRule(name="[M+NH4]+", ion_mode="positive", mass_delta=AMMONIUM_MASS, charge=1),
    AdductRule(name="[M-H2O+H]+", ion_mode="positive", mass_delta=PROTON_MASS - WATER_MASS, charge=1),
    AdductRule(name="[M-H]-", ion_mode="negative", mass_delta=-PROTON_MASS, charge=-1),
    AdductRule(name="[M+Cl]-", ion_mode="negative", mass_delta=34.968853, charge=-1),
    AdductRule(name="[M+FA-H]-", ion_mode="negative", mass_delta=44.997654, charge=-1),
)


def annotate_mz_values(
    mz_values: Iterable[float],
    *,
    library_entries: list[ChemicalLibraryEntry],
    ion_mode: str = "positive",
    ppm_tolerance: float = 5.0,
    max_candidates_per_mz: int = 5,
    query_label_prefix: str = "mz",
) -> list[ChemicalAnnotationRecord]:
    """Annotate one or more observed m/z values using a local library."""

    normalized_ion_mode = _normalize_ion_mode(ion_mode)
    relevant_adducts = [rule for rule in COMMON_ADDUCT_RULES if rule.ion_mode == normalized_ion_mode]
    annotations: list[ChemicalAnnotationRecord] = []
    for index, mz in enumerate(mz_values):
        annotations.append(
            annotate_single_mz(
                float(mz),
                library_entries=library_entries,
                adduct_rules=relevant_adducts,
                ion_mode=normalized_ion_mode,
                ppm_tolerance=ppm_tolerance,
                max_candidates=max_candidates_per_mz,
                query_label=f"{query_label_prefix}_{index + 1}",
            )
        )
    return annotations


def annotate_single_mz(
    mz: float,
    *,
    library_entries: list[ChemicalLibraryEntry],
    adduct_rules: list[AdductRule] | tuple[AdductRule, ...] = COMMON_ADDUCT_RULES,
    ion_mode: str = "positive",
    ppm_tolerance: float = 5.0,
    max_candidates: int = 5,
    query_label: str | None = None,
) -> ChemicalAnnotationRecord:
    """Annotate a single m/z value against a normalized local chemical library."""

    normalized_ion_mode = _normalize_ion_mode(ion_mode)
    matches: list[dict[str, object]] = []
    for entry in library_entries:
        for adduct_rule in adduct_rules:
            if adduct_rule.ion_mode != normalized_ion_mode:
                continue
            theoretical_mz = _compute_theoretical_mz(entry.exact_mass, adduct_rule)
            ppm_error = _compute_ppm_error(observed_mz=mz, theoretical_mz=theoretical_mz)
            if abs(ppm_error) > ppm_tolerance:
                continue
            matches.append(
                {
                    "formula": entry.formula,
                    "smiles": entry.smiles,
                    "compound_name": entry.compound_name,
                    "adduct": adduct_rule.name,
                    "charge": adduct_rule.charge,
                    "ppm_error": ppm_error,
                    "score": _score_match(ppm_error=ppm_error, has_smiles=bool(entry.smiles)),
                    "confidence_level": _confidence_from_ppm(ppm_error),
                    "source_type": "public_database",
                    "source_name": entry.source_name,
                    "source_record_id": entry.source_record_id,
                    "notes": None,
                    "class_label": entry.class_label,
                    "theoretical_mz": theoretical_mz,
                }
            )

    matches.sort(key=lambda item: (abs(float(item["ppm_error"])), -float(item["score"] or 0.0)))
    limited_matches = matches[:max_candidates]
    return ChemicalAnnotationRecord.from_mapping(
        {
            "mz": mz,
            "ppm_tolerance": ppm_tolerance,
            "ion_mode": normalized_ion_mode,
            "rank_strategy": "ppm_error_then_score",
            "query_label": query_label,
            "candidates": limited_matches or [_build_unannotated_candidate()],
        }
    )


def _compute_theoretical_mz(exact_mass: float, adduct_rule: AdductRule) -> float:
    return (exact_mass + adduct_rule.mass_delta) / abs(adduct_rule.charge)


def _compute_ppm_error(*, observed_mz: float, theoretical_mz: float) -> float:
    return ((observed_mz - theoretical_mz) / theoretical_mz) * 1_000_000.0


def _score_match(*, ppm_error: float, has_smiles: bool) -> float:
    base_score = max(0.0, 1.0 - (abs(ppm_error) / 10.0))
    if has_smiles:
        base_score += 0.05
    return round(base_score, 6)


def _confidence_from_ppm(ppm_error: float) -> str:
    error = abs(ppm_error)
    if error <= 2.0:
        return "high"
    if error <= 5.0:
        return "medium"
    return "low"


def _normalize_ion_mode(ion_mode: str) -> str:
    normalized = str(ion_mode).strip().lower()
    if normalized in {"positive", "pos", "+"}:
        return "positive"
    if normalized in {"negative", "neg", "-"}:
        return "negative"
    raise ValueError(f"Unsupported ion_mode: {ion_mode!r}")


def _build_unannotated_candidate() -> dict[str, object]:
    return {
        "formula": "unannotated",
        "smiles": None,
        "compound_name": None,
        "adduct": None,
        "charge": None,
        "ppm_error": None,
        "score": 0.0,
        "confidence_level": "unknown",
        "source_type": "unknown",
        "source_name": None,
        "source_record_id": None,
        "notes": "No candidate matched the local library under the configured ppm tolerance.",
    }
