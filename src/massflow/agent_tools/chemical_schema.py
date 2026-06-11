from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Mapping


CHEMICAL_ANNOTATION_SCHEMA_VERSION = "1.0"

CHEMICAL_CONFIDENCE_LEVELS: tuple[str, ...] = (
    "high",
    "medium",
    "low",
    "unknown",
)

CHEMICAL_SOURCE_TYPES: tuple[str, ...] = (
    "public_database",
    "local_library",
    "rule_based",
    "manual",
    "unknown",
)


class ChemicalAnnotationValidationError(ValueError):
    """Raised when a chemical annotation payload is invalid."""


@dataclass(slots=True)
class ChemicalAnnotationCandidate:
    """One candidate molecular annotation for a target m/z."""

    formula: str
    smiles: str | None = None
    compound_name: str | None = None
    adduct: str | None = None
    charge: int | None = None
    ppm_error: float | None = None
    score: float | None = None
    confidence_level: str = "unknown"
    source_type: str = "unknown"
    source_name: str | None = None
    source_record_id: str | None = None
    notes: str | None = None
    additional_fields: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class ChemicalAnnotationRecord:
    """All candidate annotations for one observed m/z target."""

    mz: float
    candidates: list[ChemicalAnnotationCandidate]
    intensity: float | None = None
    rank_strategy: str | None = None
    ppm_tolerance: float | None = None
    ion_mode: str | None = None
    query_label: str | None = None
    schema_version: str = CHEMICAL_ANNOTATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> "ChemicalAnnotationRecord":
        if not isinstance(payload, Mapping):
            raise ChemicalAnnotationValidationError("Chemical annotation record must be a mapping.")

        mz = _coerce_float(payload.get("mz"), field_name="mz", required=True)
        candidates_payload = payload.get("candidates")
        if not isinstance(candidates_payload, list) or not candidates_payload:
            raise ChemicalAnnotationValidationError(
                "Chemical annotation record requires a non-empty candidates list."
            )

        candidates = [build_chemical_annotation_candidate(item) for item in candidates_payload]
        return cls(
            mz=mz,
            candidates=candidates,
            intensity=_coerce_float(payload.get("intensity"), field_name="intensity"),
            rank_strategy=_coerce_optional_str(payload.get("rank_strategy")),
            ppm_tolerance=_coerce_float(payload.get("ppm_tolerance"), field_name="ppm_tolerance"),
            ion_mode=_coerce_optional_str(payload.get("ion_mode")),
            query_label=_coerce_optional_str(payload.get("query_label")),
        )


def build_chemical_annotation_candidate(payload: Mapping[str, Any]) -> ChemicalAnnotationCandidate:
    """Validate and normalize one candidate annotation payload."""

    if not isinstance(payload, Mapping):
        raise ChemicalAnnotationValidationError("Chemical annotation candidate must be a mapping.")

    formula = _coerce_optional_str(payload.get("formula"))
    if not formula:
        raise ChemicalAnnotationValidationError("Chemical annotation candidate requires formula.")

    confidence_level = _normalize_enum(
        payload.get("confidence_level"),
        allowed=CHEMICAL_CONFIDENCE_LEVELS,
        field_name="confidence_level",
        default="unknown",
    )
    source_type = _normalize_enum(
        payload.get("source_type"),
        allowed=CHEMICAL_SOURCE_TYPES,
        field_name="source_type",
        default="unknown",
    )
    additional_fields = _extract_additional_fields(payload)
    return ChemicalAnnotationCandidate(
        formula=formula,
        smiles=_coerce_optional_str(payload.get("smiles")),
        compound_name=_coerce_optional_str(payload.get("compound_name")),
        adduct=_coerce_optional_str(payload.get("adduct")),
        charge=_coerce_int(payload.get("charge"), field_name="charge"),
        ppm_error=_coerce_float(payload.get("ppm_error"), field_name="ppm_error"),
        score=_coerce_float(payload.get("score"), field_name="score"),
        confidence_level=confidence_level,
        source_type=source_type,
        source_name=_coerce_optional_str(payload.get("source_name")),
        source_record_id=_coerce_optional_str(payload.get("source_record_id")),
        notes=_coerce_optional_str(payload.get("notes")),
        additional_fields=additional_fields,
    )


def build_chemical_annotation_template() -> dict[str, Any]:
    """Return a JSON-friendly template for chemical annotation outputs."""

    return {
        "schema_version": CHEMICAL_ANNOTATION_SCHEMA_VERSION,
        "mz": 0.0,
        "intensity": None,
        "rank_strategy": "ppm_error_then_score",
        "ppm_tolerance": 5.0,
        "ion_mode": "positive",
        "query_label": "",
        "candidates": [
            {
                "formula": "",
                "smiles": "",
                "compound_name": "",
                "adduct": "",
                "charge": 1,
                "ppm_error": 0.0,
                "score": 0.0,
                "confidence_level": "unknown",
                "source_type": "public_database",
                "source_name": "",
                "source_record_id": "",
                "notes": "",
            }
        ],
    }


def build_chemical_annotation_context(
    annotations: list[ChemicalAnnotationRecord] | None,
) -> dict[str, Any]:
    """Build an agent-friendly context summary from chemical annotations."""

    if not annotations:
        return {
            "provided": False,
            "schema_version": CHEMICAL_ANNOTATION_SCHEMA_VERSION,
            "annotation_count": 0,
            "candidate_count": 0,
            "top_hits": [],
        }

    top_hits: list[dict[str, Any]] = []
    candidate_count = 0
    for record in annotations:
        candidate_count += len(record.candidates)
        top_candidate = record.candidates[0]
        top_hits.append(
            {
                "mz": record.mz,
                "formula": top_candidate.formula,
                "smiles": top_candidate.smiles,
                "compound_name": top_candidate.compound_name,
                "confidence_level": top_candidate.confidence_level,
                "source_name": top_candidate.source_name,
            }
        )

    return {
        "provided": True,
        "schema_version": CHEMICAL_ANNOTATION_SCHEMA_VERSION,
        "annotation_count": len(annotations),
        "candidate_count": candidate_count,
        "top_hits": top_hits,
    }


def _coerce_optional_str(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _coerce_float(value: Any, *, field_name: str, required: bool = False) -> float | None:
    if value is None or value == "":
        if required:
            raise ChemicalAnnotationValidationError(f"{field_name} is required.")
        return None
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ChemicalAnnotationValidationError(f"{field_name} must be numeric.") from exc


def _coerce_int(value: Any, *, field_name: str) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ChemicalAnnotationValidationError(f"{field_name} must be an integer.") from exc


def _normalize_enum(
    value: Any,
    *,
    allowed: tuple[str, ...],
    field_name: str,
    default: str,
) -> str:
    normalized = _coerce_optional_str(value)
    if normalized is None:
        return default
    lowered = normalized.casefold()
    for item in allowed:
        if lowered == item.casefold():
            return item
    raise ChemicalAnnotationValidationError(
        f"{field_name} must be one of {', '.join(allowed)}, got {normalized!r}."
    )


def _extract_additional_fields(payload: Mapping[str, Any]) -> dict[str, Any]:
    known_fields = {
        "formula",
        "smiles",
        "compound_name",
        "adduct",
        "charge",
        "ppm_error",
        "score",
        "confidence_level",
        "source_type",
        "source_name",
        "source_record_id",
        "notes",
    }
    return {str(key): value for key, value in payload.items() if str(key) not in known_fields}
