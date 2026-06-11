from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping


BIOLOGICAL_METADATA_SCHEMA_VERSION = "1.1"

REQUIRED_TOP_LEVEL_SECTIONS: tuple[str, ...] = ("sample", "experiment", "acquisition", "preparation", "binding")
REQUIRED_BIOLOGICAL_FIELDS: tuple[str, ...] = ("sample.sample_id", "sample.organism", "sample.organism_part")

CONTROLLED_VOCABULARY: dict[str, dict[str, str]] = {
    "ion_mode": {
        "positive": "positive",
        "pos": "positive",
        "+": "positive",
        "negative": "negative",
        "neg": "negative",
        "-": "negative",
        "unknown": "unknown",
    },
    "polarity": {
        "positive": "positive",
        "pos": "positive",
        "+": "positive",
        "negative": "negative",
        "neg": "negative",
        "-": "negative",
        "unknown": "unknown",
    },
    "sample_type": {
        "tissue": "tissue",
        "ffpe": "ffpe",
        "fresh_frozen": "fresh_frozen",
        "fresh-frozen": "fresh_frozen",
        "fresh frozen": "fresh_frozen",
        "serum": "serum",
        "plasma": "plasma",
        "cell_line": "cell_line",
        "cell line": "cell_line",
    },
}

BIOLOGICAL_METADATA_DEFAULTS: dict[str, Any] = {
    "tissue_region": "unknown",
    "sample_type": "tissue",
    "condition": "unknown",
    "group_label": "unknown",
    "ion_mode": "unknown",
    "polarity": "unknown",
    "prep_method": "unknown",
    "ionisation_source": "unknown",
    "analyzer": "unknown",
    "spectrum_mode": "unknown",
    "storage_mode": "unknown",
    "sample_stabilization": "unknown",
    "solvent": "unknown",
    "acquisition_mode": "MSI",
}

NUMERIC_FIELDS: frozenset[str] = frozenset({"pixel_size_x_um", "pixel_size_y_um"})
KNOWN_SECTION_FIELDS: dict[str, frozenset[str]] = {
    "sample": frozenset(
        {"sample_id", "organism", "organism_part", "tissue_region", "sample_type", "condition", "group_label"}
    ),
    "experiment": frozenset(
        {"polarity", "ionisation_source", "analyzer", "instrument", "spectrum_mode", "storage_mode", "acquisition_mode"}
    ),
    "acquisition": frozenset({"pixel_size_x_um", "pixel_size_y_um"}),
    "preparation": frozenset({"sample_stabilization", "solvent", "prep_method"}),
    "notes": frozenset({"submitter_notes"}),
    "binding": frozenset({"imzml_file_name", "imzml_basename"}),
}


class BiologicalMetadataValidationError(ValueError):
    """Raised when biological metadata fails validation."""


@dataclass(slots=True)
class BiologicalFileBinding:
    """Binding information that ties metadata to one imzML dataset."""

    imzml_file_name: str
    matched_by: str = "file_name"


@dataclass(slots=True)
class BiologicalMetadata:
    """Canonical schema for user-supplied biological context."""

    sample_id: str
    species: str
    organism: str
    organ: str
    organism_part: str
    tissue_region: str
    sample_type: str
    condition: str
    group_label: str
    ion_mode: str
    polarity: str
    prep_method: str
    submitter_notes: str | None = None
    instrument: str | None = None
    acquisition_mode: str | None = None
    ionisation_source: str | None = None
    analyzer: str | None = None
    spectrum_mode: str | None = None
    storage_mode: str | None = None
    pixel_size_x_um: float | None = None
    pixel_size_y_um: float | None = None
    sample_stabilization: str | None = None
    solvent: str | None = None
    schema_version: str = BIOLOGICAL_METADATA_SCHEMA_VERSION
    binding: BiologicalFileBinding | None = None
    additional_fields: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_mapping(
        cls,
        payload: Mapping[str, Any],
        *,
        imzml_path: str | None = None,
    ) -> "BiologicalMetadata":
        canonical = normalize_biological_metadata_payload(payload)
        binding = build_biological_file_binding(canonical, imzml_path=imzml_path)
        additional_fields = _extract_additional_fields(payload)
        return cls(
            sample_id=canonical["sample_id"],
            species=canonical["organism"],
            organism=canonical["organism"],
            organ=canonical["organism_part"],
            organism_part=canonical["organism_part"],
            tissue_region=canonical["tissue_region"],
            sample_type=canonical["sample_type"],
            condition=canonical["condition"],
            group_label=canonical["group_label"],
            ion_mode=canonical["polarity"],
            polarity=canonical["polarity"],
            prep_method=canonical["prep_method"],
            submitter_notes=canonical.get("submitter_notes"),
            instrument=canonical.get("instrument"),
            acquisition_mode=canonical.get("acquisition_mode"),
            ionisation_source=canonical.get("ionisation_source"),
            analyzer=canonical.get("analyzer"),
            spectrum_mode=canonical.get("spectrum_mode"),
            storage_mode=canonical.get("storage_mode"),
            pixel_size_x_um=canonical.get("pixel_size_x_um"),
            pixel_size_y_um=canonical.get("pixel_size_y_um"),
            sample_stabilization=canonical.get("sample_stabilization"),
            solvent=canonical.get("solvent"),
            binding=binding,
            additional_fields=additional_fields,
        )


def normalize_biological_metadata_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize aliases and controlled vocabulary, then validate required fields."""

    if not isinstance(payload, Mapping):
        raise BiologicalMetadataValidationError("Biological metadata payload must be a mapping.")

    missing_sections = [section_name for section_name in REQUIRED_TOP_LEVEL_SECTIONS if not isinstance(payload.get(section_name), Mapping)]
    if missing_sections:
        missing_text = ", ".join(missing_sections)
        raise BiologicalMetadataValidationError(
            f"Biological metadata must use the grouped JSON schema. Missing required sections: {missing_text}."
        )

    sample = _require_mapping_section(payload, "sample")
    experiment = _require_mapping_section(payload, "experiment")
    acquisition = _require_mapping_section(payload, "acquisition")
    preparation = _require_mapping_section(payload, "preparation")
    notes = _optional_mapping_section(payload, "notes")

    canonical: dict[str, Any] = {
        "sample_id": _normalize_scalar_value(sample.get("sample_id")),
        "organism": _normalize_scalar_value(sample.get("organism")),
        "organism_part": _normalize_scalar_value(sample.get("organism_part")),
        "tissue_region": _normalize_scalar_value(sample.get("tissue_region")),
        "sample_type": _normalize_controlled_value("sample_type", _normalize_scalar_value(sample.get("sample_type")) or "tissue"),
        "condition": _normalize_scalar_value(sample.get("condition")),
        "group_label": _normalize_scalar_value(sample.get("group_label")),
        "polarity": _normalize_controlled_value("polarity", _normalize_scalar_value(experiment.get("polarity")) or "unknown"),
        "ionisation_source": _normalize_scalar_value(experiment.get("ionisation_source")),
        "analyzer": _normalize_scalar_value(experiment.get("analyzer")),
        "instrument": _normalize_scalar_value(experiment.get("instrument")),
        "spectrum_mode": _normalize_scalar_value(experiment.get("spectrum_mode")),
        "storage_mode": _normalize_scalar_value(experiment.get("storage_mode")),
        "acquisition_mode": _normalize_scalar_value(experiment.get("acquisition_mode")),
        "pixel_size_x_um": _normalize_numeric_value(acquisition.get("pixel_size_x_um")),
        "pixel_size_y_um": _normalize_numeric_value(acquisition.get("pixel_size_y_um")),
        "sample_stabilization": _normalize_scalar_value(preparation.get("sample_stabilization")),
        "solvent": _normalize_scalar_value(preparation.get("solvent")),
        "prep_method": _normalize_scalar_value(preparation.get("prep_method")),
        "submitter_notes": _normalize_scalar_value(notes.get("submitter_notes")) if notes else None,
    }

    _apply_biological_metadata_defaults(canonical)
    missing_fields = [field_name for field_name in REQUIRED_BIOLOGICAL_FIELDS if not _lookup_required_value(canonical, field_name)]
    if missing_fields:
        missing_text = ", ".join(missing_fields)
        raise BiologicalMetadataValidationError(
            f"Missing required biological metadata fields: {missing_text}."
        )

    return canonical


def build_biological_file_binding(
    payload: Mapping[str, Any],
    *,
    imzml_path: str | None = None,
) -> BiologicalFileBinding | None:
    """Validate or infer the file binding between metadata and the target imzML."""

    binding_payload = payload.get("binding")
    declared_file_name: str | None = None
    if isinstance(binding_payload, Mapping):
        declared_file_name = _normalize_scalar_value(
            binding_payload.get("imzml_file_name") or binding_payload.get("imzml_basename")
        )

    if declared_file_name is None:
        declared_file_name = _normalize_scalar_value(
            payload.get("imzml_file_name") or payload.get("imzml_basename")
        )

    inferred_file_name = Path(imzml_path).name if imzml_path else None
    if declared_file_name and inferred_file_name and declared_file_name != inferred_file_name:
        raise BiologicalMetadataValidationError(
            "Biological metadata binding does not match the target imzML file name: "
            f"{declared_file_name} != {inferred_file_name}."
        )

    bound_file_name = declared_file_name or inferred_file_name
    if not bound_file_name:
        return None

    return BiologicalFileBinding(imzml_file_name=bound_file_name)


def build_biological_metadata_template() -> dict[str, Any]:
    """Return a JSON-serializable metadata template for upload clients."""

    return {
        "schema_version": BIOLOGICAL_METADATA_SCHEMA_VERSION,
        "sample": {
            "sample_id": "",
            "organism": "",
            "organism_part": "",
            "tissue_region": "unknown",
            "sample_type": "tissue",
            "condition": "unknown",
            "group_label": "unknown",
        },
        "experiment": {
            "polarity": "unknown",
            "ionisation_source": "unknown",
            "analyzer": "unknown",
            "instrument": "",
            "spectrum_mode": "unknown",
            "storage_mode": "unknown",
            "acquisition_mode": "MSI",
        },
        "acquisition": {
            "pixel_size_x_um": None,
            "pixel_size_y_um": None,
        },
        "preparation": {
            "sample_stabilization": "unknown",
            "solvent": "unknown",
            "prep_method": "unknown",
        },
        "notes": {
            "submitter_notes": "",
        },
        "binding": {
            "imzml_file_name": "",
        },
    }


def _find_first_present_value(payload: Mapping[str, Any], aliases: tuple[str, ...]) -> Any:
    for alias in aliases:
        if alias in payload:
            return payload[alias]
    return None


def _normalize_scalar_value(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _normalize_controlled_value(field_name: str, value: str) -> str:
    vocabulary = CONTROLLED_VOCABULARY.get(field_name)
    if not vocabulary:
        return value
    return vocabulary.get(value.casefold(), value)


def _extract_additional_fields(payload: Mapping[str, Any]) -> dict[str, Any]:
    additional_fields: dict[str, Any] = {}
    for section_name, known_fields in KNOWN_SECTION_FIELDS.items():
        section_value = payload.get(section_name)
        if isinstance(section_value, Mapping):
            extra_items = {
                str(key): value for key, value in section_value.items() if str(key) not in known_fields
            }
            if extra_items:
                additional_fields[section_name] = extra_items

    for key, value in payload.items():
        if str(key) in KNOWN_SECTION_FIELDS or str(key) == "schema_version":
            continue
        additional_fields[str(key)] = value
    return additional_fields


def _normalize_numeric_value(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise BiologicalMetadataValidationError(
            f"Biological metadata numeric field could not be parsed: {value!r}."
        ) from exc


def _apply_biological_metadata_defaults(canonical: dict[str, Any]) -> None:
    for field_name, default_value in BIOLOGICAL_METADATA_DEFAULTS.items():
        canonical.setdefault(field_name, default_value)


def _require_mapping_section(payload: Mapping[str, Any], section_name: str) -> Mapping[str, Any]:
    section_value = payload.get(section_name)
    if not isinstance(section_value, Mapping):
        raise BiologicalMetadataValidationError(
            f"Biological metadata section '{section_name}' must be a JSON object."
        )
    return section_value


def _optional_mapping_section(payload: Mapping[str, Any], section_name: str) -> Mapping[str, Any] | None:
    section_value = payload.get(section_name)
    if section_value is None:
        return None
    if not isinstance(section_value, Mapping):
        raise BiologicalMetadataValidationError(
            f"Biological metadata section '{section_name}' must be a JSON object when provided."
        )
    return section_value


def _lookup_required_value(canonical: Mapping[str, Any], dotted_name: str) -> Any:
    _, field_name = dotted_name.split(".", 1)
    return canonical.get(field_name)
