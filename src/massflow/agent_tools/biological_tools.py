from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping
import json

from .biological_schema import BiologicalMetadata


def resolve_biological_metadata(
    *,
    biological_metadata: Mapping[str, Any] | BiologicalMetadata | None = None,
    biological_metadata_json: str | None = None,
    biological_metadata_path: str | None = None,
    imzml_path: str | None = None,
) -> BiologicalMetadata | None:
    """Resolve biological metadata from object, JSON string, or JSON file path."""

    provided_count = sum(
        value is not None and value != ""
        for value in (biological_metadata, biological_metadata_json, biological_metadata_path)
    )
    if provided_count > 1:
        raise ValueError(
            "Provide biological metadata through exactly one source: "
            "biological_metadata, biological_metadata_json, or biological_metadata_path."
        )

    if isinstance(biological_metadata, BiologicalMetadata):
        return biological_metadata
    if biological_metadata is not None:
        return BiologicalMetadata.from_mapping(biological_metadata, imzml_path=imzml_path)
    if biological_metadata_json:
        return BiologicalMetadata.from_mapping(
            _parse_biological_metadata_json(biological_metadata_json, source_name="biological_metadata_json"),
            imzml_path=imzml_path,
        )
    if biological_metadata_path:
        return BiologicalMetadata.from_mapping(
            _parse_biological_metadata_file(biological_metadata_path),
            imzml_path=imzml_path,
        )
    return None


def build_biological_context(
    metadata: BiologicalMetadata | None,
    *,
    source: str | None = None,
) -> dict[str, Any]:
    """Build an agent-friendly context object from canonical biological metadata."""

    if metadata is None:
        return {
            "provided": False,
            "source": source,
            "schema_version": None,
            "sample": None,
            "experiment": None,
            "binding": None,
            "additional_fields": {},
        }

    return {
        "provided": True,
        "source": source,
        "schema_version": metadata.schema_version,
        "sample": {
            "sample_id": metadata.sample_id,
            "species": metadata.species,
            "organism": metadata.organism,
            "organ": metadata.organ,
            "organism_part": metadata.organism_part,
            "tissue_region": metadata.tissue_region,
            "sample_type": metadata.sample_type,
            "condition": metadata.condition,
            "group_label": metadata.group_label,
        },
        "experiment": {
            "ion_mode": metadata.ion_mode,
            "polarity": metadata.polarity,
            "prep_method": metadata.prep_method,
            "ionisation_source": metadata.ionisation_source,
            "analyzer": metadata.analyzer,
            "instrument": metadata.instrument,
            "acquisition_mode": metadata.acquisition_mode,
            "spectrum_mode": metadata.spectrum_mode,
            "storage_mode": metadata.storage_mode,
            "submitter_notes": metadata.submitter_notes,
        },
        "acquisition": {
            "pixel_size_x_um": metadata.pixel_size_x_um,
            "pixel_size_y_um": metadata.pixel_size_y_um,
        },
        "preparation": {
            "sample_stabilization": metadata.sample_stabilization,
            "solvent": metadata.solvent,
            "prep_method": metadata.prep_method,
        },
        "binding": None if metadata.binding is None else asdict(metadata.binding),
        "additional_fields": metadata.additional_fields,
    }


def _parse_biological_metadata_json(raw_json: str, *, source_name: str) -> Mapping[str, Any]:
    try:
        payload = json.loads(raw_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{source_name} must be valid JSON: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"{source_name} must decode to a JSON object.")
    return payload


def _parse_biological_metadata_file(file_path: str) -> Mapping[str, Any]:
    path = Path(file_path)
    payload = _parse_biological_metadata_json(
        path.read_text(encoding="utf-8-sig"),
        source_name=f"biological metadata file {path}",
    )
    return payload
