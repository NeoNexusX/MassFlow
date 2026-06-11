from massflow.agent_tools.biological_schema import (
    BiologicalMetadata,
    BiologicalMetadataValidationError,
    build_biological_metadata_template,
)


def test_biological_metadata_accepts_grouped_schema_and_binding() -> None:
    metadata = BiologicalMetadata.from_mapping(
        {
            "sample": {
                "sample_id": "Sample-001",
                "organism": "Mus musculus",
                "organism_part": "brain",
                "tissue_region": "cortex",
                "sample_type": "fresh frozen",
                "condition": "control",
                "group_label": "group-a",
            },
            "experiment": {
                "polarity": "POS",
                "ionisation_source": "MALDI",
                "analyzer": "TOF",
            },
            "preparation": {
                "prepMethod": "MALDI matrix spray",
                "sample_stabilization": "fresh_frozen",
            },
            "acquisition": {
                "pixel_size_x_um": 50,
                "pixel_size_y_um": "75",
            },
            "binding": {"imzml_file_name": "mouse_brain.imzML"},
            "operator": "alice",
        },
        imzml_path="C:/data/mouse_brain.imzML",
    )

    assert metadata.sample_id == "Sample-001"
    assert metadata.sample_type == "fresh_frozen"
    assert metadata.ion_mode == "positive"
    assert metadata.organism == "Mus musculus"
    assert metadata.organism_part == "brain"
    assert metadata.polarity == "positive"
    assert metadata.ionisation_source == "MALDI"
    assert metadata.analyzer == "TOF"
    assert metadata.pixel_size_x_um == 50.0
    assert metadata.pixel_size_y_um == 75.0
    assert metadata.binding is not None
    assert metadata.binding.imzml_file_name == "mouse_brain.imzML"
    assert metadata.additional_fields == {"operator": "alice"}


def test_biological_metadata_rejects_missing_required_fields() -> None:
    try:
        BiologicalMetadata.from_mapping(
            {
                "sample": {"sample_id": "Sample-001"},
                "experiment": {},
                "acquisition": {},
                "preparation": {},
                "binding": {"imzml_file_name": "mouse_brain.imzML"},
            }
        )
    except BiologicalMetadataValidationError as exc:
        assert "sample.organism" in str(exc)
        assert "sample.organism_part" in str(exc)
    else:
        raise AssertionError("Expected BiologicalMetadataValidationError for missing fields.")


def test_biological_metadata_rejects_old_flat_schema() -> None:
    try:
        BiologicalMetadata.from_mapping(
            {
                "sample_id": "Sample-001",
                "species": "Mus musculus",
                "organ": "brain",
            }
        )
    except BiologicalMetadataValidationError as exc:
        assert "grouped JSON schema" in str(exc)
    else:
        raise AssertionError("Expected BiologicalMetadataValidationError for old flat schema.")


def test_biological_metadata_rejects_binding_mismatch() -> None:
    try:
        BiologicalMetadata.from_mapping(
            {
                "sample": {
                    "sample_id": "Sample-001",
                    "organism": "Homo sapiens",
                    "organism_part": "liver",
                },
                "experiment": {
                    "polarity": "negative",
                },
                "preparation": {
                    "prep_method": "cryosection",
                },
                "binding": {"imzml_file_name": "other.imzML"},
            },
            imzml_path="C:/data/real.imzML",
        )
    except BiologicalMetadataValidationError as exc:
        assert "does not match" in str(exc)
    else:
        raise AssertionError("Expected BiologicalMetadataValidationError for binding mismatch.")


def test_biological_metadata_template_contains_required_keys() -> None:
    template = build_biological_metadata_template()

    assert template["schema_version"]
    assert template["sample"]["sample_id"] == ""
    assert template["experiment"]["polarity"] == "unknown"
    assert template["acquisition"]["pixel_size_x_um"] is None
    assert template["binding"]["imzml_file_name"] == ""
