from massflow.agent_tools.chemical_schema import (
    ChemicalAnnotationRecord,
    ChemicalAnnotationValidationError,
    build_chemical_annotation_context,
    build_chemical_annotation_template,
)


def test_chemical_annotation_record_builds_from_mapping() -> None:
    record = ChemicalAnnotationRecord.from_mapping(
        {
            "mz": 760.585,
            "intensity": 12345.0,
            "rank_strategy": "ppm_error_then_score",
            "ppm_tolerance": 5.0,
            "ion_mode": "positive",
            "query_label": "top_peak_1",
            "candidates": [
                {
                    "formula": "C42H82NO8P",
                    "smiles": "CCCCCCCCCCCCCCCC(=O)O",
                    "compound_name": "PC(34:1)",
                    "adduct": "[M+H]+",
                    "charge": 1,
                    "ppm_error": 1.2,
                    "score": 0.97,
                    "confidence_level": "high",
                    "source_type": "public_database",
                    "source_name": "LIPID MAPS",
                    "source_record_id": "LMGP01010001",
                }
            ],
        }
    )

    assert record.mz == 760.585
    assert record.candidates[0].formula == "C42H82NO8P"
    assert record.candidates[0].source_type == "public_database"


def test_chemical_annotation_context_summarizes_top_hits() -> None:
    record = ChemicalAnnotationRecord.from_mapping(
        {
            "mz": 184.0733,
            "candidates": [
                {
                    "formula": "C5H15NO4P",
                    "compound_name": "Phosphocholine",
                    "confidence_level": "medium",
                    "source_type": "rule_based",
                    "source_name": "diagnostic_ions",
                }
            ],
        }
    )

    context = build_chemical_annotation_context([record])
    assert context["provided"] is True
    assert context["annotation_count"] == 1
    assert context["top_hits"][0]["formula"] == "C5H15NO4P"


def test_chemical_annotation_record_rejects_invalid_candidate() -> None:
    try:
        ChemicalAnnotationRecord.from_mapping(
            {
                "mz": 500.0,
                "candidates": [
                    {
                        "compound_name": "Missing formula",
                    }
                ],
            }
        )
    except ChemicalAnnotationValidationError as exc:
        assert "formula" in str(exc)
    else:
        raise AssertionError("Expected ChemicalAnnotationValidationError for missing formula.")


def test_chemical_annotation_template_contains_database_source() -> None:
    template = build_chemical_annotation_template()
    assert template["schema_version"]
    assert template["candidates"][0]["source_type"] == "public_database"
