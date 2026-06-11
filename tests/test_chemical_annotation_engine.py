from pathlib import Path

from massflow.agent_tools.chemical_annotation_engine import annotate_mz_values, annotate_single_mz
from massflow.agent_tools.chemical_db_loader import (
    ChemicalLibraryValidationError,
    load_chemical_library_table,
    summarize_chemical_library,
)


def test_load_chemical_library_table_from_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "lipidmaps_sample.csv"
    csv_path.write_text(
        "LM_ID,NAME,FORMULA,EXACT_MASS,SMILES,INCHI_KEY,CATEGORY\n"
        "LMFA0001,Demo Lipid,C42H82NO8P,759.5777,CCCC,DEMO-INCHI,Lipids\n",
        encoding="utf-8",
    )

    entries = load_chemical_library_table(str(csv_path), source_name="LIPID MAPS")
    summary = summarize_chemical_library(entries)

    assert len(entries) == 1
    assert entries[0].source_record_id == "LMFA0001"
    assert entries[0].compound_name == "Demo Lipid"
    assert summary["entries_with_smiles"] == 1


def test_load_chemical_library_table_requires_exact_mass(tmp_path: Path) -> None:
    csv_path = tmp_path / "invalid_library.csv"
    csv_path.write_text(
        "LM_ID,NAME,FORMULA\n"
        "LMFA0001,Demo Lipid,C42H82NO8P\n",
        encoding="utf-8",
    )

    try:
        load_chemical_library_table(str(csv_path), source_name="LIPID MAPS")
    except ChemicalLibraryValidationError as exc:
        assert "exact mass" in str(exc)
    else:
        raise AssertionError("Expected ChemicalLibraryValidationError for missing exact mass.")


def test_annotate_single_mz_matches_local_library() -> None:
    entries = load_chemical_library_table_from_rows(
        [
            {
                "LM_ID": "LMFA0001",
                "NAME": "Demo Lipid",
                "FORMULA": "C42H82NO8P",
                "EXACT_MASS": "759.5777",
                "SMILES": "CCCC",
            }
        ]
    )

    record = annotate_single_mz(
        760.584976466812,
        library_entries=entries,
        ion_mode="positive",
        ppm_tolerance=5.0,
    )

    assert record.candidates[0].formula == "C42H82NO8P"
    assert record.candidates[0].adduct == "[M+H]+"
    assert record.candidates[0].source_name == "LIPID MAPS"


def test_annotate_mz_values_returns_unknown_when_no_match() -> None:
    entries = load_chemical_library_table_from_rows(
        [
            {
                "LM_ID": "LMFA0001",
                "NAME": "Demo Lipid",
                "FORMULA": "C42H82NO8P",
                "EXACT_MASS": "759.5777",
            }
        ]
    )

    annotations = annotate_mz_values(
        [100.0],
        library_entries=entries,
        ion_mode="positive",
        ppm_tolerance=5.0,
    )

    assert len(annotations) == 1
    assert annotations[0].candidates[0].formula == "unannotated"
    assert annotations[0].candidates[0].confidence_level == "unknown"


def load_chemical_library_table_from_rows(rows: list[dict[str, str]]):
    import csv
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".csv", encoding="utf-8", newline="", delete=False) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
        temp_path = handle.name

    try:
        return load_chemical_library_table(temp_path, source_name="LIPID MAPS")
    finally:
        Path(temp_path).unlink(missing_ok=True)
