from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import gzip


def build_chebi_annotation_table(
    chebi_directory: str,
    *,
    output_file_path: str | None = None,
    min_stars: int = 2,
) -> str:
    """Build one annotation table from ChEBI flat files for local m/z matching."""

    root = Path(chebi_directory)
    compounds_path = root / "compounds.tsv.gz"
    chemical_data_path = root / "chemical_data.tsv.gz"
    structures_path = root / "structures.tsv.gz"
    for path in (compounds_path, chemical_data_path, structures_path):
        if not path.exists():
            raise FileNotFoundError(f"Required ChEBI file not found: {path}")

    output_path = Path(output_file_path) if output_file_path else root / "chebi_annotation_table.tsv"
    compounds = _load_chebi_compounds(compounds_path, min_stars=min_stars)
    chemical_data = _load_chebi_chemical_data(chemical_data_path)
    structures = _load_chebi_structures(structures_path)

    fieldnames = [
        "source_record_id",
        "compound_name",
        "formula",
        "exact_mass",
        "smiles",
        "inchi_key",
        "class_label",
        "stars",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for compound_id, compound in compounds.items():
            chemical = chemical_data.get(compound_id)
            structure = structures.get(compound_id)
            if chemical is None:
                continue
            exact_mass = chemical.get("monoisotopic_mass") or chemical.get("mass")
            formula = chemical.get("formula")
            if not formula or not exact_mass:
                continue
            writer.writerow(
                {
                    "source_record_id": compound["chebi_accession"],
                    "compound_name": compound["name"],
                    "formula": formula,
                    "exact_mass": exact_mass,
                    "smiles": None if structure is None else structure.get("smiles"),
                    "inchi_key": None if structure is None else structure.get("standard_inchi_key"),
                    "class_label": "ChEBI",
                    "stars": compound["stars"],
                }
            )
    return str(output_path)


def _load_chebi_compounds(file_path: Path, *, min_stars: int) -> dict[str, dict[str, Any]]:
    compounds: dict[str, dict[str, Any]] = {}
    with gzip.open(file_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                stars = int(row.get("stars") or 0)
            except ValueError:
                continue
            if stars < min_stars:
                continue
            compound_id = str(row["id"])
            compounds[compound_id] = {
                "name": row.get("name") or row.get("ascii_name"),
                "chebi_accession": row.get("chebi_accession"),
                "stars": stars,
            }
    return compounds


def _load_chebi_chemical_data(file_path: Path) -> dict[str, dict[str, Any]]:
    chemical_data: dict[str, dict[str, Any]] = {}
    with gzip.open(file_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            compound_id = str(row["compound_id"])
            if compound_id not in chemical_data:
                chemical_data[compound_id] = row
    return chemical_data


def _load_chebi_structures(file_path: Path) -> dict[str, dict[str, Any]]:
    structures: dict[str, dict[str, Any]] = {}
    with gzip.open(file_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            compound_id = str(row["compound_id"])
            if compound_id in structures:
                continue
            structures[compound_id] = row
    return structures
