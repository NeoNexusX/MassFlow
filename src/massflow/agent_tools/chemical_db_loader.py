from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable
import csv


DEFAULT_CHEMICAL_DB_FIELD_ALIASES: dict[str, tuple[str, ...]] = {
    "source_record_id": (
        "source_record_id",
        "lm_id",
        "LM_ID",
        "accession",
        "hmdb_id",
        "HMDB_ID",
        "chebi_id",
        "CHEBI_ID",
        "id",
        "ID",
    ),
    "compound_name": (
        "compound_name",
        "name",
        "NAME",
        "Common Name",
        "common_name",
        "SYSTEMATIC_NAME",
    ),
    "formula": (
        "formula",
        "FORMULA",
        "chemical_formula",
        "Chemical Formula",
    ),
    "exact_mass": (
        "exact_mass",
        "exactmass",
        "EXACT_MASS",
        "monoisotopic_mass",
        "Monoisotopic Mass",
        "monoisotopic_molecular_weight",
    ),
    "smiles": (
        "smiles",
        "SMILES",
        "canonical_smiles",
        "Canonical SMILES",
        "SMILES (pH7)",
    ),
    "inchi_key": (
        "inchi_key",
        "INCHI_KEY",
        "inchikey",
        "InChIKey",
    ),
    "class_label": (
        "class_label",
        "CATEGORY",
        "MAIN_CLASS",
        "SUB_CLASS",
        "Class",
        "class",
    ),
}


class ChemicalLibraryValidationError(ValueError):
    """Raised when local chemical library records are malformed."""


@dataclass(slots=True)
class ChemicalLibraryEntry:
    """Normalized entry loaded from a local chemical annotation database."""

    source_name: str
    source_record_id: str | None
    compound_name: str | None
    formula: str
    exact_mass: float
    smiles: str | None = None
    inchi_key: str | None = None
    class_label: str | None = None
    additional_fields: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_chemical_library_table(
    file_path: str,
    *,
    source_name: str,
    delimiter: str | None = None,
    field_aliases: dict[str, tuple[str, ...]] | None = None,
) -> list[ChemicalLibraryEntry]:
    """Load a CSV/TSV chemical library export into normalized entries."""

    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Chemical library file not found: {path}")

    resolved_delimiter = delimiter or _infer_delimiter(path)
    aliases = field_aliases or DEFAULT_CHEMICAL_DB_FIELD_ALIASES
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=resolved_delimiter)
        rows = list(reader)
    if not rows:
        raise ChemicalLibraryValidationError("Chemical library table is empty.")

    header_index = _detect_header_row(rows, field_aliases=aliases)
    header = rows[header_index]
    data_rows = rows[header_index + 1 :]
    dict_rows = [dict(zip(header, row, strict=False)) for row in data_rows]

    entries: list[ChemicalLibraryEntry] = []
    for row in dict_rows:
        if not any(_normalize_text(value) is not None for value in row.values()):
            continue
        try:
            entries.append(
                build_chemical_library_entry(row, source_name=source_name, field_aliases=aliases)
            )
        except ChemicalLibraryValidationError:
            continue

    if not entries:
        raise ChemicalLibraryValidationError(
            "No valid chemical library entries were found after parsing the table."
        )
    return entries


def build_chemical_library_entry(
    row: dict[str, Any],
    *,
    source_name: str,
    field_aliases: dict[str, tuple[str, ...]] | None = None,
) -> ChemicalLibraryEntry:
    """Normalize one row from a downloaded chemical library table."""

    aliases = field_aliases or DEFAULT_CHEMICAL_DB_FIELD_ALIASES
    formula = _extract_first_value(row, aliases["formula"])
    exact_mass_raw = _extract_first_value(row, aliases["exact_mass"])
    if formula is None:
        raise ChemicalLibraryValidationError("Chemical library entry requires a formula column.")
    if exact_mass_raw is None:
        raise ChemicalLibraryValidationError("Chemical library entry requires an exact mass column.")

    try:
        exact_mass = float(exact_mass_raw)
    except (TypeError, ValueError) as exc:
        raise ChemicalLibraryValidationError(
            f"Chemical library entry exact mass must be numeric, got {exact_mass_raw!r}."
        ) from exc

    return ChemicalLibraryEntry(
        source_name=source_name,
        source_record_id=_extract_first_value(row, aliases["source_record_id"]),
        compound_name=_extract_first_value(row, aliases["compound_name"]),
        formula=formula,
        exact_mass=exact_mass,
        smiles=_extract_first_value(row, aliases["smiles"]),
        inchi_key=_extract_first_value(row, aliases["inchi_key"]),
        class_label=_extract_first_value(row, aliases["class_label"]),
        additional_fields=_extract_additional_fields(row, aliases),
    )


def summarize_chemical_library(entries: Iterable[ChemicalLibraryEntry]) -> dict[str, Any]:
    """Return a lightweight summary for a loaded local chemical library."""

    entries_list = list(entries)
    source_names = sorted({entry.source_name for entry in entries_list})
    formulas = sorted({entry.formula for entry in entries_list if entry.formula})
    return {
        "entry_count": len(entries_list),
        "source_names": source_names,
        "unique_formula_count": len(formulas),
        "entries_with_smiles": sum(1 for entry in entries_list if entry.smiles),
    }


def load_chemical_library_directory(
    directory_path: str,
    *,
    source_name: str = "auto",
    patterns: tuple[str, ...] = ("*.tsv", "*.csv"),
    recursive: bool = True,
) -> list[ChemicalLibraryEntry]:
    """Load and combine all supported chemical library tables under one directory."""

    root = Path(directory_path)
    if not root.exists():
        raise FileNotFoundError(f"Chemical library directory not found: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Chemical library path is not a directory: {root}")

    files: list[Path] = []
    for pattern in patterns:
        if recursive:
            files.extend(sorted(root.rglob(pattern)))
        else:
            files.extend(sorted(root.glob(pattern)))
    if not files:
        raise ChemicalLibraryValidationError(
            f"No supported chemical library files were found under {root}."
        )

    entries: list[ChemicalLibraryEntry] = []
    for file_path in files:
        entries.extend(
            load_chemical_library_table(
                str(file_path),
                source_name=_resolve_source_name(source_name, file_path),
            )
        )
    if not entries:
        raise ChemicalLibraryValidationError(
            f"No valid chemical library entries were loaded from directory {root}."
        )
    return entries


def _infer_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".tsv":
        return "\t"
    return ","


def _detect_header_row(
    rows: list[list[str]],
    *,
    field_aliases: dict[str, tuple[str, ...]],
) -> int:
    required_aliases = set(field_aliases["formula"]) | set(field_aliases["exact_mass"])
    for index, row in enumerate(rows[:10]):
        if required_aliases.intersection(set(row)):
            return index
    raise ChemicalLibraryValidationError(
        "Failed to detect a supported header row in the first 10 lines of the chemical library table."
    )


def _extract_first_value(row: dict[str, Any], aliases: tuple[str, ...]) -> str | None:
    for alias in aliases:
        if alias in row:
            normalized = _normalize_text(row[alias])
            if normalized is not None:
                return normalized
    return None


def _normalize_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _extract_additional_fields(
    row: dict[str, Any],
    aliases: dict[str, tuple[str, ...]],
) -> dict[str, Any]:
    known_fields = {alias for options in aliases.values() for alias in options}
    return {key: value for key, value in row.items() if key not in known_fields}


def _resolve_source_name(source_name: str, file_path: Path) -> str:
    if source_name != "auto":
        return source_name

    normalized = str(file_path).lower()
    if "lipidmaps" in normalized or file_path.stem.upper().startswith("LM"):
        return "LIPID MAPS"
    if "chebi" in normalized:
        return "ChEBI"
    if "hmdb" in normalized:
        return "HMDB"
    return file_path.parent.name or file_path.stem
