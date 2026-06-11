from pathlib import Path

from massflow.agent_tools.chemical_tools import annotate_mz_with_local_library


def test_annotate_mz_with_local_library_directory(tmp_path: Path) -> None:
    library_dir = tmp_path / "lipidmaps"
    library_dir.mkdir()
    (library_dir / "LMGP.tsv").write_text(
        "06/02/26\n"
        "regno\tlm_id\tname\texactmass\tformula\tsmiles\n"
        "1\tLMGP04010070\tDemo Lipid\t759.5777\tC42H82NO8P\tCCCC\n",
        encoding="utf-8",
    )

    result = annotate_mz_with_local_library(
        [760.584976466812],
        library_dir=str(library_dir),
        source_name="LIPID MAPS",
        ion_mode="positive",
        ppm_tolerance=5.0,
    )

    assert result.ok is True
    assert result.details["chemical_context"]["provided"] is True
    assert result.details["annotations"][0]["candidates"][0]["formula"] == "C42H82NO8P"
