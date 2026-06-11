from __future__ import annotations

import json
from typing import Any

from massflow.agent_tools import (
    analyze_imzml_end_to_end,
    apply_spatial_png_mask,
    build_spatial_preview,
    dataset_qc_summary,
    load_imzml,
    recommend_preprocessing_plan,
    run_preprocessing,
    validate_integrated_conclusion,
)
from massflow.agent_tools._utils import make_json_safe


def register_massflow_tools(server: Any) -> Any:
    """Register the stable MassFlow MCP tools on a FastMCP server instance."""

    @server.tool(
        name="load_imzml",
        description=(
            "Load imzML header and dataset index only. "
            "Use this first to inspect dataset size, coordinates, occupancy, "
            "and optional structured biological metadata."
        ),
    )
    def load_imzml_tool(
        imzml_path: str,
        output_dir: str | None = None,
        biological_metadata_json: str | None = None,
        biological_metadata_path: str | None = None,
        max_threads: int = 8,
    ) -> dict[str, Any]:
        result = load_imzml(
            imzml_path,
            output_dir=output_dir,
            biological_metadata_json=biological_metadata_json,
            biological_metadata_path=biological_metadata_path,
            max_threads=max_threads,
        )
        return make_json_safe(result.to_dict())

    @server.tool(
        name="dataset_qc_summary",
        description=(
            "Sample spectra from one imzML dataset and return lightweight QC metrics, "
            "including TIC, peak counts, SNR summaries, QC flags, and optional "
            "biological context."
        ),
    )
    def dataset_qc_summary_tool(
        imzml_path: str,
        output_dir: str | None = None,
        biological_metadata_json: str | None = None,
        biological_metadata_path: str | None = None,
        sample_size: int = 16,
        noise_method: str = "sd",
        max_threads: int = 8,
    ) -> dict[str, Any]:
        result = dataset_qc_summary(
            imzml_path,
            output_dir=output_dir,
            biological_metadata_json=biological_metadata_json,
            biological_metadata_path=biological_metadata_path,
            sample_size=sample_size,
            noise_method=noise_method,
            max_threads=max_threads,
        )
        return make_json_safe(result.to_dict())

    @server.tool(
        name="recommend_preprocessing_plan",
        description=(
            "Inspect one imzML dataset, run lightweight QC, and recommend which preprocessing "
            "operations should be applied before integrated conclusion generation."
        ),
    )
    def recommend_preprocessing_plan_tool(
        imzml_path: str,
        output_dir: str | None = None,
        biological_metadata_json: str | None = None,
        biological_metadata_path: str | None = None,
        sample_size: int = 16,
        noise_method: str = "sd",
        max_threads: int = 8,
    ) -> dict[str, Any]:
        result = recommend_preprocessing_plan(
            imzml_path,
            output_dir=output_dir,
            biological_metadata_json=biological_metadata_json,
            biological_metadata_path=biological_metadata_path,
            sample_size=sample_size,
            noise_method=noise_method,
            max_threads=max_threads,
        )
        return make_json_safe(result.to_dict())

    @server.tool(
        name="run_preprocessing",
        description=(
            "Run the stable preprocessing pipeline for one imzML dataset. "
            "Returns before/after QC comparison and persisted processed imzML paths."
        ),
    )
    def run_preprocessing_tool(
        imzml_path: str,
        output_dir: str | None = None,
        operations_json: str | None = None,
        batch_size: int = 256,
        max_threads: int = 8,
        qc_sample_size: int = 12,
    ) -> dict[str, Any]:
        operations = _parse_json_argument(
            operations_json,
            expected_type=list,
            argument_name="operations_json",
        )
        result = run_preprocessing(
            imzml_path,
            output_dir=output_dir,
            operations=operations,
            batch_size=batch_size,
            max_threads=max_threads,
            qc_sample_size=qc_sample_size,
        )
        return make_json_safe(result.to_dict())

    @server.tool(
        name="build_spatial_preview",
        description=(
            "Generate basic spatial preview images for one imzML dataset, including TIC and "
            "base-peak PNG files. Use these images to manually create a PNG mask for later "
            "matrix/background removal."
        ),
    )
    def build_spatial_preview_tool(
        imzml_path: str,
        output_dir: str | None = None,
        max_threads: int = 8,
        colormap: str = "magma",
    ) -> dict[str, Any]:
        result = build_spatial_preview(
            imzml_path,
            output_dir=output_dir,
            max_threads=max_threads,
            colormap=colormap,
        )
        return make_json_safe(result.to_dict())

    @server.tool(
        name="apply_spatial_png_mask",
        description=(
            "Apply a manually edited PNG mask to one imzML dataset and return the retained "
            "pixel coordinates. White pixels are treated as tissue by default."
        ),
    )
    def apply_spatial_png_mask_tool(
        imzml_path: str,
        mask_png_path: str,
        output_dir: str | None = None,
        max_threads: int = 8,
        threshold: float = 0.5,
        keep_mode: str = "white_foreground",
    ) -> dict[str, Any]:
        result = apply_spatial_png_mask(
            imzml_path,
            mask_png_path,
            output_dir=output_dir,
            max_threads=max_threads,
            threshold=threshold,
            keep_mode=keep_mode,
        )
        return make_json_safe(result.to_dict())

    @server.tool(
        name="analyze_imzml_end_to_end",
        description=(
            "Run the current stable workflow: load, QC, optional preprocessing, "
            "chemical annotation, and integrated conclusion generation. "
            "Supports optional structured biological metadata input and local chemical annotation."
        ),
    )
    def analyze_imzml_end_to_end_tool(
        imzml_path: str,
        output_dir: str | None = None,
        biological_metadata_json: str | None = None,
        biological_metadata_path: str | None = None,
        chemical_library_dir: str | None = None,
        chemical_source_name: str = "auto",
        chemical_ion_mode: str | None = None,
        chemical_ppm_tolerance: float = 5.0,
        chemical_max_mz_values: int = 12,
        chemical_max_candidates_per_mz: int = 5,
        run_preprocessing_first: bool = True,
        preprocessing_operations_json: str | None = None,
    ) -> dict[str, Any]:
        preprocessing_operations = _parse_json_argument(
            preprocessing_operations_json,
            expected_type=list,
            argument_name="preprocessing_operations_json",
        )
        result = analyze_imzml_end_to_end(
            imzml_path,
            output_dir=output_dir,
            biological_metadata_json=biological_metadata_json,
            biological_metadata_path=biological_metadata_path,
            chemical_library_dir=chemical_library_dir,
            chemical_source_name=chemical_source_name,
            chemical_ion_mode=chemical_ion_mode,
            chemical_ppm_tolerance=chemical_ppm_tolerance,
            chemical_max_mz_values=chemical_max_mz_values,
            chemical_max_candidates_per_mz=chemical_max_candidates_per_mz,
            run_preprocessing_first=run_preprocessing_first,
            preprocessing_operations=preprocessing_operations,
        )
        return make_json_safe(result.to_dict())

    @server.tool(
        name="validate_integrated_conclusion",
        description=(
            "Run the integrated workflow and compare its conclusion with a legacy statistical-only baseline. "
            "Useful for regression testing and validating the evidence shift from statistics-first to integrated reasoning."
        ),
    )
    def validate_integrated_conclusion_tool(
        imzml_path: str,
        output_dir: str | None = None,
        biological_metadata_json: str | None = None,
        biological_metadata_path: str | None = None,
        chemical_library_dir: str | None = None,
        chemical_source_name: str = "auto",
        chemical_ion_mode: str | None = None,
        chemical_ppm_tolerance: float = 5.0,
        chemical_max_mz_values: int = 12,
        chemical_max_candidates_per_mz: int = 5,
        run_preprocessing_first: bool = False,
        preprocessing_operations_json: str | None = None,
    ) -> dict[str, Any]:
        preprocessing_operations = _parse_json_argument(
            preprocessing_operations_json,
            expected_type=list,
            argument_name="preprocessing_operations_json",
        )
        result = validate_integrated_conclusion(
            imzml_path,
            output_dir=output_dir,
            biological_metadata_json=biological_metadata_json,
            biological_metadata_path=biological_metadata_path,
            chemical_library_dir=chemical_library_dir,
            chemical_source_name=chemical_source_name,
            chemical_ion_mode=chemical_ion_mode,
            chemical_ppm_tolerance=chemical_ppm_tolerance,
            chemical_max_mz_values=chemical_max_mz_values,
            chemical_max_candidates_per_mz=chemical_max_candidates_per_mz,
            run_preprocessing_first=run_preprocessing_first,
            preprocessing_operations=preprocessing_operations,
        )
        return make_json_safe(result.to_dict())

    return server


def _parse_json_argument(
    raw_value: str | None,
    *,
    expected_type: type,
    argument_name: str,
) -> Any:
    if raw_value is None or raw_value == "":
        return None

    try:
        parsed = json.loads(raw_value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{argument_name} must be valid JSON: {exc}") from exc

    if not isinstance(parsed, expected_type):
        raise ValueError(
            f"{argument_name} must decode to {expected_type.__name__}, "
            f"got {type(parsed).__name__}."
        )
    return parsed
