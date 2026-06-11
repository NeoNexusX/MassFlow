"""Agent-friendly tool layer for MassFlow."""

from .biological_schema import (
    BIOLOGICAL_METADATA_SCHEMA_VERSION,
    REQUIRED_BIOLOGICAL_FIELDS,
    BiologicalMetadata,
    BiologicalMetadataValidationError,
    build_biological_metadata_template,
)
from .biological_tools import build_biological_context, resolve_biological_metadata
from .chemical_annotation_engine import (
    COMMON_ADDUCT_RULES,
    annotate_mz_values,
    annotate_single_mz,
)
from .chemical_db_builders import build_chebi_annotation_table
from .chemical_db_loader import (
    DEFAULT_CHEMICAL_DB_FIELD_ALIASES,
    ChemicalLibraryEntry,
    ChemicalLibraryValidationError,
    build_chemical_library_entry,
    load_chemical_library_directory,
    load_chemical_library_table,
    summarize_chemical_library,
)
from .chemical_schema import (
    CHEMICAL_ANNOTATION_SCHEMA_VERSION,
    CHEMICAL_CONFIDENCE_LEVELS,
    CHEMICAL_SOURCE_TYPES,
    ChemicalAnnotationCandidate,
    ChemicalAnnotationRecord,
    ChemicalAnnotationValidationError,
    build_chemical_annotation_candidate,
    build_chemical_annotation_context,
    build_chemical_annotation_template,
)
from .chemical_tools import annotate_mz_with_local_library
from .conclusion_engine import build_integrated_conclusion
from .evidence_tools import build_statistical_evidence, build_unified_evidence
from .mask_tools import (
    apply_spatial_png_mask,
    binarize_mask_image,
    load_mask_png,
    normalize_mask_image,
    select_coordinates_from_mask,
    summarize_selected_coordinates,
)
from .preprocessing_recommender import recommend_preprocessing_plan
from .schemas import ToolResult
from .spatial_tools import build_spatial_preview, rasterize_coordinate_values, summarize_raster
from .validation_tools import (
    build_comparison_summary,
    build_legacy_statistical_baseline,
    build_validation_checks,
    validate_integrated_conclusion,
)

__all__ = [
    "BIOLOGICAL_METADATA_SCHEMA_VERSION",
    "REQUIRED_BIOLOGICAL_FIELDS",
    "BiologicalMetadata",
    "BiologicalMetadataValidationError",
    "build_biological_metadata_template",
    "build_biological_context",
    "resolve_biological_metadata",
    "COMMON_ADDUCT_RULES",
    "annotate_mz_values",
    "annotate_single_mz",
    "build_chebi_annotation_table",
    "DEFAULT_CHEMICAL_DB_FIELD_ALIASES",
    "ChemicalLibraryEntry",
    "ChemicalLibraryValidationError",
    "build_chemical_library_entry",
    "load_chemical_library_directory",
    "load_chemical_library_table",
    "summarize_chemical_library",
    "CHEMICAL_ANNOTATION_SCHEMA_VERSION",
    "CHEMICAL_CONFIDENCE_LEVELS",
    "CHEMICAL_SOURCE_TYPES",
    "ChemicalAnnotationCandidate",
    "ChemicalAnnotationRecord",
    "ChemicalAnnotationValidationError",
    "build_chemical_annotation_candidate",
    "build_chemical_annotation_context",
    "build_chemical_annotation_template",
    "annotate_mz_with_local_library",
    "build_integrated_conclusion",
    "build_statistical_evidence",
    "build_unified_evidence",
    "apply_spatial_png_mask",
    "load_mask_png",
    "normalize_mask_image",
    "binarize_mask_image",
    "select_coordinates_from_mask",
    "summarize_selected_coordinates",
    "recommend_preprocessing_plan",
    "build_spatial_preview",
    "rasterize_coordinate_values",
    "summarize_raster",
    "build_comparison_summary",
    "build_legacy_statistical_baseline",
    "build_validation_checks",
    "validate_integrated_conclusion",
    "ToolResult",
    "load_imzml",
    "dataset_qc_summary",
    "run_preprocessing",
    "analyze_imzml_end_to_end",
    "run_self_check",
]


def __getattr__(name: str):
    if name in {"load_imzml", "dataset_qc_summary"}:
        from .qc_tools import dataset_qc_summary, load_imzml

        return {
            "load_imzml": load_imzml,
            "dataset_qc_summary": dataset_qc_summary,
        }[name]
    if name == "run_preprocessing":
        from .preprocess_tools import run_preprocessing

        return run_preprocessing
    if name == "recommend_preprocessing_plan":
        from .preprocessing_recommender import recommend_preprocessing_plan

        return recommend_preprocessing_plan
    if name == "analyze_imzml_end_to_end":
        from .workflow_tools import analyze_imzml_end_to_end

        return analyze_imzml_end_to_end
    if name == "run_self_check":
        from .doctor import run_self_check

        return run_self_check
    if name == "annotate_mz_with_local_library":
        from .chemical_tools import annotate_mz_with_local_library

        return annotate_mz_with_local_library
    raise AttributeError(name)
