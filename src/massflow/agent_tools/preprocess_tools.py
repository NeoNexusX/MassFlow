from __future__ import annotations

from pathlib import Path
from typing import Any

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor
from massflow.tools import get_logger

from ._utils import copy_artifact, ensure_directory, write_result_json
from .qc_tools import _build_dataset_summary, _build_qc_metrics, _sample_qc_rows
from .schemas import ToolResult

logger = get_logger("massflow.agent_tools.preprocess")

DEFAULT_PREPROCESS_OPERATIONS: list[dict[str, Any]] = [
    {"name": "baseline_correction", "params": {"method": "locmin_numba"}},
    {"name": "noise_reduction", "params": {"method": "gaussian_numba"}},
    {"name": "normalization", "params": {"method": "tic_numba"}},
]

_SUPPORTED_PREPROCESS_OPERATIONS = {
    "baseline_correction",
    "noise_reduction",
    "normalization",
    "peak_pick",
    "peak_align",
}


def run_preprocessing(
    imzml_path: str,
    *,
    operations: list[dict[str, Any]] | None = None,
    batch_size: int = 256,
    max_threads: int = 8,
    keep_order: bool = False,
    noise_method: str = "sd",
    qc_sample_size: int = 12,
    output_dir: str | None = None,
) -> ToolResult:
    """Execute a preprocessing pipeline and persist the processed imzML pair."""
    tool_name = "run_preprocessing"
    warnings: list[str] = []

    normalized_operations = _normalize_operations(operations)
    if not normalized_operations:
        return ToolResult.failure(
            tool_name,
            "No preprocessing operations were provided.",
            imzml_path,
            errors=["At least one preprocessing operation is required."],
        )

    output_root = ensure_directory(output_dir, temp_prefix="massflow_preprocess")
    processed_data_manager = None

    try:
        with MSDataManagerImzML(filepath=imzml_path, max_threads=max_threads) as data_manager:
            data_manager.load_head_data()
            before_dataset = _build_dataset_summary(data_manager, warnings)
            before_samples = _sample_qc_rows(
                data_manager,
                sample_size=qc_sample_size,
                noise_method=noise_method,
                top_peak_limit=5,
                warnings=warnings,
            )
            before_qc = _build_qc_metrics(before_samples)

            preprocessor = Preprocessor(
                data_manager,
                batch_size=batch_size,
                temp_dir=str(output_root / "staging"),
                keep_order=keep_order,
            )
            _apply_operations(preprocessor, normalized_operations)
            processed_data_manager = preprocessor.start()

        if processed_data_manager is None:
            raise RuntimeError("Preprocessor returned no output data manager.")

        try:
            after_dataset = _build_dataset_summary(processed_data_manager, warnings)
            after_samples = _sample_qc_rows(
                processed_data_manager,
                sample_size=qc_sample_size,
                noise_method=noise_method,
                top_peak_limit=5,
                warnings=warnings,
            )
            after_qc = _build_qc_metrics(after_samples)
            persisted_artifacts = _persist_processed_dataset(
                processed_data_manager,
                imzml_path=imzml_path,
                destination_dir=output_root / "processed",
            )
        finally:
            processed_data_manager.close()

        comparison = _build_preprocess_comparison(before_qc, after_qc)
        metrics = {
            "operations_count": len(normalized_operations),
            "output_spectra_count": after_dataset["coordinates"]["spectra_count"],
            "snr_median_delta": comparison["snr_median_delta"],
            "tic_median_delta": comparison["tic_median_delta"],
        }
        result = ToolResult.success(
            tool_name,
            "Preprocessing pipeline completed and processed imzML was persisted.",
            imzml_path,
            artifacts=persisted_artifacts,
            metrics=metrics,
            details={
                "operations": normalized_operations,
                "before": {
                    "dataset": before_dataset,
                    "qc": before_qc,
                    "sampled_spectra": before_samples,
                },
                "after": {
                    "dataset": after_dataset,
                    "qc": after_qc,
                    "sampled_spectra": after_samples,
                },
                "comparison": comparison,
            },
            warnings=warnings,
        )
    except Exception as exc:  # pylint: disable=broad-except
        logger.exception("run_preprocessing failed for %s", imzml_path)
        return ToolResult.failure(
            tool_name,
            "Failed to execute preprocessing pipeline.",
            imzml_path,
            errors=[str(exc)],
            warnings=warnings,
        )

    artifact_path = write_result_json(result, str(output_root), f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def _normalize_operations(
    operations: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    source = operations or DEFAULT_PREPROCESS_OPERATIONS
    normalized: list[dict[str, Any]] = []

    for item in source:
        name = str(item.get("name", "")).strip()
        params = item.get("params", {})
        if name not in _SUPPORTED_PREPROCESS_OPERATIONS:
            raise ValueError(f"Unsupported preprocessing operation: {name}")
        if params is None:
            params = {}
        if not isinstance(params, dict):
            raise TypeError(f"Operation params for {name} must be a dict.")
        normalized.append({"name": name, "params": dict(params)})

    return normalized


def _apply_operations(preprocessor: Preprocessor, operations: list[dict[str, Any]]) -> None:
    for operation in operations:
        name = operation["name"]
        params = operation["params"]
        method = getattr(preprocessor, name, None)
        if method is None:
            raise AttributeError(f"Preprocessor does not expose operation {name}")
        method(**params)


def _persist_processed_dataset(
    processed_data_manager: MSDataManagerImzML,
    *,
    imzml_path: str,
    destination_dir: Path,
) -> dict[str, str]:
    dataset_name = Path(imzml_path).stem
    imzml_copy = copy_artifact(
        processed_data_manager.imzml_filepath,
        destination_dir / f"{dataset_name}_processed.imzML",
    )
    ibd_copy = copy_artifact(
        processed_data_manager.ibd_filepath,
        destination_dir / f"{dataset_name}_processed.ibd",
    )

    artifacts: dict[str, str] = {}
    if imzml_copy:
        artifacts["processed_imzml"] = imzml_copy
    if ibd_copy:
        artifacts["processed_ibd"] = ibd_copy
    return artifacts


def _build_preprocess_comparison(
    before_qc: dict[str, Any],
    after_qc: dict[str, Any],
) -> dict[str, Any]:
    return {
        "snr_median_before": before_qc["snr_summary"]["median"],
        "snr_median_after": after_qc["snr_summary"]["median"],
        "snr_median_delta": _safe_delta(
            before_qc["snr_summary"]["median"],
            after_qc["snr_summary"]["median"],
        ),
        "tic_median_before": before_qc["tic_summary"]["median"],
        "tic_median_after": after_qc["tic_summary"]["median"],
        "tic_median_delta": _safe_delta(
            before_qc["tic_summary"]["median"],
            after_qc["tic_summary"]["median"],
        ),
        "peak_count_mean_before": before_qc["peak_count_summary"]["mean"],
        "peak_count_mean_after": after_qc["peak_count_summary"]["mean"],
        "peak_count_mean_delta": _safe_delta(
            before_qc["peak_count_summary"]["mean"],
            after_qc["peak_count_summary"]["mean"],
        ),
    }


def _safe_delta(before_value: Any, after_value: Any) -> float | None:
    if before_value is None or after_value is None:
        return None
    return float(after_value) - float(before_value)
