from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.tools import get_logger

from ._utils import ensure_directory, write_result_json
from .schemas import ToolResult

logger = get_logger("massflow.agent_tools.analysis")

MSIPL_ENABLED_DEFAULT = False


def run_msipl_analysis(
    imzml_path: str,
    *,
    mode: str = "train",
    weights_path: str | None = None,
    quick_train: bool = True,
    enable_msipl: bool = MSIPL_ENABLED_DEFAULT,
    output_dir: str | None = None,
    config_overrides: dict[str, Any] | None = None,
) -> ToolResult:
    """Run msiPL train/infer and return an agent-friendly summary."""
    tool_name = "run_msipl_analysis"
    warnings: list[str] = []
    output_root = ensure_directory(output_dir, temp_prefix="massflow_msipl")

    # Keep the msiPL-facing tool reserved until the local training/inference
    # stack is ready. This avoids exposing a half-configured capability to Agent.
    if not enable_msipl:
        result = ToolResult.failure(
            tool_name,
            "msiPL analysis is currently disabled in agent_tools. Enable it after the local msiPL pipeline is ready.",
            imzml_path,
            warnings=["msiPL stage is reserved and not enabled by default."],
            errors=["Set enable_msipl=True after msiPL is fully configured."],
        )
        artifact_path = write_result_json(result, str(output_root), f"{tool_name}_{result.job_id}.json")
        if artifact_path:
            result.artifacts["summary_json"] = artifact_path
        return result

    try:
        from massflow.dl.msipl import MsiplConfig, infer_msipl, train_msipl
    except Exception as exc:  # pylint: disable=broad-except
        logger.exception("Failed to import msiPL modules.")
        return ToolResult.failure(
            tool_name,
            "msiPL modules are unavailable. Install the optional msiPL dependencies first.",
            imzml_path,
            errors=[str(exc)],
        )

    normalized_mode = str(mode).strip().lower()
    if normalized_mode not in {"train", "infer"}:
        return ToolResult.failure(
            tool_name,
            "Unsupported msiPL mode.",
            imzml_path,
            errors=[f"mode must be 'train' or 'infer', got {mode!r}"],
        )

    if normalized_mode == "infer" and not weights_path:
        return ToolResult.failure(
            tool_name,
            "Inference mode requires a weights file.",
            imzml_path,
            errors=["weights_path is required when mode='infer'."],
        )

    cfg = _build_msipl_config(
        MsiplConfig,
        mode=normalized_mode,
        quick_train=quick_train,
        output_dir=str(output_root),
        weights_path=weights_path,
        config_overrides=config_overrides,
    )

    try:
        with MSDataManagerImzML(filepath=imzml_path) as data_manager:
            data_manager.load_head_data()
            if normalized_mode == "train":
                raw_result = train_msipl(data_manager, cfg)
            else:
                raw_result = infer_msipl(data_manager, weights_path=weights_path or "", config=cfg)
    except Exception as exc:  # pylint: disable=broad-except
        logger.exception("run_msipl_analysis failed for %s", imzml_path)
        return ToolResult.failure(
            tool_name,
            "Failed to execute msiPL analysis.",
            imzml_path,
            errors=[str(exc)],
        )

    metrics = _build_msipl_metrics(raw_result)
    details = {
        "mode": normalized_mode,
        "config": dict(cfg.__dict__),
        "result_summary": _build_msipl_result_summary(raw_result),
    }
    result = ToolResult.success(
        tool_name,
        "msiPL analysis completed.",
        imzml_path,
        artifacts=dict(raw_result.artifacts),
        metrics=metrics,
        details=details,
        warnings=warnings,
    )

    artifact_path = write_result_json(result, str(output_root), f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def _build_msipl_config(
    config_cls,
    *,
    mode: str,
    quick_train: bool,
    output_dir: str,
    weights_path: str | None,
    config_overrides: dict[str, Any] | None,
):
    base_kwargs: dict[str, Any] = {
        "backend": "lightning",
        "use_gpu": False,
        "show_progress": False,
        "quick_train": quick_train,
        "quick_epochs": 10,
        "output_dir": output_dir,
        "export_arrays": True,
        "export_summary": True,
        "export_report": True,
        "run_clustering": True,
        "run_peak_learning": True,
        "report_cluster_k": 5,
        "report_clusters": 5,
        "report_peaks_per_cluster": 2,
        "save_weights": mode == "train",
    }
    if mode == "infer":
        base_kwargs["save_weights"] = False
    if weights_path:
        base_kwargs["weights_path"] = weights_path
    if config_overrides:
        base_kwargs.update(config_overrides)
    return config_cls(**base_kwargs)


def _build_msipl_metrics(raw_result) -> dict[str, Any]:
    latent_shape = list(raw_result.latent_z.shape) if raw_result.latent_z is not None else None
    reconstructed_shape = (
        list(raw_result.reconstructed.shape) if raw_result.reconstructed is not None else None
    )
    cluster_count = 0
    if raw_result.cluster_labels is not None and raw_result.cluster_labels.size > 0:
        cluster_count = int(len(np.unique(raw_result.cluster_labels)))

    learned_peak_count = 0
    if raw_result.learned_mz_peaks is not None:
        learned_peak_count = int(np.size(raw_result.learned_mz_peaks))

    return {
        "mse": raw_result.mse,
        "latent_shape": latent_shape,
        "reconstructed_shape": reconstructed_shape,
        "cluster_count": cluster_count,
        "learned_peak_count": learned_peak_count,
    }


def _build_msipl_result_summary(raw_result) -> dict[str, Any]:
    history_keys = sorted(list(getattr(raw_result, "history", {}).keys()))
    summary = {
        "output_dir": raw_result.output_dir,
        "backend": raw_result.backend,
        "latent_dim": raw_result.latent_dim,
        "metadata": dict(raw_result.metadata),
        "weights_path": getattr(raw_result, "weights_path", None),
        "history_keys": history_keys,
    }
    if raw_result.cluster_labels is not None and raw_result.cluster_labels.size > 0:
        labels = raw_result.cluster_labels.astype(int)
        summary["cluster_histogram"] = {
            str(label): int(np.sum(labels == label)) for label in np.unique(labels)
        }
    if raw_result.learned_mz_peaks is not None:
        summary["learned_mz_peaks_shape"] = list(np.shape(raw_result.learned_mz_peaks))
    return summary
