from __future__ import annotations

from typing import Any

import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.tools import get_logger

from .biological_tools import build_biological_context, resolve_biological_metadata
from .evidence_tools import build_statistical_evidence
from ._utils import build_even_sample_indices, summarize_coordinates, summarize_numeric, write_result_json
from .schemas import ToolResult

logger = get_logger("massflow.agent_tools.qc")

_META_FIELDS = (
    "name",
    "continuous",
    "processed",
    "centroid_spectrum",
    "profile_spectrum",
    "instrument_model",
    "ms1_spectrum",
    "msn_spectrum",
    "pixel_size_x",
    "pixel_size_y",
    "max_count_of_pixels_x",
    "max_count_of_pixels_y",
    "spectrum_count_num",
)


def load_imzml(
    imzml_path: str,
    *,
    target_locs: tuple[tuple[int, int], tuple[int, int]] | None = None,
    biological_metadata: dict[str, Any] | None = None,
    biological_metadata_json: str | None = None,
    biological_metadata_path: str | None = None,
    max_threads: int = 8,
    output_dir: str | None = None,
) -> ToolResult:
    """Load only the imzML header/data index and return a compact dataset summary."""
    tool_name = "load_imzml"
    warnings: list[str] = []

    try:
        resolved_biological_metadata = resolve_biological_metadata(
            biological_metadata=biological_metadata,
            biological_metadata_json=biological_metadata_json,
            biological_metadata_path=biological_metadata_path,
            imzml_path=imzml_path,
        )
        with MSDataManagerImzML(
            filepath=imzml_path,
            target_locs=target_locs,
            max_threads=max_threads,
        ) as data_manager:
            data_manager.load_head_data()
            dataset_summary = _build_dataset_summary(data_manager, warnings)
            dataset_summary["biological_context"] = build_biological_context(
                resolved_biological_metadata,
                source=_biological_metadata_source(
                    biological_metadata=biological_metadata,
                    biological_metadata_json=biological_metadata_json,
                    biological_metadata_path=biological_metadata_path,
                ),
            )
            metrics = {
                "spectra_count": dataset_summary["coordinates"]["spectra_count"],
                "shared_mz_length": dataset_summary["shared_mz_length"],
                "occupancy_ratio": dataset_summary["occupancy"]["occupancy_ratio"],
                "biological_metadata_provided": resolved_biological_metadata is not None,
            }
            result = ToolResult.success(
                tool_name,
                "Loaded imzML header and dataset index.",
                imzml_path,
                metrics=metrics,
                details=dataset_summary,
                warnings=warnings,
            )
    except Exception as exc:  # pylint: disable=broad-except
        logger.exception("load_imzml failed for %s", imzml_path)
        return ToolResult.failure(
            tool_name,
            "Failed to load imzML header.",
            imzml_path,
            errors=[str(exc)],
        )

    artifact_path = write_result_json(result, output_dir, f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def dataset_qc_summary(
    imzml_path: str,
    *,
    sample_size: int = 16,
    noise_method: str = "sd",
    top_peak_limit: int = 5,
    target_locs: tuple[tuple[int, int], tuple[int, int]] | None = None,
    biological_metadata: dict[str, Any] | None = None,
    biological_metadata_json: str | None = None,
    biological_metadata_path: str | None = None,
    max_threads: int = 8,
    output_dir: str | None = None,
) -> ToolResult:
    """Run a lightweight QC pass using evenly sampled spectra."""
    tool_name = "dataset_qc_summary"
    warnings: list[str] = []

    try:
        resolved_biological_metadata = resolve_biological_metadata(
            biological_metadata=biological_metadata,
            biological_metadata_json=biological_metadata_json,
            biological_metadata_path=biological_metadata_path,
            imzml_path=imzml_path,
        )
        with MSDataManagerImzML(
            filepath=imzml_path,
            target_locs=target_locs,
            max_threads=max_threads,
        ) as data_manager:
            data_manager.load_head_data()
            dataset_summary = _build_dataset_summary(data_manager, warnings)
            biological_context = build_biological_context(
                resolved_biological_metadata,
                source=_biological_metadata_source(
                    biological_metadata=biological_metadata,
                    biological_metadata_json=biological_metadata_json,
                    biological_metadata_path=biological_metadata_path,
                ),
            )
            dataset_summary["biological_context"] = biological_context
            sample_rows = _sample_qc_rows(
                data_manager,
                sample_size=sample_size,
                noise_method=noise_method,
                top_peak_limit=top_peak_limit,
                warnings=warnings,
            )
            qc_metrics = _build_qc_metrics(sample_rows)
            qc_flags = _build_qc_flags(dataset_summary, qc_metrics)
            statistical_evidence = build_statistical_evidence(
                dataset_summary=dataset_summary,
                qc_metrics=qc_metrics,
                qc_flags=qc_flags,
                sampled_spectra=sample_rows,
            )
            warnings.extend(flag for flag in qc_flags if flag not in warnings)

            result = ToolResult.success(
                tool_name,
                "Generated lightweight QC summary from sampled spectra.",
                imzml_path,
                metrics=qc_metrics,
                details={
                    "dataset": dataset_summary,
                    "sampled_spectra": sample_rows,
                    "qc_flags": qc_flags,
                    "statistical_evidence": statistical_evidence,
                },
                warnings=warnings,
            )
    except Exception as exc:  # pylint: disable=broad-except
        logger.exception("dataset_qc_summary failed for %s", imzml_path)
        return ToolResult.failure(
            tool_name,
            "Failed to generate QC summary.",
            imzml_path,
            errors=[str(exc)],
        )

    artifact_path = write_result_json(result, output_dir, f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def _build_dataset_summary(
    data_manager: MSDataManagerImzML,
    warnings: list[str],
) -> dict[str, Any]:
    meta = data_manager.ms.meta
    coordinates = summarize_coordinates(
        (int(spectrum.x), int(spectrum.y), int(spectrum.z)) for spectrum in data_manager.ms
    )
    occupancy = _build_occupancy_summary(data_manager, warnings)

    meta_summary = {field: getattr(meta, field, None) for field in _META_FIELDS}
    meta_summary["target_locs"] = data_manager.target_locs
    meta_summary["imzml_path"] = data_manager.imzml_filepath
    meta_summary["ibd_path"] = data_manager.ibd_filepath

    shared_mz_length = 0
    if data_manager.ms.shared_mz_list is not None:
        shared_mz_length = int(len(data_manager.ms.shared_mz_list))

    return {
        "meta": meta_summary,
        "coordinates": coordinates,
        "shared_mz_length": shared_mz_length,
        "occupancy": occupancy,
    }


def _build_occupancy_summary(
    data_manager: MSDataManagerImzML,
    warnings: list[str],
) -> dict[str, float | int | None]:
    meta = data_manager.ms.meta
    width = int(meta.max_count_of_pixels_x) if meta.max_count_of_pixels_x else 0
    height = int(meta.max_count_of_pixels_y) if meta.max_count_of_pixels_y else 0
    expected_pixels = width * height if width and height else None

    try:
        data_manager.create_ms_meta_mask()
    except Exception as exc:  # pylint: disable=broad-except
        warnings.append(f"Failed to build occupancy mask: {exc}")
        return {
            "expected_pixels": expected_pixels,
            "occupied_pixels": None,
            "occupancy_ratio": None,
        }

    mask = meta.mask
    if mask is None:
        warnings.append("Occupancy mask was not generated.")
        return {
            "expected_pixels": expected_pixels,
            "occupied_pixels": None,
            "occupancy_ratio": None,
        }

    occupied_pixels = int(np.sum(mask))
    occupancy_ratio = None
    if expected_pixels:
        occupancy_ratio = float(occupied_pixels / expected_pixels)

    return {
        "expected_pixels": expected_pixels,
        "occupied_pixels": occupied_pixels,
        "occupancy_ratio": occupancy_ratio,
    }


def _sample_qc_rows(
    data_manager: MSDataManagerImzML,
    *,
    sample_size: int,
    noise_method: str,
    top_peak_limit: int,
    warnings: list[str],
) -> list[dict[str, Any]]:
    indices = build_even_sample_indices(len(data_manager.ms), sample_size)
    sample_rows: list[dict[str, Any]] = []

    for index in indices:
        spectrum = data_manager.ms[index]
        try:
            mz_array = np.asarray(spectrum.mz_list, dtype=np.float64)
            intensity_array = np.asarray(spectrum.intensity, dtype=np.float64)
            tic = float(np.sum(intensity_array, dtype=np.float64))
            base_peak = float(np.max(intensity_array)) if intensity_array.size else 0.0
            snr_value = _safe_snr(spectrum, noise_method, warnings, index)
            sample_rows.append(
                {
                    "index": int(index),
                    "coordinate": [int(spectrum.x), int(spectrum.y), int(spectrum.z)],
                    "peak_count": int(intensity_array.size),
                    "mz_min": float(np.min(mz_array)) if mz_array.size else None,
                    "mz_max": float(np.max(mz_array)) if mz_array.size else None,
                    "tic": tic,
                    "base_peak_intensity": base_peak,
                    "snr": snr_value,
                    "top_peaks": _extract_top_peaks(
                        mz_array,
                        intensity_array,
                        top_peak_limit=top_peak_limit,
                    ),
                }
            )
        except Exception as exc:  # pylint: disable=broad-except
            warnings.append(f"Failed to inspect sampled spectrum {index}: {exc}")
        finally:
            spectrum.clear_data()

    return sample_rows


def _extract_top_peaks(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    *,
    top_peak_limit: int,
) -> list[dict[str, float]]:
    if top_peak_limit <= 0 or mz_array.size == 0 or intensity_array.size == 0:
        return []

    limit = min(int(top_peak_limit), int(intensity_array.size))
    top_indices = np.argsort(intensity_array)[-limit:][::-1]
    peaks: list[dict[str, float]] = []
    for peak_index in top_indices.tolist():
        peaks.append(
            {
                "mz": float(mz_array[peak_index]),
                "intensity": float(intensity_array[peak_index]),
            }
        )
    return peaks


def _safe_snr(
    spectrum,
    noise_method: str,
    warnings: list[str],
    sample_index: int,
) -> float | None:
    try:
        snr_value = float(SpectrumPreprocess.calculate_snr_spectrum(spectrum, method=noise_method))
    except Exception as exc:  # pylint: disable=broad-except
        warnings.append(f"Failed to calculate SNR for sampled spectrum {sample_index}: {exc}")
        return None

    if np.isfinite(snr_value):
        return snr_value

    warnings.append(f"Sampled spectrum {sample_index} produced a non-finite SNR.")
    return None


def _build_qc_metrics(sample_rows: list[dict[str, Any]]) -> dict[str, Any]:
    tic_values = [float(row["tic"]) for row in sample_rows]
    base_peak_values = [float(row["base_peak_intensity"]) for row in sample_rows]
    peak_counts = [float(row["peak_count"]) for row in sample_rows]
    snr_values = [float(row["snr"]) for row in sample_rows if row["snr"] is not None]

    return {
        "sampled_spectra_count": len(sample_rows),
        "tic_summary": summarize_numeric(tic_values),
        "base_peak_summary": summarize_numeric(base_peak_values),
        "peak_count_summary": summarize_numeric(peak_counts),
        "snr_summary": summarize_numeric(snr_values),
    }


def _build_qc_flags(
    dataset_summary: dict[str, Any],
    qc_metrics: dict[str, Any],
) -> list[str]:
    flags: list[str] = []
    occupancy_ratio = dataset_summary["occupancy"]["occupancy_ratio"]
    snr_median = qc_metrics["snr_summary"]["median"]
    peak_count_mean = qc_metrics["peak_count_summary"]["mean"]

    if occupancy_ratio is not None and occupancy_ratio < 0.1:
        flags.append("Occupancy ratio is low; verify ROI or sparse sampling assumptions.")
    if snr_median is not None and snr_median < 3.0:
        flags.append("Sampled median SNR is low; denoising or baseline correction may be needed.")
    if peak_count_mean is not None and peak_count_mean < 50:
        flags.append("Sampled spectra appear short; confirm centroid/processed data expectations.")

    return flags


def _biological_metadata_source(
    *,
    biological_metadata: dict[str, Any] | None,
    biological_metadata_json: str | None,
    biological_metadata_path: str | None,
) -> str | None:
    if biological_metadata is not None:
        return "object"
    if biological_metadata_json:
        return "json"
    if biological_metadata_path:
        return "file"
    return None
