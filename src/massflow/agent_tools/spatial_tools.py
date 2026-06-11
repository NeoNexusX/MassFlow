from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.tools import get_logger

from ._utils import ensure_directory, write_result_json
from .schemas import ToolResult

logger = get_logger("massflow.agent_tools.spatial")


def build_spatial_preview(
    imzml_path: str,
    *,
    max_threads: int = 8,
    output_dir: str | None = None,
    colormap: str = "magma",
) -> ToolResult:
    """Generate basic spatial preview images from one imzML dataset."""

    tool_name = "build_spatial_preview"
    output_root = ensure_directory(output_dir, temp_prefix="massflow_spatial")
    warnings: list[str] = []

    try:
        with MSDataManagerImzML(filepath=imzml_path, max_threads=max_threads) as data_manager:
            data_manager.load_head_data()
            tic_rows: list[tuple[int, int, float]] = []
            base_peak_rows: list[tuple[int, int, float]] = []

            for spectrum in data_manager.ms:
                intensity = np.asarray(spectrum.intensity, dtype=np.float64)
                tic_value = float(np.sum(intensity, dtype=np.float64)) if intensity.size else 0.0
                base_peak_value = float(np.max(intensity)) if intensity.size else 0.0
                tic_rows.append((int(spectrum.x), int(spectrum.y), tic_value))
                base_peak_rows.append((int(spectrum.x), int(spectrum.y), base_peak_value))
                spectrum.clear_data()

        tic_image, image_meta = rasterize_coordinate_values(tic_rows)
        base_peak_image, _ = rasterize_coordinate_values(base_peak_rows)

        tic_png_path = str(output_root / "tic_preview.png")
        base_peak_png_path = str(output_root / "base_peak_preview.png")
        np.save(output_root / "tic_preview.npy", tic_image)
        np.save(output_root / "base_peak_preview.npy", base_peak_image)

        _save_image(tic_image, tic_png_path, colormap=colormap)
        _save_image(base_peak_image, base_peak_png_path, colormap=colormap)

        result = ToolResult.success(
            tool_name,
            "Generated basic spatial preview images.",
            imzml_path,
            artifacts={
                "tic_png": tic_png_path,
                "base_peak_png": base_peak_png_path,
                "tic_npy": str(output_root / "tic_preview.npy"),
                "base_peak_npy": str(output_root / "base_peak_preview.npy"),
            },
            metrics={
                "pixel_count": image_meta["pixel_count"],
                "image_width": image_meta["width"],
                "image_height": image_meta["height"],
                "tic_max": image_meta["max_value"],
                "tic_nonzero_pixels": image_meta["nonzero_pixels"],
            },
            details={
                "image_meta": image_meta,
                "base_peak_summary": summarize_raster(base_peak_image),
                "tic_summary": summarize_raster(tic_image),
                "recommended_next_action": (
                    "Use tic_preview.png or base_peak_preview.png to manually remove matrix/background "
                    "and then upload a PNG mask for region-restricted analysis."
                ),
            },
            warnings=warnings,
        )
    except Exception as exc:  # pylint: disable=broad-except
        logger.exception("build_spatial_preview failed for %s", imzml_path)
        return ToolResult.failure(
            tool_name,
            "Failed to generate spatial preview images.",
            imzml_path,
            errors=[str(exc)],
        )

    artifact_path = write_result_json(result, str(output_root), f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def rasterize_coordinate_values(
    rows: list[tuple[int, int, float]],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Rasterize x/y/value triples into an image array."""

    if not rows:
        empty_image = np.zeros((0, 0), dtype=np.float32)
        return empty_image, {
            "width": 0,
            "height": 0,
            "min_x": None,
            "max_x": None,
            "min_y": None,
            "max_y": None,
            "pixel_count": 0,
            "nonzero_pixels": 0,
            "max_value": None,
        }

    xs = np.asarray([row[0] for row in rows], dtype=np.int32)
    ys = np.asarray([row[1] for row in rows], dtype=np.int32)
    values = np.asarray([row[2] for row in rows], dtype=np.float32)

    min_x = int(np.min(xs))
    max_x = int(np.max(xs))
    min_y = int(np.min(ys))
    max_y = int(np.max(ys))
    width = max_x - min_x + 1
    height = max_y - min_y + 1

    image = np.full((height, width), np.nan, dtype=np.float32)
    for x_value, y_value, pixel_value in rows:
        image[int(y_value) - min_y, int(x_value) - min_x] = float(pixel_value)

    return image, {
        "width": width,
        "height": height,
        "min_x": min_x,
        "max_x": max_x,
        "min_y": min_y,
        "max_y": max_y,
        "pixel_count": int(len(rows)),
        "nonzero_pixels": int(np.sum(np.nan_to_num(image, nan=0.0) > 0)),
        "max_value": float(np.nanmax(image)) if not np.isnan(image).all() else None,
    }


def summarize_raster(image: np.ndarray) -> dict[str, Any]:
    if image.size == 0:
        return {"count": 0, "min": None, "max": None, "mean": None}

    valid_values = image[~np.isnan(image)]
    if valid_values.size == 0:
        return {"count": 0, "min": None, "max": None, "mean": None}

    return {
        "count": int(valid_values.size),
        "min": float(np.min(valid_values)),
        "max": float(np.max(valid_values)),
        "mean": float(np.mean(valid_values)),
    }


def _save_image(image: np.ndarray, output_path: str, *, colormap: str) -> None:
    plot_image = np.nan_to_num(image, nan=0.0)
    plt.imsave(output_path, plot_image, cmap=colormap)
