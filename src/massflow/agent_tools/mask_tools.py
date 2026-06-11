from __future__ import annotations

from pathlib import Path
from typing import Any
import json

import matplotlib.pyplot as plt
import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.tools import get_logger

from ._utils import ensure_directory, write_result_json
from .schemas import ToolResult
from .spatial_tools import rasterize_coordinate_values

logger = get_logger("massflow.agent_tools.mask")


def apply_spatial_png_mask(
    imzml_path: str,
    mask_png_path: str,
    *,
    max_threads: int = 8,
    output_dir: str | None = None,
    threshold: float = 0.5,
    keep_mode: str = "white_foreground",
) -> ToolResult:
    """Map a PNG mask onto imzML coordinates and return the retained pixels."""

    tool_name = "apply_spatial_png_mask"
    output_root = ensure_directory(output_dir, temp_prefix="massflow_mask")

    try:
        with MSDataManagerImzML(filepath=imzml_path, max_threads=max_threads) as data_manager:
            data_manager.load_head_data()
            coordinate_rows = [(int(spectrum.x), int(spectrum.y), 1.0) for spectrum in data_manager.ms]

        _, image_meta = rasterize_coordinate_values(coordinate_rows)
        mask_image = load_mask_png(mask_png_path)
        binary_mask = binarize_mask_image(mask_image, threshold=threshold, keep_mode=keep_mode)
        _validate_mask_shape(binary_mask, image_meta)

        selected_coordinates = select_coordinates_from_mask(coordinate_rows, binary_mask, image_meta)
        selection_summary = summarize_selected_coordinates(selected_coordinates)
        selected_ratio = (
            float(len(selected_coordinates) / len(coordinate_rows)) if coordinate_rows else 0.0
        )

        binary_png_path = str(output_root / "mask_binary.png")
        npy_path = str(output_root / "mask_binary.npy")
        coordinates_json_path = str(output_root / "selected_coordinates.json")
        bbox_json_path = str(output_root / "selected_bbox.json")

        plt.imsave(binary_png_path, binary_mask.astype(np.float32), cmap="gray")
        np.save(npy_path, binary_mask.astype(np.uint8))
        Path(coordinates_json_path).write_text(
            json.dumps({"selected_coordinates": selected_coordinates}, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        Path(bbox_json_path).write_text(
            json.dumps(selection_summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        result = ToolResult.success(
            tool_name,
            "Applied PNG mask to dataset coordinates.",
            imzml_path,
            artifacts={
                "mask_binary_png": binary_png_path,
                "mask_binary_npy": npy_path,
                "selected_coordinates_json": coordinates_json_path,
                "selected_bbox_json": bbox_json_path,
            },
            metrics={
                "total_coordinate_count": len(coordinate_rows),
                "selected_coordinate_count": len(selected_coordinates),
                "selected_ratio": selected_ratio,
                "image_width": image_meta["width"],
                "image_height": image_meta["height"],
            },
            details={
                "mask_png_path": mask_png_path,
                "image_meta": image_meta,
                "selection_summary": selection_summary,
                "selected_coordinates_preview": selected_coordinates[:20],
                "recommended_next_action": (
                    "Use selected_coordinates_json as the ROI input for downstream key m/z extraction "
                    "and region-restricted chemical annotation."
                ),
            },
        )
    except Exception as exc:  # pylint: disable=broad-except
        logger.exception("apply_spatial_png_mask failed for %s", imzml_path)
        return ToolResult.failure(
            tool_name,
            "Failed to apply the PNG mask to the dataset.",
            imzml_path,
            errors=[str(exc)],
        )

    artifact_path = write_result_json(result, str(output_root), f"{tool_name}_{result.job_id}.json")
    if artifact_path:
        result.artifacts["summary_json"] = artifact_path
    return result


def load_mask_png(mask_png_path: str) -> np.ndarray:
    image = plt.imread(mask_png_path)
    return normalize_mask_image(image)


def normalize_mask_image(mask_image: np.ndarray) -> np.ndarray:
    array = np.asarray(mask_image)
    if array.ndim == 2:
        grayscale = array.astype(np.float32)
    elif array.ndim == 3 and array.shape[2] >= 3:
        grayscale = np.mean(array[..., :3].astype(np.float32), axis=2)
    else:
        raise ValueError("Mask PNG must be a grayscale or RGB/RGBA image.")

    if grayscale.size == 0:
        raise ValueError("Mask PNG is empty.")

    max_value = float(np.max(grayscale))
    if max_value > 1.0:
        grayscale = grayscale / 255.0
    return grayscale


def binarize_mask_image(
    mask_image: np.ndarray,
    *,
    threshold: float,
    keep_mode: str,
) -> np.ndarray:
    normalized_mode = keep_mode.strip().lower()
    if normalized_mode not in {"white_foreground", "black_foreground"}:
        raise ValueError("keep_mode must be 'white_foreground' or 'black_foreground'.")

    if normalized_mode == "white_foreground":
        return mask_image >= float(threshold)
    return mask_image <= float(threshold)


def select_coordinates_from_mask(
    coordinate_rows: list[tuple[int, int, float]],
    binary_mask: np.ndarray,
    image_meta: dict[str, Any],
) -> list[list[int]]:
    min_x = int(image_meta["min_x"])
    min_y = int(image_meta["min_y"])
    selected: list[list[int]] = []
    for x_value, y_value, _ in coordinate_rows:
        row_index = int(y_value) - min_y
        col_index = int(x_value) - min_x
        if bool(binary_mask[row_index, col_index]):
            selected.append([int(x_value), int(y_value), 1])
    return selected


def summarize_selected_coordinates(selected_coordinates: list[list[int]]) -> dict[str, Any]:
    if not selected_coordinates:
        return {
            "selected_count": 0,
            "min_x": None,
            "max_x": None,
            "min_y": None,
            "max_y": None,
        }

    xs = np.asarray([row[0] for row in selected_coordinates], dtype=np.int32)
    ys = np.asarray([row[1] for row in selected_coordinates], dtype=np.int32)
    return {
        "selected_count": int(len(selected_coordinates)),
        "min_x": int(np.min(xs)),
        "max_x": int(np.max(xs)),
        "min_y": int(np.min(ys)),
        "max_y": int(np.max(ys)),
    }


def _validate_mask_shape(binary_mask: np.ndarray, image_meta: dict[str, Any]) -> None:
    expected_shape = (int(image_meta["height"]), int(image_meta["width"]))
    if tuple(binary_mask.shape) != expected_shape:
        raise ValueError(
            "Mask PNG dimensions do not match the dataset raster size: "
            f"{tuple(binary_mask.shape)} != {expected_shape}."
        )
