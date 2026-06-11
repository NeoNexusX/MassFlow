"""Save presentation-ready TIC, UMAP, and KMeans images."""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from massflow.data_manager import MSDataManagerImzML
from massflow.segmentation.umap_kmeans import FeatureAnalysisResult
from massflow.tools.logger import get_logger
from massflow.tools.plot_ion_image import read_tic_image

logger = get_logger("massflow.segmentation.plot_image")


# ---------------------------------------------------------------------------
# Image transformation helpers
# ---------------------------------------------------------------------------


def _scale_image(image: np.ndarray, output_scale: int) -> np.ndarray:
    """Scale an image by repeating pixels without interpolation."""
    return np.repeat(np.repeat(image, output_scale, axis=0), output_scale, axis=1)


def _map_umap_colors(
    image: np.ndarray,
    gamma: float | tuple[float, float, float],
) -> np.ndarray:
    """Apply independent gamma mappings to the UMAP RGB channels."""
    gamma_array = np.asarray(gamma, dtype=np.float32)
    if gamma_array.ndim == 0:
        gamma_array = np.repeat(gamma_array, 3)
    if gamma_array.shape != (3,) or np.any(gamma_array <= 0):
        raise ValueError(
            "umap_gamma must be positive or contain three positive values."
        )

    normalized = image.astype(np.float32) / 255.0
    mapped = np.power(normalized, gamma_array.reshape(1, 1, 3))
    return np.clip(mapped * 255, 0, 255).astype(np.uint8)


def _colorize_tic(image: np.ndarray, cmap_name: str) -> np.ndarray:
    """Map TIC intensities to an RGB image using positive-value contrast."""
    positive = image[np.isfinite(image) & (image > 0)]
    vmax = float(np.percentile(positive, 95)) if positive.size else 1.0
    normalized = np.clip(image / vmax, 0.0, 1.0)
    rgba = matplotlib.colormaps[cmap_name](normalized, bytes=True)
    return np.asarray(rgba[..., :3], dtype=np.uint8)


def _save_png(image: np.ndarray, output_path: Path, *, dpi: int) -> None:
    """Save an RGB image as PNG."""
    plt.imsave(output_path, image, dpi=dpi)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def plot_front_images(
    data_manager: MSDataManagerImzML,
    source: FeatureAnalysisResult,
    output_dir: str | Path,
    dpi: int = 300,
    tic_cmaps: str = "viridis",
    output_scale: int = 2,
    umap_gamma: float | tuple[float, float, float] = 0.8,
) -> None:
    """Save scaled TIC, UMAP, and KMeans PNG images.

    ``umap_gamma`` accepts one gamma value for all RGB channels or three values
    for independent red, green, and blue mappings. Values below 1 brighten a
    channel; values above 1 darken it.
    """
    if dpi <= 0:
        raise ValueError("dpi must be positive.")
    if output_scale < 1:
        raise ValueError("output_scale must be a positive integer.")

    expected_shape = (
        int(data_manager.ms.meta.max_count_of_pixels_y), # type: ignore
        int(data_manager.ms.meta.max_count_of_pixels_x), # type: ignore
    )
    if source.kmeans_label_image.shape != expected_shape:
        raise ValueError(
            "Feature-analysis result and data manager must have matching dimensions."
        )

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    tic_image = _colorize_tic(read_tic_image(data_manager), tic_cmaps)
    umap_image = _map_umap_colors(source.umap_image, umap_gamma)

    _save_png(
        _scale_image(tic_image, output_scale),
        out / "tic_image.png",
        dpi=dpi,
    )
    _save_png(
        _scale_image(umap_image, output_scale),
        out / "umap_image.png",
        dpi=dpi,
    )
    _save_png(
        _scale_image(source.kmeans_image, output_scale),
        out / "kmeans_image.png",
        dpi=dpi,
    )

    logger.info(f"Front images saved to: {out}")


if __name__ == "__main__":
    dm = MSDataManagerImzML(
        filepath="/Users/dre/Desktop/dre/MassFlow/tests/brain-1-imzml/segmentation_masked.imzML"
    )
    dm.load_head_data()

    uk_result = FeatureAnalysisResult.from_zarr(
        "/Users/dre/Desktop/dre/MassFlow/tests/output/umap_kmeans_result.zarr"
    )

    plot_front_images(
        data_manager=dm,
        source=uk_result,
        output_dir="/Users/dre/Desktop/dre/MassFlow/tests/front_images",
        dpi=300,
        tic_cmaps="magma",
        output_scale=2,
        umap_gamma=(1, 1, 1),
    )
