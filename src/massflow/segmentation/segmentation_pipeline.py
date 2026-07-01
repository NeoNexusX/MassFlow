"""Segmentation pipeline and its UMAP/KMeans stage result."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
import zarr

from massflow.segmentation.region_mask import DataManager
from massflow.tools.logger import get_logger

logger = get_logger("massflow.segmentation.segmentation_pipeline")

TAB20_COLOR_NAMES = (
    "blue",
    "light_blue",
    "orange",
    "light_orange",
    "green",
    "light_green",
    "red",
    "light_red",
    "purple",
    "light_purple",
    "brown",
    "light_brown",
    "pink",
    "light_pink",
    "gray",
    "light_gray",
    "yellow_green",
    "light_yellow_green",
    "cyan",
    "light_cyan",
)


# ---------------------------------------------------------------------------
# UMAP/KMeans stage result
# ---------------------------------------------------------------------------

def _read_zarr_array(result_group: zarr.Group, name: str) -> np.ndarray:
    """Read a required Zarr array as a NumPy array."""
    if name not in result_group:
        raise ValueError(f"Zarr group must contain '{name}'.")
    array = result_group[name]
    if not isinstance(array, zarr.Array):
        raise ValueError(f"Zarr group item '{name}' must be an array.")
    return np.asarray(array)


@dataclass(frozen=True, slots=True)
class UmapKmeansResult:
    """Spatial UMAP/KMeans outputs and their human-readable label colors."""

    umap_image: np.ndarray = field(repr=False)
    kmeans_image: np.ndarray = field(repr=False)
    kmeans_label_image: np.ndarray = field(repr=False)
    label_colors: dict[int, str]

    def __post_init__(self) -> None:
        if self.kmeans_label_image.ndim != 2:
            raise ValueError("kmeans_label_image must be a two-dimensional array.")
        expected_image_shape = (*self.kmeans_label_image.shape, 3)
        if self.umap_image.shape != expected_image_shape:
            raise ValueError("umap_image must have shape (height, width, 3).")
        if self.kmeans_image.shape != expected_image_shape:
            raise ValueError("kmeans_image must have shape (height, width, 3).")

    @classmethod
    def from_zarr(cls, zarr_path: str | Path) -> UmapKmeansResult:
        """Load a feature-analysis result from a saved Zarr group."""
        result_group = zarr.open_group(str(zarr_path), mode="r")
        label_image = _read_zarr_array(result_group, "kmeans_label_image")
        labels = np.unique(label_image[label_image >= 0])
        saved_colors = result_group.attrs.get("label_colors", {})
        if saved_colors:
            if not isinstance(saved_colors, dict):
                raise ValueError("Zarr attribute 'label_colors' must be an object.")
            label_colors = {
                int(label): str(color_name)
                for label, color_name in saved_colors.items()
            }
        else:
            if labels.size and int(labels[-1]) >= len(TAB20_COLOR_NAMES):
                raise ValueError(
                    f"KMeans labels must be in the range [0, {len(TAB20_COLOR_NAMES) - 1}] "
                    "when Zarr attribute 'label_colors' is missing."
                )
            label_colors = {
                int(label): TAB20_COLOR_NAMES[int(label)]
                for label in labels
            }
        return cls(
            umap_image=_read_zarr_array(result_group, "umap_image"),
            kmeans_image=_read_zarr_array(result_group, "kmeans_image"),
            kmeans_label_image=label_image,
            label_colors=label_colors,
        )

    @classmethod
    def from_npz(cls, npz_path: str | Path) -> UmapKmeansResult:
        """Load a feature-analysis result from a saved NPZ file."""
        with np.load(npz_path, allow_pickle=False) as result_file:
            label_colors = {
                int(label): str(color_name)
                for label, color_name in zip(
                    result_file["label_color_labels"],
                    result_file["label_color_names"],
                )
            }
            return cls(
                umap_image=result_file["umap_image"],
                kmeans_image=result_file["kmeans_image"],
                kmeans_label_image=result_file["kmeans_label_image"],
                label_colors=label_colors,
            )

    def print_label_colors(self) -> None:
        """Log the KMeans label-to-color mapping."""
        mapping = ", ".join(
            f"{label}: {color_name}"
            for label, color_name in sorted(self.label_colors.items())
        )
        logger.info(f"KMeans label colors: {mapping}")


# ---------------------------------------------------------------------------
# Segmentation workflow
# ---------------------------------------------------------------------------

class SegmentationPipeline:
    """Run and retain the state of a complete segmentation workflow."""

    def __init__(
        self,
        data_manager: DataManager,
        *,
        output_dir: str | Path = "./segmentation_output",
        result: UmapKmeansResult | None = None,
    ) -> None:
        self.data_manager = data_manager
        self.output_dir = Path(output_dir)
        self.result = result
        self.mask: np.ndarray | None = None

    def _require_result(self) -> UmapKmeansResult:
        if self.result is None:
            raise RuntimeError("run_umap_kmeans() must be called before cluster segmentation.")
        return self.result

    def _require_mask(self) -> np.ndarray:
        if self.mask is None:
            raise RuntimeError(
                "A cluster segmentation method must be called before apply_region_mask()."
            )
        return self.mask

    def run_umap_kmeans(
        self,
        *,
        sample_ratio: float = 0.1,
        transform_batch_size: int | None = 1024,
        dpi: int = 600,
        save_images: bool = True,
        save_matrix: Literal["npz", "zarr"] | None = None,
        n_components: int = 3,
        umap_metric: str = "cosine",
        n_clusters: int = 10,
        random_state: int | None = None,
    ) -> UmapKmeansResult:
        """Run UMAP/KMeans, retain its result, and clear any previous mask."""
        from massflow.segmentation.umap_kmeans import plot_umap_kmeans_image

        self.result = plot_umap_kmeans_image(
            self.data_manager,
            output_dir=self.output_dir / "umap_kmeans",
            sample_ratio=sample_ratio,
            transform_batch_size=transform_batch_size,
            dpi=dpi,
            save_images=save_images,
            save_matrix=save_matrix,
            n_components=n_components,
            umap_metric=umap_metric,
            n_clusters=n_clusters,
            random_state=random_state,
        )
        self.mask = None
        return self.result

    def remove_clusters(self, labels_to_remove: Sequence[int]) -> np.ndarray:
        """Replace the current mask by removing the specified clusters."""
        from massflow.segmentation.cluster_segmentation import remove_clusters

        self.mask = remove_clusters(
            self._require_result(),
            labels_to_remove,
            output_dir=self.output_dir / "cluster_segmentation",
        )
        return self.mask

    def remove_background(
        self,
        labels_to_remove: Sequence[int] | None = None,
        *,
        remove_num: int = 1,
    ) -> np.ndarray:
        """Replace the current mask by removing exterior cluster layers."""
        from massflow.segmentation.cluster_segmentation import remove_background

        self.mask = remove_background(
            self._require_result(),
            labels_to_remove,
            remove_num=remove_num,
            output_dir=self.output_dir / "cluster_segmentation",
        )
        return self.mask

    def keep_clusters(self, labels_to_keep: Sequence[int]) -> np.ndarray:
        """Replace the current mask by retaining the specified clusters."""
        from massflow.segmentation.cluster_segmentation import keep_clusters

        self.mask = keep_clusters(
            self._require_result(),
            labels_to_keep,
            output_dir=self.output_dir / "cluster_segmentation",
        )
        return self.mask

    def apply_region_mask(
        self,
        output_path: str | Path,
        *,
        batch_size: int = 256,
        temp_dir: str | Path | None = None,
    ) -> DataManager:
        """Apply the current mask to the source data manager and return the new manager."""
        from massflow.segmentation.region_mask import apply_region_mask

        return apply_region_mask(
            self.data_manager,
            self._require_mask(),
            output_path,
            batch_size=batch_size,
            temp_dir=temp_dir,
        )
