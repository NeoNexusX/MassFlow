# UMAP Dimensionality Reduction and Visualization

This page introduces the role of UMAP in MassFlow and how to interpret its results. UMAP is used for dimensionality reduction and visualization of high-dimensional mass-spectrometry features.



Authors: Leqiang Dong, Jingwei Dong

**Feedback:** If you find any problems while using MassFlow or reading the documentation, you can submit feedback through [GitHub Issues](https://github.com/NeoNexusX/MassFlow/issues), or contact us by email at 707244007@qq.com.

## Overview

In an MSI dataset, each pixel contains a large number of *m/z* positions and corresponding intensities. Directly comparing these high-dimensional mass-spectrometry features is difficult.

UMAP compresses the high-dimensional spectral features of each pixel into a smaller number of feature dimensions, while trying to keep pixels with similar overall spectral characteristics close to one another in the low-dimensional space. MassFlow maps three UMAP features to the RGB color channels and places them back at the pixels' spatial coordinates to form a spatial UMAP image.

The UMAP image shows whether the overall mass-spectrometry characteristics of different pixels are similar and how those differences are distributed in tissue space.

UMAP itself does not assign biological names to regions, and a particular color cannot directly be interpreted as tumor, necrotic, or normal tissue.

## Current API

The current workflow calculates UMAP and KMeans together through `SegmentationPipeline.run_umap_kmeans()`. Although the API returns both results at the same time, they serve different purposes:

- UMAP: dimensionality reduction and visualization of continuous features;
- KMeans: discrete class assignment based on the low-dimensional features.

### SegmentationPipeline.run_umap_kmeans

```python
massflow.segmentation.SegmentationPipeline.run_umap_kmeans(
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
    show: bool = False,
) -> UmapKmeansResult
```

- Function: runs UMAP and KMeans using the mass-spectrometry features of each pixel and generates spatial clustering results.

- Parameters:

  - `sample_ratio`: proportion of pixels used to build the UMAP model, default `0.1`. Larger values use more pixels for model fitting but also increase time and memory requirements.

  - `umap_metric`: method for calculating distances between pixel-level mass-spectrometry features, default `"cosine"`. Keep the default for your first use.
    - `“cosine”`: Focuses more on the relative directions of two spectral features rather than their overall amplitude. This is typically the preferred choice for MSI data that has already undergone intensity normalization and where the goal is to highlight spectral pattern differences;
  
    - `“euclidean”`: Calculates distance based on the absolute numerical differences between features. This can be used for comparison when the differences in feature amplitude themselves are intended to contribute to region differentiation;
  
  - `n_components`: the current spatial color image uses exactly three features, corresponding to the RGB channels.
  
  - `random_state`: fixing the random seed makes repeated results easier to compare, although the UMAP implementation may run with less parallelism when a fixed seed is used.

This page explains only the UMAP part. For tissue-region segmentation, please see [KMeans Tissue Region Segmentation](./segmentation.md).

```python
from massflow.data_manager import MSDataManagerImzML
from massflow.segmentation import SegmentationPipeline


dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()

pipeline = SegmentationPipeline(
    dm,
    output_dir="segmentation_output",
)

result = pipeline.run_umap_kmeans(
    sample_ratio=0.1,
    n_clusters=10,
    save_images=True,
    show=False,
)

dm.close()
```

The run produces both UMAP and KMeans results. When reviewing UMAP, focus on `umap_image.png`.

The following UMAP images were obtained with three different `sample_ratio` values:

`sample_ratio=0.1`:

![Snipaste_2026-08-090391211](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-090391211.jpg_view)

`sample_ratio=0.5`:

![image-20260815174917455](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/20260815174917515.jpg_view)

`sample_ratio=1.0`:

![image-20260815211057250](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/20260815211057364.jpg_view)

## How to Read a UMAP Image

- Pixels with similar colors in space have similar UMAP features.
- Locations where the color changes abruptly may correspond to spatial boundaries where the mass-spectrometry composition changes more clearly.
- Continuous colors represent continuous differences; they do not mean that tissue classes have already been segmented.
- The same color has no fixed biological meaning across different datasets or parameter settings.

UMAP results should be interpreted together with TIC images, known ion images, H&E images, or other experimental information.

## Output

When image saving is enabled by default, the result directory contains:

```text
segmentation_output/umap_kmeans/
├── umap_image.png
└── kmeans_image.png
```

`umap_image.png` is the result discussed on this page.

## References

- [`segmentation/segmentation_pipeline.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/segmentation/segmentation_pipeline.py): `SegmentationPipeline` and `UmapKmeansResult`
- [`segmentation/umap_kmeans.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/segmentation/umap_kmeans.py): UMAP/KMeans analysis and spatial-image output
