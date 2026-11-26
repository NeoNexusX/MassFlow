# UMAP Dimensionality Reduction and Visualization

This page introduces the role of UMAP in MassFlow and how to interpret its results. UMAP is used for dimensionality reduction and visualization of high-dimensional mass-spectrometry features.



Authors: Leqiang Dong, Jingwei Dong

**Feedback:** If you find any problems while using MassFlow or reading the documentation, you can submit feedback through [GitHub Issues](https://github.com/NeoNexusX/MassFlow/issues), or contact us by email at 707244007@qq.com.

## Overview

In an MSI dataset, each pixel contains a large number of *m/z* positions and corresponding intensities. Directly comparing these high-dimensional mass-spectrometry features is difficult.

UMAP compresses the high-dimensional spectral features of each pixel into a smaller number of feature dimensions, while trying to keep pixels with similar overall spectral characteristics close to one another in the low-dimensional space. MassFlow maps the first three UMAP features to the RGB color channels and places them back at the pixels' spatial coordinates to form a spatial UMAP image.

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
    sample_ratio: float = 1.0,
    transform_batch_size: int | None = 1024,
    dpi: int = 600,
    save_images: bool = True,
    save_matrix: Literal["npz", "zarr"] | None = None,
    n_components: int = 3,
    umap_metric: str = "cosine",
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    spread: float = 1.0,
    n_epochs: int | None = None,
    learning_rate: float = 1.0,
    umap_init: str | np.ndarray = "spectral",
    repulsion_strength: float = 1.0,
    negative_sample_rate: int = 5,
    n_clusters: int = 10,
    kmeans_init: KMeansInit = "k-means++",
    n_init: int | Literal["auto"] = "auto",
    max_iter: int = 300,
    tol: float = 1e-4,
    algorithm: Literal["lloyd", "elkan"] = "lloyd",
    random_state: int | None = None,
    show: bool = True,
) -> UmapKmeansResult
```

- Function: runs UMAP and KMeans using the mass-spectrometry features of each pixel and generates spatial clustering results.

- Parameters:

  - `sample_ratio`: proportion of pixels used to build the UMAP model, default `1.0`. Larger values use more pixels for model fitting but also increase time and memory requirements.

  - `umap_metric`: method for calculating distances between pixel-level mass-spectrometry features, default `"cosine"`. Keep the default for your first use.
    - `“cosine”`: Focuses more on the relative directions of two spectral features rather than their overall amplitude. This is typically the preferred choice for MSI data that has already undergone intensity normalization and where the goal is to highlight spectral pattern differences;
  
    - `“euclidean”`: Calculates distance based on the absolute numerical differences between features. This can be used for comparison when the differences in feature amplitude themselves are intended to contribute to region differentiation;
  
  - `n_components`: number of feature dimensions produced by UMAP, default `3`. It must be a positive integer. KMeans uses all dimensions, while the spatial UMAP image maps the first three to the RGB channels. With one dimension the image is grayscale; with two dimensions only the red and green channels are used.
  
  - `n_neighbors` controls how many neighbors UMAP considers when building the neighborhood graph. The default is `15`. Smaller values emphasize local structure, while larger values capture more global structure.
  - `min_dist` and `spread` shape the distances in the embedding. `min_dist` defaults to `0.1`, and smaller values allow nearby points to pack more tightly. `spread` defaults to `1.0` and sets the distance scale. `min_dist` must be between `0` and `spread`, inclusive.
  - `n_epochs` sets the number of optimization epochs. The default `None` lets UMAP select a value based on the number of fitted pixels. `learning_rate` sets the initial learning rate and defaults to `1.0`.
  - `umap_init` defaults to `"spectral"`. It also accepts `"random"`, `"pca"`, `"tswspectral"`, or a NumPy array of shape `(number of fitted pixels, n_components)`.
  - `repulsion_strength` weights negative samples during optimization and defaults to `1.0`. `negative_sample_rate` sets the number of negatives per positive sample and defaults to `5`. Higher rates increase optimization cost.
  - `random_state` defaults to `None`. Passing an integer seeds both UMAP and KMeans and makes UMAP fitting single-threaded.

To generate a spatial UMAP image on its own, pass the same UMAP model parameters to `massflow.segmentation.plot_umap_image()`. KMeans parameters are documented separately. [KMeans Tissue Region Segmentation](./segmentation.md)

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
    sample_ratio=1.0,
    n_neighbors=15,
    min_dist=0.1,
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
- When `n_components` is greater than 3, the UMAP image colors reflect only the first three dimensions while KMeans divides pixels using all dimensions, so pixels with similar colors may still fall into different clusters.
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
