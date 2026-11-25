# UMAP 降维与可视化

本文介绍 MassFlow 中 UMAP 的作用和结果解释。UMAP 是用于高维质谱特征降维和可视化的方法。



作者：董乐强，董京伟

**问题反馈：** 如果在使用 MassFlow 或阅读文档时发现问题，可通过 [GitHub Issues](https://github.com/NeoNexusX/MassFlow/issues) 提交反馈，或联系邮箱：707244007@qq.com。

## 概述

一张 MSI 数据中，每个像素都包含大量 *m/z* 位置及对应强度。直接比较这些高维质谱特征比较困难。

UMAP 会把每个像素的高维质谱特征压缩到较少的特征维度，使质谱整体特征相近的像素在低维空间中也尽量保持接近。MassFlow 将 UMAP 特征的前三维映射到 RGB 颜色，并重新放回像素的空间坐标中，从而形成 UMAP 空间图。

UMAP 图展示的是不同像素的整体质谱特征是否相似，以及这种差异在组织空间中怎样分布。

它本身不会给区域命名，也不能直接说明某种颜色就是肿瘤、坏死或正常组织。

## 当前 API

当前流程为通过 `SegmentationPipeline.run_umap_kmeans()` 一次性计算 UMAP 和 KMeans。虽然 API 同时返回两种结果，但二者作用不同：

- UMAP：降维和连续特征可视化；
- KMeans：基于低维特征进行离散类别划分。

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

- 功能：根据每个像素的质谱特征运行 UMAP 和 KMeans，并生成空间聚类结果。
- 参数

  - `sample_ratio`：用于建立 UMAP 模型的像素比例，默认 `1.0`。数值越大，参与建模的像素越多，但时间和内存开销也会增加；
  - `umap_metric`：像素间质谱特征距离的计算方式，默认 `"cosine"`。第一次使用建议保持默认；
    - `"cosine"`：更关注两个像素质谱特征的相对组成方向，而不是整体幅度大小。对于已经做过强度归一化、希望突出谱型差异的 MSI 数据，通常可以优先使用；
    - `"euclidean"`：根据各特征之间的绝对数值差计算距离。当特征幅度差异本身也希望参与区域区分时可以比较使用；
  - `n_components`：UMAP 降维后的特征维度，默认 `3`，必须为正整数。KMeans 聚类使用全部维度；UMAP 空间图取前三维对应 RGB 三个通道，维度为一时输出灰度图，维度为二时只使用 R、G 两个通道；
  - `n_neighbors` 控制构建邻域图时考虑的近邻数量，默认 `15`。较小的值更关注局部结构，较大的值更关注整体结构。
  - `min_dist` 和 `spread` 共同控制嵌入中点的分布。`min_dist` 默认 `0.1`，较小的值允许相近的点更紧密地排列。`spread` 默认 `1.0`，控制嵌入的距离尺度。`min_dist` 必须满足 `0 <= min_dist <= spread`。
  - `n_epochs` 控制优化轮数，默认 `None`，由 UMAP 根据拟合数据量选择。`learning_rate` 控制初始学习率，默认 `1.0`。
  - `umap_init` 默认 `"spectral"`，也支持 `"random"`、`"pca"`、`"tswspectral"` 或 NumPy 数组。数组的形状应为 `(参与拟合的像素数, n_components)`。
  - `repulsion_strength` 控制优化中负样本的权重，默认 `1.0`。`negative_sample_rate` 控制每个正样本对应的负样本数，默认 `5`，提高该值会增加优化开销。
  - `random_state` 默认 `None`。传入整数会同时为 UMAP 和 KMeans 设置随机种子，UMAP 拟合会使用单线程。

单独生成 UMAP 空间图时，可通过 `massflow.segmentation.plot_umap_image()` 设置相同的 UMAP 模型参数。KMeans 参数见 [KMeans 组织区域划分](./segmentation.md)。

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

运行后会同时得到 UMAP 和 KMeans 结果。查看 UMAP 时重点关注 `umap_image.png`。

以下为设置三种不同sample_ratio所得到的UMAP图：

sample_ratio=0.1：

![Snipaste_2026-08-090391211](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-090391211.jpg_view)

sample_ratio=0.5

![image-20260815174917455](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/20260815174917515.jpg_view)

sample_ratio=1.0

![image-20260815211057250](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/20260815211057364.jpg_view)

## UMAP 图怎样看

- 空间上颜色接近的像素，表示它们的 UMAP 特征较接近；
- 颜色突然变化的位置可能对应质谱组成发生较明显变化的空间边界；
- 连续颜色表示连续差异，并不是已经完成了组织类别划分；
- 当 `n_components` 大于 3 时，UMAP 图的颜色只反映前三维，而 KMeans 划分基于全部维度，因此颜色相近的像素也可能被分到不同类别；
- 不同数据集或不同参数下，相同颜色没有固定生物学含义。

UMAP 结果应结合 TIC、已知离子图、H&E 或其他实验信息解释。

## 输出结果

默认保存图像时，结果目录中会包含：

```text
segmentation_output/umap_kmeans/
├── umap_image.png
└── kmeans_image.png
```

其中 `umap_image.png` 是本页关注的结果。 

## 参考

- [`segmentation/segmentation_pipeline.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/segmentation/segmentation_pipeline.py)：`SegmentationPipeline` 和 `UmapKmeansResult`
- [`segmentation/umap_kmeans.py`](https://github.com/NeoNexusX/MassFlow/blob/main/src/massflow/segmentation/umap_kmeans.py)：UMAP/KMeans 分析和空间图输出
