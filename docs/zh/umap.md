# UMAP 降维与可视化

本文介绍 MassFlow 中 UMAP 的作用和结果解释。UMAP 是用于高维质谱特征降维和可视化的方法。



作者：董乐强，董京伟

**问题反馈：** 如果在使用 MassFlow 或阅读文档时发现问题，可通过 [GitHub Issues](https://github.com/NeoNexusX/MassFlow/issues) 提交反馈，或联系邮箱：707244007@qq.com。

## 概述

一张 MSI 数据中，每个像素都包含大量 *m/z* 位置及对应强度。直接比较这些高维质谱特征比较困难。

UMAP 会把每个像素的高维质谱特征压缩到较少的特征维度，使质谱整体特征相近的像素在低维空间中也尽量保持接近。MassFlow 当前将 3 个 UMAP 特征映射到 RGB 颜色，并重新放回像素的空间坐标中，从而形成 UMAP 空间图。

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
    n_clusters: int = 10,
    random_state: int | None = None,
    show: bool = False,
) -> UmapKmeansResult
```

- 功能：根据每个像素的质谱特征运行 UMAP 和 KMeans，并生成空间聚类结果。
- 参数

  - `sample_ratio`：用于建立 UMAP 模型的像素比例，默认 `1.0`。数值越大，参与建模的像素越多，但时间和内存开销也会增加；
  - `umap_metric`：像素间质谱特征距离的计算方式，默认 `"cosine"`。第一次使用建议保持默认；
    - `"cosine"`：更关注两个像素质谱特征的相对组成方向，而不是整体幅度大小。对于已经做过强度归一化、希望突出谱型差异的 MSI 数据，通常可以优先使用；
    - `"euclidean"`：根据各特征之间的绝对数值差计算距离。当特征幅度差异本身也希望参与区域区分时可以比较使用；
  - `n_components`：当前空间彩色图固定使用 3 个特征，对应 RGB 三个通道；
  - `random_state`：固定随机种子可便于重复比较结果，但当前 UMAP 在固定种子时可能降低并行速度。

本页只解释 UMAP部分；组织区域划分见 [KMeans 组织区域划分](./segmentation.md)。

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
