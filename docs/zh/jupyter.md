# 最佳实践

本文介绍在使用 MassFlow 进行 KMeans 组织区域划分前，需要掌握的 Jupyter Notebook 基本操作。

[TOC]

## 概述

普通 Python 脚本通常会从第一行一直运行到最后一行。

Jupyter Notebook 则可以把代码拆分成多个代码单元格，每个单元格都可以单独运行。前面单元格中已经创建的数据和变量，会保留在当前 Notebook 中，后面的单元格可以直接继续使用。

例如，第一个单元格可以先读取 MSI 数据：

```python
from massflow.data_manager import MSDataManagerImzML

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()
```

运行完成后，`dm` 会保留在当前 Notebook 中。

后面的单元格可以直接继续使用：

```
print(dm)
```

不需要重新读取数据。



## 1.准备 Jupyter Notebook

Jupyter Notebook 应与 MassFlow 使用同一个 Python 环境。

如果当前 MassFlow 环境中还未安装 Jupyter，可以在 MassFlow 项目目录的终端中安装：

```
uv pip install jupyter
```

安装完成后运行：

```
uv run jupyter notebook
```

如果是使用 Jupyterlab 也可以运行：

```
uv run jupyter lab
```



启动成功后，会在浏览器中打开 Jupyter 页面。

如图所示：

![Snipaste_2026-08-12_13-21-49](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_13-21-49.jpg_view)

## 2.新建Notebook

打开 Jupyter 后，新建一个 Python Notebook。

建议将 Notebook 放在 MassFlow 项目目录中，例如：

```
MassFlow/
├── data/
│   └── example.imzML
├── src/
├── tests/
└── segmentation.ipynb
```

步骤如图所示：

![Snipaste_2026-08-12_20-08-05](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_20-08-05.jpg_view)

## 3. 运行代码单元格

Jupyter Notebook 中每一个代码框都是一个代码单元格，前面单元格中已经创建的变量，可以在后面的单元格中继续使用。

例如，第一个单元格中创建：

```
from massflow.data_manager import MSDataManagerImzML
from massflow.segmentation import SegmentationPipeline

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()

pipeline = SegmentationPipeline(dm)
```

运行后，当前 Notebook 中会创建两个变量：`dm` 和 `pipeline`。

所以下一个单元格可以直接运行：

```
result = pipeline.run_umap_kmeans(
    n_clusters=10,
)
```

不需要再次读取数据，也不需要重新创建 `pipeline`。

后面如果已经确认需要处理的类别，也可以继续直接运行：

```
pipeline.remove_clusters([0])
```

## 4. 需要按顺序运行

后面的代码通常依赖前面的运行结果。因此第一次使用时应从上到下依次运行。

例如：

```
result = pipeline.run_umap_kmeans(
    sample_ratio=0.1,
    n_clusters=10,
    random_state=0,
    save_images=True,
    save_matrix="zarr",
    show=False,
)
```

需要前面已经创建好 `pipeline`。

如果没有先运行读取数据的单元格，就直接运行这段代码，会出现报错：

![Snipaste_2026-08-12_20-20-34](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_20-20-34.jpg_view)

这表示当前 Notebook 中还没有创建 `pipeline`。

因此，第一次运行 Notebook 时，建议按照代码单元格依次运行。

## 5. 修改参数后要重新运行

修改一个代码单元格中的参数后，需要重新运行该单元格，新的设置才会生效。

例如原来是：

```
result = pipeline.run_umap_kmeans(
    n_clusters=5,
)
```

修改为：

```
result = pipeline.run_umap_kmeans(
    n_clusters=8,
)
```

修改代码本身不会自动重新计算。

需要再次运行这个单元格，然后重新查看新的 KMeans 结果。

如果后面的区域选择依赖新的聚类结果，也应重新检查并运行相关的后续步骤。

## 6. Notebook 重启后变量会被清除

Notebook 中创建的`dm`,`pipeline`,`result`等变量只保存在当前运行环境中。

如果关闭 Jupyter 或重启 Kernel，这些变量会被清除。

即使页面中的代码仍然存在，也需要重新运行前面的代码单元格。

例如出现：

![Snipaste_2026-08-12_20-20-34](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_20-20-34.jpg_view)

可能是因为：

- 创建 `pipeline` 的单元格还没有运行；
- Jupyter 的 Kernel 已经重新启动。

此时从前面的数据读取步骤开始重新运行即可。

## 下一步

掌握以上 Jupyter Notebook 基本操作后，可以继续阅读：

[KMeans 组织区域划分](./segmentation.md)
