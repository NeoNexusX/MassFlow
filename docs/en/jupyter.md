# Best Practices

This page introduces the basic Jupyter Notebook operations that are useful before performing KMeans tissue-region segmentation with MassFlow.



If you  want to try MSI data processing online or search for MSI datasets of interest, visit https://sxo.bionet.pro and register for free.

<img src="https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/20260815222635996.img_view" alt="image-20260815222635943" style="zoom: 25%;" />

## Overview

A normal Python script usually runs continuously from the first line to the last.

Jupyter Notebook can divide code into multiple code cells, and each cell can be run separately. Data and variables created in earlier cells remain available in the Notebook and can be used directly by later cells.

For example, the first cell can load MSI data:

```python
from massflow.data_manager import MSDataManagerImzML

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()
```

After this cell runs, `dm` remains available in the Notebook.

A later cell can use it directly:

```
print(dm)
```

There is no need to read the data again.



## 1. Prepare Jupyter Notebook

Jupyter Notebook should use the same Python environment as MassFlow.

If Jupyter is not installed in the MassFlow environment, install it from a terminal opened in the MassFlow project directory:

```
uv pip install jupyter
```

After installation, run:

```
uv run jupyter notebook
```

If you prefer JupyterLab, you can run:

```
uv run jupyter lab
```



After Jupyter starts successfully, its page opens in your browser.

Example:

![Snipaste_2026-08-12_13-21-49](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_13-21-49.jpg_view)

## 2. Create a New Notebook

After opening Jupyter, create a new Python Notebook.

It is recommended to place the Notebook inside the MassFlow project directory, for example:

```
MassFlow/
├── data/
│   └── example.imzML
├── src/
├── tests/
└── segmentation.ipynb
```

The steps are shown below:

![Snipaste_2026-08-12_20-08-05](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_20-08-05.jpg_view)

## 3. Run Code Cells

Each code box in a Jupyter Notebook is a code cell. Variables created in earlier cells can continue to be used in later cells.

For example, create the following objects in the first cell:

```
from massflow.data_manager import MSDataManagerImzML
from massflow.segmentation import SegmentationPipeline

dm = MSDataManagerImzML(filepath="data/example.imzML")
dm.load_head_data()

pipeline = SegmentationPipeline(dm)
```

After this cell runs, two variables exist in the Notebook: `dm` and `pipeline`.

The next cell can run:

```
result = pipeline.run_umap_kmeans(
    n_clusters=10,
)
```

There is no need to read the data again or recreate `pipeline`.

Later, after you have confirmed which class needs to be processed, you can also run:

```
pipeline.remove_clusters([0])
```

## 4. Run Cells in Order

Later code usually depends on results created by earlier cells. When using a Notebook for the first time, run the cells from top to bottom in order.

For example:

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

This code requires `pipeline` to have already been created in an earlier cell.

If you run this code before running the data-loading cell, an error appears:

![Snipaste_2026-08-12_20-20-34](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_20-20-34.jpg_view)

This means that `pipeline` has not yet been created in the Notebook.

Therefore, when running a Notebook for the first time, it is recommended to execute the code cells in sequence.

## 5. Rerun a Cell After Changing Parameters

After changing a parameter in a code cell, you must rerun that cell for the new setting to take effect.

For example, suppose the original code is:

```
result = pipeline.run_umap_kmeans(
    n_clusters=5,
)
```

and you change it to:

```
result = pipeline.run_umap_kmeans(
    n_clusters=8,
)
```

Editing the code alone does not automatically recompute the result.

Run the cell again, then inspect the new KMeans result.

If later region-selection steps depend on the new clustering result, review and rerun the relevant downstream cells as well.

## 6. Variables Are Cleared When the Notebook Restarts

Variables such as `dm`, `pipeline`, and `result` exist only in the current running environment.

If Jupyter is closed or the Kernel is restarted, these variables are cleared.

Even though the code remains visible on the page, you still need to rerun the earlier code cells.

For example, if you see:

![Snipaste_2026-08-12_20-20-34](https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/Snipaste_2026-08-12_20-20-34.jpg_view)

possible reasons include:

- the cell that creates `pipeline` has not been run yet;
- the Kernel has been restarted.

In either case, rerun the Notebook starting from the data-loading steps.

## Next Step

After becoming familiar with these basic Jupyter Notebook operations, continue with:

[KMeans Tissue Region Segmentation](./segmentation.md)





Authors: Leqiang Dong, Jingwei Dong
