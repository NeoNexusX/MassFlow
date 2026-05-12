"""Zarr-backend ion image dataset with multi-worker safety."""

from __future__ import annotations

import numpy as np
import zarr

from massflow.morch.ion_image_binary_dataset import IonImageIndexDataset
from massflow.morch.ion_image_dataset import IonImageDataset


class IonImageZarrDataset(IonImageDataset):
    """Random-access dataset backed by Zarr ion image stores.

    Each item reads one dense ``(H, W)`` ion image from the ``ion_images`` array.
    Multi-worker safe: each worker lazily opens its own Zarr group.
    """

    def __init__(self, data_manager) -> None:
        self._path = str(data_manager.filepath)
        self._height = int(data_manager.height)
        self._width = int(data_manager.width)
        self._n = data_manager.n_ions
        self._np_dtype = np.dtype(getattr(data_manager, "intensity_dtype", np.float32))
        self._group = None
        self._dset = None

    @property
    def height(self) -> int:
        return self._height

    @property
    def width(self) -> int:
        return self._width

    @property
    def n_ions(self) -> int:
        return self._n

    def __len__(self) -> int:
        return self._n

    def index_dataset(self) -> IonImageIndexDataset:
        """Return an index-only dataset for batch-aware DataLoader usage."""
        return IonImageIndexDataset(self._n)

    def __getitem__(self, idx: int):  # type: ignore[override]
        import torch  # type: ignore[import]

        if idx < 0:
            idx += self._n
        if idx < 0 or idx >= self._n:
            raise IndexError("ion index out of range")

        image = np.asarray(self._ensure_dset()[idx])
        return torch.from_numpy(image).unsqueeze(0).clone()

    def get_batch(self, indices):
        """Read a batch of Zarr ion images and return ``(B, 1, H, W)``."""
        import torch  # type: ignore[import]

        idx = np.asarray(indices, dtype=np.int64)
        if idx.ndim == 0:
            idx = idx.reshape(1)
        if idx.ndim != 1 or idx.size == 0:
            raise ValueError("indices must be a 1-D sequence of ion indices")

        idx = idx.copy()
        idx[idx < 0] += self._n
        if int(idx.min()) < 0 or int(idx.max()) >= self._n:
            raise IndexError("ion index out of range")

        dset = self._ensure_dset()
        batch = np.empty((idx.size, 1, self._height, self._width), dtype=dset.dtype)
        for out_i, ion_idx in enumerate(idx):
            batch[out_i, 0] = np.asarray(dset[int(ion_idx)])
        return torch.from_numpy(batch)

    def _ensure_dset(self) -> zarr.Array:
        if self._dset is None:
            self._group = zarr.open_group(self._path, mode="r")
            self._dset = self._group["ion_images"]
        return self._dset # type: ignore[return-value]
