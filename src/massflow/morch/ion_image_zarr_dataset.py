"""PyTorch dataset for ion-major MassFlow MSI Zarr stores."""

from __future__ import annotations

import numpy as np

from massflow.morch.ion_image_binary_dataset import IonImageIndexDataset
from massflow.morch.ion_image_dataset import IonImageDataset
from massflow.tools.zarr_reader import MSIZarrReader


class IonImageZarrDataset(IonImageDataset):
    """Random-access ion images reconstructed from compact ion vectors."""

    def __init__(self, data_manager) -> None:
        self._path = str(data_manager.filepath)
        self._height = int(data_manager.height)
        self._width = int(data_manager.width)
        self._n = int(data_manager.n_ions)
        self._np_dtype = np.dtype(data_manager.intensity_dtype)
        coordinates = np.asarray(data_manager.coordinates, dtype=np.int64)
        self._flat_idx = (
            (coordinates[:, 1] - 1) * self._width + coordinates[:, 0] - 1
        )
        self._reader: MSIZarrReader | None = None

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
        return IonImageIndexDataset(self._n)

    def __getitem__(self, idx: int):  # type: ignore[override]
        import torch  # type: ignore[import]

        idx = self._normalize_index(idx)
        values = self._ensure_reader().read_ion_vector(idx)
        image = np.zeros((self._height, self._width), dtype=self._np_dtype)
        image.ravel()[self._flat_idx] = values
        return torch.from_numpy(image).unsqueeze(0)

    def get_batch(self, indices):
        import torch  # type: ignore[import]

        idx = np.asarray(indices, dtype=np.int64)
        if idx.ndim == 0:
            idx = idx.reshape(1)
        if idx.ndim != 1 or idx.size == 0:
            raise ValueError("indices must be a non-empty 1-D sequence.")
        idx = np.asarray([self._normalize_index(int(value)) for value in idx])
        reader = self._ensure_reader()
        rows = np.vstack([reader.read_ion_vector(int(value)) for value in idx])
        batch = np.zeros(
            (idx.size, 1, self._height, self._width),
            dtype=self._np_dtype,
        )
        batch[:, 0].reshape(idx.size, -1)[:, self._flat_idx] = rows
        return torch.from_numpy(batch)

    def _ensure_reader(self) -> MSIZarrReader:
        if self._reader is None:
            self._reader = MSIZarrReader(self._path)
        return self._reader

    def _normalize_index(self, idx: int) -> int:
        if idx < 0:
            idx += self._n
        if idx < 0 or idx >= self._n:
            raise IndexError("ion index out of range")
        return idx
