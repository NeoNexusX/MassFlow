"""Binary-backend ion image dataset."""

from __future__ import annotations

import numpy as np

from massflow.data_manager import IonImageDataManagerBinary
from massflow.morch.ion_image_dataset import IonImageDataset


class IonImageBinaryDataset(IonImageDataset):
    """Fast random-access 2D ion image loader backed by the binary format.

    Computes scatter indices once from coordinates, then scatters ion
    vectors into 2D images on each ``__getitem__`` with a single
    vectorised assignment.

    Multi-worker safe: each worker lazily opens its own memmap-backed
    data manager via ``_ensure_dm()``, avoiding pickle issues with file
    handles.

    Usage::

        from massflow.data_manager import IonImageDataManagerBinary
        from massflow.morch import IonImageDatasetBinary
        from torch.utils.data import DataLoader

        dm = IonImageDataManagerBinary("output.iib", mode="r")
        ds = IonImageDatasetBinary(dm)

        loader = DataLoader(ds, batch_size=32, shuffle=True,
                            num_workers=4, pin_memory=True)
        for batch in loader:
            pred = model(batch.cuda())  # batch: (B, 1, H, W)
    """

    def __init__(self, data_manager: IonImageDataManagerBinary) -> None:
        if data_manager.height is None or data_manager.width is None:
            raise ValueError("data_manager height/width must be set")
        if data_manager.coordinates is None:
            raise ValueError("data_manager coordinates must be loaded")

        self._filepath = str(data_manager.filepath)
        self._height = int(data_manager.height)
        self._width = int(data_manager.width)
        self._np_dtype = data_manager.intensity_dtype
        self._n_ions = data_manager.n_ions

        coords = np.array(data_manager.coordinates)
        xs = coords[:, 0].astype(np.int64, copy=False)
        ys = coords[:, 1].astype(np.int64, copy=False)

        if xs.size > 0 and int(xs.min()) != 0:
            xs = xs - 1
            ys = ys - 1

        valid = (xs >= 0) & (xs < self._width) & (ys >= 0) & (ys < self._height)
        self._flat_idx = (ys[valid] * self._width + xs[valid])
        self._valid_px = np.where(valid)[0]

        self._dm: IonImageDataManagerBinary | None = None
        self._buffer: np.ndarray | None = None

    # -- IonImageDataset interface -------------------------------------------

    @property
    def height(self) -> int:
        return self._height

    @property
    def width(self) -> int:
        return self._width

    @property
    def n_ions(self) -> int:
        return self._n_ions

    def __len__(self) -> int:
        return self._n_ions

    def __getitem__(self, idx: int): #type: ignore[override]
        import torch #type: ignore[import]

        if idx < 0:
            idx += self._n_ions

        dm = self._ensure_dm()

        if self._buffer is None:
            self._buffer = np.zeros(self._height * self._width, dtype=self._np_dtype)

        self._buffer.fill(0)
        self._buffer[self._flat_idx] = dm.intensity_matrix[idx, self._valid_px]

        return torch.from_numpy(self._buffer).reshape(1, self._height, self._width).clone()

    # -- internals -----------------------------------------------------------

    def _ensure_dm(self) -> IonImageDataManagerBinary:
        if self._dm is None:
            self._dm = IonImageDataManagerBinary(self._filepath, mode="r")
        return self._dm
