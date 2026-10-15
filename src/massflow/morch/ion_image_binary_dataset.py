"""Binary-backend ion image dataset."""

from __future__ import annotations

import numpy as np

from massflow.data_manager import IonImageDataManagerBinary
from massflow.morch.ion_image_dataset import IonImageDataset
from massflow.tools.funs import zero_based_xy


class IonImageIndexDataset:
    """Small DataLoader dataset that yields ion indices only."""

    def __init__(self, size: int) -> None:
        self._size = int(size)

    def __len__(self) -> int:
        return self._size

    def __getitem__(self, idx: int) -> int:
        return int(idx)


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
        from massflow.morch import IonImageBinaryDataset
        from torch.utils.data import DataLoader

        dm = IonImageDataManagerBinary("output.iib", mode="r")
        ds = IonImageBinaryDataset(dm)

        loader = DataLoader(ds, batch_size=32, shuffle=True,
                            num_workers=4, pin_memory=True)
        for batch in loader:
            pred = model(batch.cuda())  # batch: (B, 1, H, W)
    """

    def __init__(self, data_manager: IonImageDataManagerBinary) -> None:
        # Properties auto-resolve via _open_for_read() and raise ValueError if
        # still uninitialised, so no manual None checks are needed here.
        self._filepath = str(data_manager.filepath)
        self._height = int(data_manager.height)
        self._width = int(data_manager.width)
        self._np_dtype = data_manager.intensity_dtype
        self._n_ions = data_manager.n_ions

        xs, ys = zero_based_xy(data_manager.coordinates)
        valid = (xs >= 0) & (xs < self._width) & (ys >= 0) & (ys < self._height)
        self._flat_idx = (ys[valid] * self._width + xs[valid])
        self._valid_px = np.where(valid)[0]
        self._all_pixels_valid = bool(np.all(valid))

        self._dm: IonImageDataManagerBinary | None = None

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

    def index_dataset(self) -> IonImageIndexDataset:
        """Return an index-only dataset for batch-aware DataLoader usage."""
        return IonImageIndexDataset(self._n_ions)

    def __getitem__(self, idx: int): #type: ignore[override]
        import torch #type: ignore[import]

        if idx < 0:
            idx += self._n_ions

        dm = self._ensure_dm()

        image = np.zeros((self._height, self._width), dtype=self._np_dtype)
        image.ravel()[self._flat_idx] = dm.intensity_matrix[idx, self._valid_px]
        return torch.from_numpy(image).unsqueeze(0)

    def get_batch(self, indices):
        """Read a batch of ion images and return ``(B, 1, H, W)``.

        This is intended for use as a PyTorch ``DataLoader`` ``collate_fn``.
        The loader can yield integer indices, and this method performs the
        actual batched memmap read plus scatter in one pass.
        """
        import torch #type: ignore[import]

        idx = np.asarray(indices, dtype=np.int64)
        if idx.ndim == 0:
            idx = idx.reshape(1)
        if idx.ndim != 1:
            raise ValueError("indices must be a 1-D sequence of ion indices")

        if idx.size == 0:
            empty = np.empty((0, 1, self._height, self._width), dtype=self._np_dtype)
            return torch.from_numpy(empty)

        idx = idx.copy()
        idx[idx < 0] += self._n_ions
        if int(idx.min()) < 0 or int(idx.max()) >= self._n_ions:
            raise IndexError("ion index out of range")

        dm = self._ensure_dm()

        rows = np.asarray(dm.intensity_matrix[idx])
        if not self._all_pixels_valid:
            rows = rows[:, self._valid_px]

        batch = np.zeros((idx.size, 1, self._height, self._width), dtype=self._np_dtype)
        batch[:, 0].reshape(idx.size, -1)[:, self._flat_idx] = rows
        return torch.from_numpy(batch)

    # -- internals -----------------------------------------------------------

    def _ensure_dm(self) -> IonImageDataManagerBinary:
        if self._dm is None:
            self._dm = IonImageDataManagerBinary(self._filepath, mode="r")
        return self._dm
