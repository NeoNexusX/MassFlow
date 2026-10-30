"""HDF5-backend ion image dataset with multi-worker safety."""

from __future__ import annotations

import h5py
import numpy as np

from massflow.morch.ion_image_binary_dataset import IonImageIndexDataset
from massflow.morch.ion_image_dataset import IonImageDataset


class IonImageH5Dataset(IonImageDataset):
    """Random-access dataset backed by bitshuffle+lz4 compressed H5 files.

    Each ``__getitem__`` reads a single H5 chunk (one ion image), decompresses
    it, and returns a ``(1, H, W)`` tensor. Multi-worker safe: each worker
    lazily opens its own ``h5py.File`` on first access.
    """

    def __init__(self, data_manager) -> None:
        # Properties auto-resolve via _open_for_read() and raise ValueError if
        # still uninitialised, so no manual None checks are needed here.
        self._path = str(data_manager.filepath)
        self._height = int(data_manager.height)
        self._width = int(data_manager.width)
        self._n = data_manager.n_ions
        self._np_dtype = np.dtype(getattr(data_manager, "intensity_dtype", np.float32))
        self._file: h5py.File | None = None
        self._dset: h5py.Dataset | None = None

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

    def __getitem__(self, idx: int):#type: ignore[override]
        import torch  # type: ignore[import]

        if idx < 0:
            idx += self._n
        dset = self._ensure_dset()
        img = dset[idx]  # (H, W) single chunk, decompressed on read
        return torch.from_numpy(np.asarray(img)).unsqueeze(0).clone()

    def get_batch(self, indices):
        """Read a batch of H5 ion images and return ``(B, 1, H, W)``.

        Intended for use as a PyTorch ``DataLoader`` ``collate_fn``. The
        DataLoader still controls sampling order; this method only converts the
        sampled indices into a ready-to-train batch tensor.
        """
        import torch  # type: ignore[import]

        idx = np.asarray(indices, dtype=np.int64)
        if idx.ndim == 0:
            idx = idx.reshape(1)
        if idx.ndim != 1 or idx.size == 0:
            raise ValueError("indices must be a 1-D sequence of ion indices")

        dset = self._ensure_dset()

        batch = np.empty((idx.size, 1, self._height, self._width), dtype=dset.dtype)
        for out_i, ion_idx in enumerate(idx):
            dset.read_direct(
                batch,
                source_sel=np.s_[int(ion_idx), :, :],
                dest_sel=np.s_[out_i, 0, :, :],
            )
        return torch.from_numpy(batch)

    def _ensure_dset(self) -> h5py.Dataset:
        if self._file is None:
            self._file = h5py.File(self._path, "r")
            self._dset = self._file["ion_images"]
        assert self._dset is not None
        return self._dset

    def __getstate__(self):
        """Drop live h5py handles so the dataset survives pickling (spawn-based
        DataLoader workers). Each worker reopens its own file lazily via
        ``_ensure_dset`` on first access."""
        state = self.__dict__.copy()
        state["_file"] = None
        state["_dset"] = None
        return state
