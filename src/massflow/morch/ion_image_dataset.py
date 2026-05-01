"""Abstract ion image dataset interface for PyTorch training."""

from __future__ import annotations

from abc import ABC, abstractmethod


class IonImageDataset(ABC):
    """Abstract interface for fast random-access 2D ion image loading.

    Compatible with ``torch.utils.data.DataLoader`` via duck-typing
    (``__len__`` + ``__getitem__``).  Backends implement the actual
    storage and scatter logic.

    Pixels not covered by any coordinate default to 0.
    """

    @abstractmethod
    def __len__(self) -> int:
        """Total number of ion images."""

    @abstractmethod
    def __getitem__(self, idx: int):
        """Return ion image ``idx`` as a ``(1, H, W)`` torch tensor."""

    @property
    @abstractmethod
    def height(self) -> int:
        """Image height in pixels."""

    @property
    @abstractmethod
    def width(self) -> int:
        """Image width in pixels."""

    @property
    @abstractmethod
    def n_ions(self) -> int:
        """Number of ion images."""
