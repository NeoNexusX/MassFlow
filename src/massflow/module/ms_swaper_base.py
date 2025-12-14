from dataclasses import dataclass
from abc import ABC, abstractmethod
from typing import Optional
import numpy as np


@dataclass
class SpectrumSwaperBase(ABC):
    spectrum: object  # SpectrumBaseModule
    offset: int = 0
    length: int = -1

    @property
    def swap_file_path(self) -> str:
        return getattr(self.spectrum, "swap_file_path")

    @property
    # design target:
    # 1.dynamic length not bind with spectrum data, only get first time when length <=0
    def length(self) -> int:
        if self.length <= 0:
            mz = getattr(self.spectrum, "mz_list", None)
            return 0 if mz is None else len(mz)
        return self.length

    @property
    def mz_dtype(self) -> Optional[np.dtype]:
        mz = getattr(self.spectrum, "mz_list", None)
        return None if mz is None else np.dtype(mz.dtype)

    @property
    def intensity_dtype(self) -> Optional[np.dtype]:
        inten = getattr(self.spectrum, "intensity", None)
        return None if inten is None else np.dtype(inten.dtype)

    @abstractmethod
    def intensity_offset(self) -> int:
        pass

    @abstractmethod
    def swap_out2disk(self, file_obj):
        pass
