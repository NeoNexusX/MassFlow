from dataclasses import dataclass
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional,TYPE_CHECKING
import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

if TYPE_CHECKING:
    from massflow.module.ms_module import SpectrumBaseModule

@dataclass
class SpectrumSwaperBase(ABC):
    spectrum: "SpectrumBaseModule"
    # continuous 数据下考虑，还是每个spectrum单独存储mz/int的offset更合理
    mz_offset: int = -1
    int_offset: int = -1
    _length: int = -1
    _mz_dtype: Optional[np.dtype] = None
    _intensity_dtype: Optional[np.dtype] = None

    @property
    def swap_file_path(self) -> Path:
        return getattr(self.spectrum, "swap_file_path")

    @property
    # design target:
    # 1.dynamic length not bind with spectrum data, only get first time when length <=0
    def length(self) -> int:
        if self._length <= 0:
            mz = getattr(self.spectrum, "mz_list", None)
            self._length = 0 if mz is None else len(mz)
        return  self._length
    
    @length.setter
    def length(self, value):
        self._length = value
    
    @property
    def mz_dtype(self) -> Optional[np.dtype]:
        if self._mz_dtype is None:
            mz = getattr(self.spectrum, "mz_list", None)
            self._mz_dtype = None if mz is None else mz.dtype
        return self._mz_dtype
    
    @property
    def intensity_dtype(self) -> Optional[np.dtype]:
        if self._intensity_dtype is None:
            inten = getattr(self.spectrum, "intensity", None)
            self._intensity_dtype = None if inten is None else np.dtype(inten.dtype)
        return self._intensity_dtype

    @abstractmethod
    def intensity_offset(self) -> int:
        pass

    @abstractmethod
    def swap_out2disk(self, writer: ImzMLWriter):
        """Swap spectrum data out to disk using a raw file descriptor.

        Implementations should use low-level OS APIs (e.g. os.lseek/os.write)
        and operate on bytes only.
        """
        pass
