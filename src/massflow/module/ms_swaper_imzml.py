from dataclasses import dataclass
from typing import Any
from .ms_swaper_base import SpectrumSwaperBase
from massflow.logger import get_logger
from pyimzml.ImzMLWriter import ImzMLWriter
import numpy as np

logger = get_logger(__name__)

@dataclass
class SpectrumSwaperImzml(SpectrumSwaperBase):

    def swap_out2disk(self, writer: ImzMLWriter):

        if (
            self.spectrum is None
            or self.spectrum.mz_list is None
            or self.spectrum.intensity is None
        ):
            logger.error("Both spectrum mz_list and intensity cannot be None")
            raise ValueError("Both spectrum mz_list and intensity cannot be None")

        assert len(self.spectrum.mz_list) == self.length

        f_mz = self.spectrum.mz_list
        f_intensitys = self.spectrum.intensity

        coords_tuple = (
            self.spectrum.coordinates.x,
            self.spectrum.coordinates.y,
            self.spectrum.coordinates.z
        )

        writer.addSpectrum(f_mz, f_intensitys, coords_tuple)

        # 获取刚写入的谱图信息
        current_spec = writer.spectra[-1]
        self.mz_offset = current_spec.mz_offset
        self.int_offset = current_spec.int_offset
        self.length = current_spec.mz_len

        self.spectrum.mz_list = None
        self.spectrum.intensity = None

    def read_block(self):

        ibd_path = self.swap_file_path.with_suffix(".ibd")

        mz_itemsize = self.mz_dtype.itemsize
        int_itemsize = self.intensity_dtype.itemsize

        mz_bytes = self.length * mz_itemsize
        int_bytes = self.length * int_itemsize

        with open(ibd_path, "rb") as f:
            
            # 读取mz数据
            f.seek(self.mz_offset)
            mz_buffer = f.read(mz_bytes)
            mz = np.frombuffer(mz_buffer, dtype=self.mz_dtype)

            # 读取intensity数据
            f.seek(self.int_offset)
            int_buffer = f.read(int_bytes)
            intensity = np.frombuffer(int_buffer, dtype=self.intensity_dtype)

        # 更新数据
        self.spectrum.mz_list = mz.copy() # copy 确保拥有内存所有权
        self.spectrum.intensity = intensity.copy()

        # mz_block = np.memmap(
        #     ibd_path,
        #     dtype=np.uint8,
        #     mode="r+",
        #     offset=self.mz_offset,
        #     shape=(mz_bytes,),
        # )

        # int_block = np.memmap(
        #     ibd_path,
        #     dtype=np.uint8,
        #     mode="r+",
        #     offset=self.int_offset,
        #     shape=(int_bytes,),
        # )
        
        # # Extract mz
        # mz = np.frombuffer(mz_block[:], dtype=self.mz_dtype, count=self.length)

        # # Extract intensity
        # intensity = np.frombuffer(
        #     int_block[:],
        #     dtype=self.intensity_dtype,
        #     count=self.length,
        # )
        
        # # update the spectrum
        # self.spectrum.mz_list = mz 
        # self.spectrum.intensity = intensity