from dataclasses import dataclass
from venv import logger
import numpy as np
from massflow.module.ms_swaper_base import SpectrumSwaperBase
from massflow.logger import get_logger

logger = get_logger(__name__)

@dataclass
class SpectrumSwaperBinary(SpectrumSwaperBase):
    # Inherited swap_file_path, offset, length fields
    # Inherited mz_offset method

    def intensity_offset(self) -> int:
        """Calculate intensity start address based on mz dtype"""
        if self.mz_dtype is None:
            raise ValueError("spectrum is empty or mz_dtype is None")
        return self.offset + self.length * np.dtype(self.mz_dtype).itemsize

    # plan a
    def read_block(self, mz_dtype, intensity_dtype):
        """Read block at once and split into two arrays (zero copy)"""

        if self.mz_dtype is None or self.intensity_dtype is None:
            raise ValueError("spectrum is empty or dtype is None")
        if self.length <= 0:
            raise ValueError("length <= 0, cannot read block")

        mz_itemsize = np.dtype(mz_dtype).itemsize
        int_itemsize = np.dtype(intensity_dtype).itemsize

        mz_bytes = self.length * mz_itemsize
        int_bytes = self.length * int_itemsize
        total_bytes = mz_bytes + int_bytes

        # Map the entire block at once
        block = np.memmap(
            self.swap_file_path,
            dtype=np.uint8,
            mode="r",
            offset=self.offset,
            shape=(total_bytes,),
        )

        # Extract mz
        mz = np.frombuffer(block[:mz_bytes], dtype=mz_dtype, count=self.length)

        # Extract intensity
        intensity = np.frombuffer(
            block[mz_bytes : mz_bytes + int_bytes],
            dtype=intensity_dtype,
            count=self.length,
        )
        
        # update the spectrum
        self.spectrum.mz_list = mz
        self.spectrum.intensity = intensity


    # plan b
    # def read_block(self, mz_dtype, intensity_dtype):
    #     """Read block at once and split into two arrays (zero copy)"""

    #     intensity = np.memmap(
    #         self.swap_file_path,
    #         dtype=intensity_dtype,
    #         mode="r",
    #         offset=self.intensity_offset(mz_dtype),
    #         shape=(self.length,)
    #     )
    #     mz = np.memmap(
    #         self.swap_file_path,
    #         dtype=mz_dtype,
    #         mode="r",
    #         offset=self.offset,
    #         shape=(self.length,)
    #     )
    #     return mz, intensity

    def swap_out2disk(
        self, mz_array, intensity_array, mz_dtype, intensity_dtype, file_obj
    ):
        """
        Write block (file write position controlled by external file_obj, append or specific seek)
        Update offset and length after writing
        """
        if mz_array is not None and intensity_array is not None:
            assert len(mz_array) == len(intensity_array)
            self.length = len(mz_array)
            mz_dtype = mz_array.dtype
            intensity_dtype = intensity_array.dtype

        else:
            logger.error()
        # Determine the start point for block writing
        self.offset = file_obj.tell()

        # Write mz
        f_mz = np.asarray(mz_array, dtype=mz_dtype).tobytes()
        file_obj.write(f_mz)

        # Write intensity
        f_int = np.asarray(intensity_array, dtype=intensity_dtype).tobytes()
        file_obj.write(f_int)
        
        # update 
        # file_obj.flush()