from dataclasses import dataclass
import os
import numpy as np
from massflow.module.spectrum_swaper_base import SpectrumSwaperBase
from massflow.logger import get_logger

logger = get_logger(__name__)

@dataclass
class SpectrumSwaperBinary(SpectrumSwaperBase):
    # Inherited offset
    # Inherited properties/methods from base class

    def intensity_offset(self) -> int:
        """Calculate intensity start address based on mz dtype"""
        if self.mz_dtype is None:
            raise ValueError("spectrum is empty or mz_dtype is None")
        return self.offset + self.length * np.dtype(self.mz_dtype).itemsize

    # Plan: read the entire spectrum block at once and split into mz and intensity arrays
    def read_block(self):
        """Read the entire block and split into mz and intensity arrays (zero-copy).

        This method memory-maps the region and should produce two views or copies
        depending on implementation. It expects a contiguous block of mz data followed
        by intensity data.
        """

        if self.mz_dtype is None or self.intensity_dtype is None:
            logger.error("spectrum is empty or dtype is None")
            raise ValueError("spectrum is empty or dtype is None")
        if self.length <= 0:
            logger.error("length <= 0, cannot read block")
            raise ValueError("length <= 0, cannot read block")
        if self.swap_file_path is None:
            logger.error("swap_file_path is None, cannot read block")
            raise ValueError("swap_file_path is None, cannot read block")

        mz_itemsize = self.mz_dtype.itemsize
        int_itemsize = self.intensity_dtype.itemsize

        mz_bytes = self.length * mz_itemsize
        int_bytes = self.length * int_itemsize
        total_bytes = mz_bytes + int_bytes

        # Map the entire block at once
        block = np.memmap(
            self.swap_file_path,
            dtype=np.uint8,
            mode="r+",
            offset=self.offset,
            shape=(total_bytes,),
        )

    def swap_out2disk(self, fd: int):
        """Write spectrum block to disk using a raw file descriptor.

        The caller controls the lifetime and open mode of the fd. This method
        only performs low-level seeks and writes of raw bytes.
        """

        # pre check for data consistency
        if (
            self.spectrum is None
            or self.spectrum.mz_list is None
            or self.spectrum.intensity is None
        ):
            logger.error("Both spectrum mz_list and intensity cannot be None")
            raise ValueError("Both spectrum mz_list and intensity cannot be None")

        assert len(self.spectrum.mz_list) == self.length

        # Determine the start point for block writing (current file position of fd)
        self.offset = os.lseek(fd, 0, os.SEEK_CUR)

        logger.debug(
            f"Swapping out spectrum to disk at offset: {self.offset}, length: {self.length}"
        )

        # Prepare bytes
        f_mz = np.asarray(self.spectrum.mz_list, dtype=self.mz_dtype).tobytes()
        f_int = np.asarray(self.spectrum.intensity, dtype=self.intensity_dtype).tobytes()

        # Write mz and intensity using os.write
        os.write(fd, f_mz)
        os.write(fd, f_int)

        # free memory
        self.spectrum.mz_list = None  # free memory
        self.spectrum.intensity = None  # free memory