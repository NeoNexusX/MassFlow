"""
Mass Spectrometry Module for MassFlow Framework

This module provides core classes and functionality for handling mass spectrometry (MS) data,
particularly for Mass Spectrometry Imaging (MSI) applications. It includes support for
lazy loading, efficient data management, and visualization capabilities.

Classes:
    Spectrum: Base class for mass spectrum data (imported from `massflow.module.spectrum`)
    SpectrumImzML: Specialized class for handling ImzML format with lazy loading

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""

from pathlib import Path
from massflow.module.spectrum import Spectrum
from typing import Optional
from massflow.logger import get_logger
logger = get_logger("ms_module")

class SpectrumImzML(Spectrum):
    """
    Specialized mass spectrum class for ImzML format with lazy loading capabilities.

    This class extends `Spectrum` to provide efficient handling of ImzML (Imaging Mass
    Spectrometry Markup Language) data. It implements lazy loading to minimize
    memory usage by loading spectrum data only when accessed.

    The class holds a reference to a spectrum reader and an index, loading the actual
    m/z and intensity arrays on-demand when the properties are first accessed.

    Attributes:
        _reader: Spectrum reader bound to the underlying ImzML parser implementation.
        _ibd_path (str): Path to the corresponding binary .ibd file.
        _index (int): Index of the spectrum within the ImzML file.

    Inherited Attributes (from `Spectrum`):
        coordinates (List[int]): Coordinates [x, y, z] of the spectrum
        x, y, z (int): Individual coordinate components
        sort_by_mz (bool): Whether the spectrum is sorted by m/z values.

    Properties:
        mz_list (array-like): Lazily loaded array of m/z values
        intensity (array-like): Lazily loaded array of intensity values

    Notes:
        - Data loading is deferred until the first property access.
        - Both `mz_list` and `intensity` are loaded together for efficiency.
        - Visualization/manipulation methods are inherited from `Spectrum`.
    """

    def __init__(self,
                 reader,
                 ibd_path, 
                 index: int, 
                 coordinates,
                 mz_list=None,
                 intensity=None,
                 sort_by_mz: bool = True,
                 swap_file_path: Optional[Path] = None):

        # ImzML reader part
        self._reader = reader
        self._ibd_path = ibd_path
        self._index = int(index)

        # Super init part
        super().__init__(mz_list=mz_list,
                         intensity=intensity,
                         coordinates=coordinates,
                         sort_by_mz=sort_by_mz,
                         swap_file_path=swap_file_path)

    def _resolve_data(self):
        """
        Internal hook to resolve/load data.
        Base implementation handles loading from swap file if available.
        """
        # first try swaper
        if self.swaper.offset >= 0:
            self.swaper.read_block()

        # then try imzml reader
        elif self._ibd_path is not None:
            with open(self._ibd_path, 'rb') as f:
                mz, intensity = self._reader.read_spectrum_from_file(f, self._index)
                self._mz_list = mz.copy()
                self._intensity = intensity.copy()

        # finally error
        else:
            logger.error("mz_list is None, can not resolve data.")
            raise ValueError("mz_list is None, can not resolve data.")
