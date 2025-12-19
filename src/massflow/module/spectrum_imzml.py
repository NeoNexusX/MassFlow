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

from massflow.module.spectrum import Spectrum
from pyimzml.ImzMLWriter import ImzMLWriter
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
                 index: int, 
                 coordinates,
                 reader=None,
                 ibd_path=None, 
                 mz_list=None,
                 intensity=None,
                 sort_by_mz: bool = True,):

        # ImzML reader part
        self._reader = reader
        self._ibd_path = ibd_path
        self._index = int(index)

        # Super init part
        super().__init__(mz_list=mz_list,
                         intensity=intensity,
                         coordinates=coordinates,
                         sort_by_mz=sort_by_mz,)

    def _resolve_data(self):
        """
        Internal hook to resolve/load data.
        Base implementation handles loading from swap file if available.
        """

        # then try imzml reader
        if self._ibd_path is not None:
            with open(self._ibd_path, 'rb') as f:
                mz, intensity = self._reader.read_spectrum_from_file(f, self._index)
                self._mz_list = mz.copy()
                self._intensity = intensity.copy()

        # finally error
        else:
            logger.error("mz_list is None, can not resolve data.")
            raise ValueError("mz_list is None, can not resolve data.")
        
    def swap_out2disk(self, writer: ImzMLWriter):
        """Write spectrum data to disk using ImzMLWriter."""

        if self._reader is None:

            if self.mz_list is None or self.intensity is None:
                logger.error("Both mz_list and intensity cannot be None")
                raise ValueError("Both mz_list and intensity cannot be None")

            assert len(self.mz_list) == len(self.intensity)

            f_mz = self.mz_list
            f_intensity = self.intensity

            coords_tuple = (self.x, self.y, self.z)

            writer.addSpectrum(f_mz, f_intensity, coords_tuple)

        else:
            self._mz_list = None
            self._intensity = None
