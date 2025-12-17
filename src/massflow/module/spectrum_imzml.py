"""
Mass Spectrometry ImzML Spectrum Module for MassFlow Framework

This module provides a specialized class for handling ImzML format with lazy loading.

Classes:
    SpectrumImzML: Specialized class for handling ImzML format with lazy loading

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""

from typing import Optional
import numpy as np
from massflow.module.spectrum_base import SpectrumBaseModule
from massflow.logger import get_logger

logger = get_logger("SpectrumImzML")


class SpectrumImzML(SpectrumBaseModule):
    """
    Specialized mass spectrum class for ImzML format with lazy loading capabilities.

    This class extends SpectrumBaseModule to provide efficient handling of ImzML (Imaging Mass
    Spectrometry Markup Language) format data. It implements lazy loading to minimize
    memory usage by loading spectrum data only when accessed.

    The class holds a reference to a portable spectrum reader and an index, loading the actual
    m/z and intensity data on-demand when the properties are first accessed.

    Attributes:
        _reader: Portable spectrum reader bound to the underlying ImzMLParser.
        _ibd_path (str): Path to the corresponding binary .ibd file.
        _index (int): Index of the spectrum within the ImzML file.

    Inherited Attributes:
        coordinates (List[int]): 3D coordinates [x, y, z] of the spectrum
        x, y, z (int): Individual coordinate components
        sort_by_mz (bool): Flag inherited from base indicating if data is sorted by m/z values.

    Properties:
        mz_list (np.ndarray): Lazily loaded array of m/z values
        intensity (np.ndarray): Lazily loaded array of intensity values

    Note:
        - Data loading is deferred until first property access
        - Both mz_list and intensity are loaded together for efficiency
        - Inherits all visualization and manipulation methods from SpectrumBaseModule
    """

    def __init__(self,
                 reader,
                 ibd_path: str,
                 index: int,
                 coordinates,
                 mz_list: Optional[np.ndarray] = None,
                 intensity: Optional[np.ndarray] = None,
                 sort_by_mz: bool = True):
        """
        Initialize a lazy-loading ImzML spectrum.

        Args:
            reader: Portable spectrum reader from `ImzMLParser.portable_spectrum_reader()`.
            ibd_path (str): Path to the paired `.ibd` binary file.
            index (int): Spectrum index in the ImzML file.
            coordinates (PixelCoordinates | Sequence[int]): 3D coordinates [x, y, z].
            mz_list (np.ndarray | None): Preloaded m/z array; None to lazy-load.
            intensity (np.ndarray | None): Preloaded intensity array; None to lazy-load.
            sort_by_mz (bool): Whether the data is already sorted by m/z.
        """

        self._reader = reader
        self._ibd_path = ibd_path
        self._index = int(index)

        super().__init__(mz_list=mz_list,
                    intensity=intensity,
                    coordinates=coordinates,
                    sort_by_mz=sort_by_mz)

    @property
    def mz_list(self):
        """
        Lazily load and return m/z values for the spectrum.

        Returns:
            np.ndarray: Array of m/z values.

        Notes:
            - Triggers loading from the .ibd file via the portable spectrum reader on first access.
            - Also initializes `self._intensity` to keep arrays in sync.
        """
        if self._mz_list is None:
            _ = self.intensity  # Trigger loading of both mz and intensity
        return self._mz_list

    @mz_list.setter
    def mz_list(self, value):
        self._mz_list = value

    @property
    def intensity(self):
        """
        Lazily load and return intensity values for the spectrum.

        Returns:
            np.ndarray: Array of intensity values.

        Notes:
            - Accessing this property triggers loading of both m/z and intensity values from the .ibd file,
              so both arrays remain synchronized.
        """
        if self._intensity is None:
            with open(self._ibd_path, 'rb') as f:
                mz, intensity = self._reader.read_spectrum_from_file(f, self._index)
            self._intensity = intensity.copy()
            if self._mz_list is None:
                self._mz_list = mz.copy()
        return self._intensity

    @intensity.setter
    def intensity(self, value):
        self._intensity = value
