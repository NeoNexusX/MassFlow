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
from contextlib import nullcontext
from typing import Optional
from pyimzml.ImzMLParser import PortableSpectrumReader
from massflow.module.spectrum import Spectrum
from massflow.tools.stream_imzml_writer import ImzMLWriter
from massflow.tools.logger import get_logger

logger = get_logger("massflow.module.spectrum_imzml")

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
                 coordinates,
                 index: Optional[int] = None,
                 reader=None,
                 ibd_path=None,
                 mz_list=None,
                 intensity=None,
                 shared_mz_list=None,
                 sort_by_mz: bool = True,):

        # ImzML reader part
        self._reader :Optional[PortableSpectrumReader]= reader
        self._ibd_path = ibd_path
        self._index = index

        # Super init part
        super().__init__(mz_list=mz_list,
                         intensity=intensity,
                         coordinate=coordinates,
                         shared_mz_list=shared_mz_list,
                         sort_by_mz=sort_by_mz,)

    def resolve_data(self, **kwargs):
        """
        Internal hook to resolve/load data.
        Base implementation handles loading from swap file if available.
        """
        if self._intensity is not None:
            return

        if self._reader is None or self._ibd_path is None or self._index is None:
            logger.error("ImzML spectrum reader or ibd_path or index not properly set.")
            raise ValueError("ImzML spectrum reader or ibd_path or index not properly set.")

        file_handle = kwargs.get("file_handle")

        with nullcontext(file_handle) if file_handle is not None else open(self._ibd_path, "rb") as f:
            mz_list, intensity = self._reader.read_spectrum_from_file(f, self._index)

        if self.shared_mz_list is not None:
            self._mz_list = self.shared_mz_list
        elif self._mz_list is None:
            self._mz_list = mz_list.copy()

        self._intensity = intensity.copy()

    def swap_out2disk(self, writer: ImzMLWriter):
        """Write spectrum data to disk using ImzMLWriter."""

        if self._reader is None:

            if self.mz_list is None or self.intensity is None or self.coordinate is None:
                logger.error("Both mz_list and intensity cannot be None")
                raise ValueError("Both mz_list and intensity cannot be None")

            # add spectrum to writer
            writer.add_spectrum(self.mz_list, self.intensity, self.coordinate.get_tuple())

        # clear loaded data to save memory
        self.clear_data()
