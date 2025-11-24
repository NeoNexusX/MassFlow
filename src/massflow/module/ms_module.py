"""
Mass Spectrometry Module for MassFlow Framework

This module provides core classes and functionality for handling mass spectrometry (MS) data,
particularly for Mass Spectrometry Imaging (MSI) applications. It includes support for
lazy loading, efficient data management, and visualization capabilities.

Classes:
    SpectrumBaseModule: Base class for mass spectrum data with coordinates
    SpectrumImzML: Specialized class for handling ImzML format with lazy loading
    MS: Collection class for managing multiple mass spectra

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""

from typing import List, Tuple, Union, Optional, Sequence,Literal
import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser
from massflow.logger import get_logger
from massflow.module.meta_data import MetaDataFileBase
logger = get_logger("ms_module")


class PixelCoordinates:
    """
    Coordinate holder for MSI pixels with optional zero-based view.

    This class stores raw `x`, `y`, `z` coordinates and provides a `zero_based`
    switch that, when enabled, makes the public getters return values offset by `-1`
    (useful when converting 1-based coordinates to 0-based indices).

    Attributes:
        x (int): X coordinate (getter returns adjusted value if `zero_based` is True).
        y (int): Y coordinate (getter returns adjusted value if `zero_based` is True).
        z (int): Z coordinate (getter returns adjusted value if `zero_based` is True).
        zero_based (bool): Whether to expose coordinates in zero-based form.
    """

    def __init__(self, x: int, y: int, z: int, zero_based: bool = False):
        """
        Initialize pixel coordinates with optional zero-based adjustment.

        Args:
            x (int): X coordinate value.
            y (int): Y coordinate value.
            z (int): Z coordinate value.
            zero_based (bool): If True, getters will return (value - 1)
                to convert 1-based coordinates to 0-based.

        Returns:
            None

        Raises:
            TypeError: If any of x, y, or z cannot be converted to int.
        """

        # internal storage
        self._zero_based = None
        # Store raw coordinates internally; getters apply adjustment when required.
        self._x : Optional[int] = None
        self._y : Optional[int] = None
        self._z : Optional[int] = None

        self.x = x
        self.y = y
        self.z = z
        self.zero_based = zero_based

    @property
    def x(self) -> Optional[int]:
        """
        Get the X coordinate.

        Returns:
            int: X coordinate, adjusted by `-1` when `zero_based` is True.
        """
        if self._x is not None:
            return self._x - 1 if self.zero_based else self._x
        else:
            logger.error("X coordinate is not set.")
            return None

    @x.setter
    def x(self, value: int):
        self._x = value

    @property
    def y(self) -> Optional[int]:
        """
        Get the Y coordinate.

        Returns:
            int: Y coordinate, adjusted by `-1` when `zero_based` is True.
        """
        if self._y is  not None:
            return self._y - 1 if self.zero_based else self._y
        else:
            logger.error("Y coordinate is not set.")
            return None

    @y.setter
    def y(self, value: int):
        self._y = value

    @property
    def z(self) -> Optional[int]:
        """
        Get the Z coordinate.

        Returns:
            int: Z coordinate, adjusted by `-1` when `zero_based` is True.
        """
        if self._z is not None:
            return self._z - 1 if self.zero_based else self._z
        else:
            logger.error("Z coordinate is not set.")
            return None

    @z.setter
    def z(self, value: int):
        self._z = value

    @property
    def zero_based(self) -> Optional[bool]:
        """
        Flag indicating whether getters return zero-based coordinates.

        Returns:
            bool: True if coordinates are zero-based; False otherwise.
        """
        if self._zero_based is not None:
            return self._zero_based
        else:
            logger.error("zero_based flag is not set.")
            return None

    @zero_based.setter
    def zero_based(self, value: bool):
        """
        Set the flag controlling zero-based view of coordinates.

        Args:
            value (bool): True for zero-based; False for one-based.
        """
        self._zero_based = bool(value)

    def __eq__(self, other: object) -> bool:
        """
        Compare equality based on adjusted coordinates.

        Args:
            other (PixelCoordinates): Another PixelCoordinates instance to compare.

        Returns:
            bool: True if all adjusted coordinates match; False otherwise.

        Raises:
            TypeError: If `other` is not an instance of PixelCoordinates.
        """
        if not isinstance(other, PixelCoordinates):
            raise TypeError("other must be a PixelCoor")
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self):
        """
        Compute hash based on adjusted coordinates.

        Returns:
            int: Hash value derived from adjusted (x, y, z).

        Raises:
            None
        """
        return hash((self.x, self.y, self.z))

    def lefter(self, other: "PixelCoordinates") -> Optional[bool]:
        """
        Determine if this pixel is to the left of another pixel.

        Args:
            other (PixelCoordinates): The other pixel to compare against.

        Returns:
            bool: True if this pixel's X is less than the other's X; otherwise False.
        """
        if self.x is None or other.x is None:
            logger.error("X coordinate is None, cannot compare.")
            return None
        else:
            return self.x < other.x

    def righter(self, other: "PixelCoordinates") -> Optional[bool]:
        """
        Determine if this pixel is to the right of another pixel.

        Args:
            other (PixelCoordinates): The other pixel to compare against.

        Returns:
            bool: True if this pixel's X is greater than the other's X; otherwise False.
        """
        if self.x is None or other.x is None:
            logger.error("X coordinate is None, cannot compare.")
            return None
        else:
            return self.x > other.x

    def upper(self, other: "PixelCoordinates") -> Optional[bool]:
        """
        Determine if this pixel is above another pixel.

        Args:
            other (PixelCoordinates): The other pixel to compare against.

        Returns:
            bool: True if this pixel's Y is greater than the other's Y; otherwise False.
        """
        if self.y is None or other.y is None:
            logger.error("Y coordinate is None, cannot compare.")
            return None
        else:
            return self.y > other.y

    def __repr__(self) -> str:
        """
        Return a readable string representation of the coordinates.

        Returns:
            str: A string formatted as '(x, y, z)'.

        Raises:
            None
        """
        return f"({self.x}, {self.y}, {self.z})"

    def __len__(self) -> int:
        """
        Return the number of dimensions.

        Returns:
            int: Returns 3 when all coordinates are initialized; otherwise 0.

        Raises:
            None
        """
        if self._x is None or self._y is None or self._z is None:
            return 0
        return 3


class SpectrumBaseModule:
    """
    Base class for mass spectrum data with spatial coordinates.

    This class represents a single mass spectrum with associated m/z values, intensities,
    and 3D spatial coordinates. It supports lazy loading for memory efficiency and provides
    basic operations for mass spectrometry data manipulation and visualization.

    Attributes:
        coordinates (List[int]): 3D coordinates [x, y, z] of the spectrum
        x (int): X coordinate
        y (int): Y coordinate
        z (int): Z coordinate
        sorted_by_mz (bool): Flag indicating if data is sorted by m/z values

    Properties:
        mz_list (np.ndarray): Array of m/z values
        intensity (np.ndarray): Array of intensity values corresponding to m/z values

    Note:
        - Coordinates must be a list of exactly three or two integers
        - m/z and intensity arrays must have the same length
        - Lazy loading is supported through None initialization of mz_list and intensity
    """

    def __init__(
        self,
        mz_list: Optional[np.ndarray],
        intensity: Optional[np.ndarray],
        coordinates: Union[PixelCoordinates, Sequence[int]],
        sorted_by_mz_fun: bool = False,
    ):
        """
        Initialize a mass spectrum with m/z values, intensities, and coordinates.

        Args:
            mz_list (Optional[np.ndarray]): Array of m/z values. Can be None for lazy loading.
            intensity (Optional[np.ndarray]): Array of intensity values. Can be None for lazy loading.
            coordinates (Union[PixelCoordinates, List[int], Tuple[int,int,int]]): Spatial coordinates;
                if a list/tuple is provided, it will be converted to a PixelCoordinates object.
            sorted_by_mz_fun (bool, optional): Whether the data is already sorted by m/z. Defaults to False.

        Raises:
            AssertionError: If `coordinates` does not have length 3.
            TypeError: If `coordinates` is neither `PixelCoordinates` nor a list/tuple of three ints.
        """

        assert len(coordinates) == 3, "Coordinates must be a list of three integers."

        # Lazy loading
        self._mz_list = mz_list
        self._intensity = intensity

        # Initialize coordinates to PixelCoordinates instance
        if isinstance(coordinates, PixelCoordinates):
            self.coordinates = coordinates
        elif isinstance(coordinates, (list, tuple)) and len(coordinates) == 3:
            x, y, z = coordinates
            self.coordinates = PixelCoordinates(int(x), int(y), int(z))
        else:
            raise TypeError(
                "coordinates must be PixelCoordinates or a list/tuple of three ints."
            )

        self.sorted_by_mz = False
        if sorted_by_mz_fun:
            self.sort_by_mz()

    # lazy load properties
    @property
    def mz_list(self):
        """
        Get the mz values of the MSI data.

        Returns:
            np.ndarray: An array of mz values.
        """
        return self._mz_list

    @mz_list.setter
    def mz_list(self, value):
        self._mz_list = value

    @property
    def intensity(self):
        """
        Get the intensity values of the MSI data.

        Returns:
            np.ndarray: An array of intensity values.
        """
        return self._intensity

    @intensity.setter
    def intensity(self, value):
        self._intensity = value

    def get_coordinates(self) -> PixelCoordinates:
        """
        Get the coordinates of the MSI data.

        Returns:
            PixelCoordinates: Object containing x, y, and z coordinates.

        Raises:
            None
        """
        return self.coordinates

    def __len__(self):
        """
        Return the number of m/z peaks in the spectrum.

        Returns:
            int: Number of peaks (length of mz_list array)
        """
        if self.mz_list is not None:
            return len(self.mz_list)

    def __eq__(self, other):
        """
        Check equality based on coordinates.

        Two SpectrumBaseModule instances are considered equal if they have the same coordinates.

        Args:
            other: Object to compare with

        Returns:
            bool: True if coordinates are equal, False otherwise
        """
        if not isinstance(other, SpectrumBaseModule):
            return False
        return self.coordinates == other.coordinates

    def __getitem__(self, index):
        """
        Get m/z and intensity values at specified index.

        Args:
            index (int): Index of the peak to retrieve

        Returns:
            Tuple[float, float]: (m/z, intensity) values at the given index

        Raises:
            IndexError: If index is out of range
        """
        if self.mz_list is None or self.intensity is None:
            logger.error("mz_list or intensity is None, can not get item.")
            raise IndexError("mz_list or intensity is None, can not get item.")
        return self.mz_list[index], self.intensity[index]

    def sort_by_mz(self):
        """
        Sort the mz_list and intensity arrays by m/z values in ascending order.

        This method sorts both arrays simultaneously to maintain correspondence
        between m/z values and their intensities. The operation is performed
        in-place and updates the `sorted_by_mz` flag.

        Note:
            - Only sorts if data is not already sorted (`sorted_by_mz` is False)
            - Logs a warning if mz_list or intensity is None
            - After sorting, `sorted_by_mz` flag is set to True
        """
        if (self.sorted_by_mz is False and self.mz_list is not None and self.intensity is not None):

            sorted_indices = np.argsort(self.mz_list)
            self._mz_list = self.mz_list[sorted_indices]
            self._intensity = self.intensity[sorted_indices]
            self.sorted_by_mz = True

        elif self.mz_list is None or self.intensity is None:
            logger.warning("mz_list or intensity is None, can not sort by mz.")

    def crop_range(
        self,
        xr: Optional[Sequence[float]] = None,
        sort_by_mz: bool = True,
    ):
        """
        Crop spectrum within a given m/z range and return a shallow copy.

        Args:
            xr (Optional[Tuple[float, float] | List[float]]): Lower and upper m/z bounds.
            sort_by_mz (bool): Whether to sort by m/z prior to range search.

        Returns:
            SpectrumBaseModule: New spectrum object sharing attributes with sliced arrays.

        Raises:
            ValueError: If `mz_list` is None or `xr` is invalid.
        """
        if self.mz_list is not None and self.intensity is not None:
            if xr is not None and len(xr) == 2:
                mz_c = None
                inten_c = None
                # without sort and dont want to sort
                if not self.sorted_by_mz and not sort_by_mz:
                    mask_xr = np.ones_like(self.mz_list, dtype=bool)

                    if xr is not None and self.mz_list is not None:
                        mask_xr &= (self.mz_list >= xr[0]) & (self.mz_list <= xr[1])
                    else:
                        logger.error("xr or mz_list is None, can not crop by mz.")
                        raise ValueError("xr or mz_list is None, can not crop by mz.")

                    mz_c = self.mz_list[mask_xr]
                    inten_c = self.intensity[mask_xr]

                # with sort and want to sort  or want sort but no sort
                elif sort_by_mz:
                    if not self.sorted_by_mz:
                        self.sort_by_mz()
                    start_index = np.searchsorted(self.mz_list, xr[0], side="left")
                    end_index = np.searchsorted(self.mz_list, xr[1], side="right")

                    # build the data
                    mz_c = self.mz_list[start_index:end_index]
                    inten_c = self.intensity[start_index:end_index]
            else:
                logger.info("xr is empty, can not crop by mz. use full mz range")
                return self
        else:
            logger.error(
                "mz_list is None, can not crop by mz. xr must be a tuple of two float numbers."
            )
            raise ValueError(
                "mz_list is None, can not crop by mz. xr must be a tuple of two float numbers."
            )

        # ---build a new object ：same class、copy attributes ---
        new_obj = self.__class__.__new__(self.__class__)
        new_obj.__dict__.update(self.__dict__)  # Copy all attributes
        new_obj.mz_list = mz_c
        new_obj.intensity = inten_c
        return new_obj

    @property
    def x(self):
        """
        Get the X coordinate from the underlying PixelCoordinates.

        Returns:
            int: The X coordinate.

        Raises:
            AttributeError: If 'coordinates' is not initialized.
        """
        return self.coordinates.x

    @property
    def y(self):
        """
        Get the Y coordinate from the underlying PixelCoordinates.

        Returns:
            int: The Y coordinate.

        Raises:
            AttributeError: If 'coordinates' is not initialized.
        """
        return self.coordinates.y

    @property
    def z(self):
        """
        Get the Z coordinate from the underlying PixelCoordinates.

        Returns:
            int: The Z coordinate.

        Raises:
            AttributeError: If 'coordinates' is not initialized.
        """
        return self.coordinates.z


class SpectrumImzML(SpectrumBaseModule):
    """
    Specialized mass spectrum class for ImzML format with lazy loading capabilities.

    This class extends SpectrumBaseModule to provide efficient handling of ImzML (Imaging Mass
    Spectrometry Markup Language) format data. It implements lazy loading to minimize
    memory usage by loading spectrum data only when accessed.

    The class holds a reference to an ImzMLParser and an index, loading the actual
    m/z and intensity data on-demand when the properties are first accessed.

    Attributes:
        _parser (ImzMLParser): Parser instance for reading ImzML data
        _index (int): Index of the spectrum within the ImzML file

    Inherited Attributes:
        coordinates (List[int]): 3D coordinates [x, y, z] of the spectrum
        x, y, z (int): Individual coordinate components
        sorted_by_mz_fun (bool): Flag indicating if data is sorted by m/z values

    Properties:
        mz_list (np.ndarray): Lazily loaded array of m/z values
        intensity (np.ndarray): Lazily loaded array of intensity values

    Note:
        - Data loading is deferred until first property access
        - Both mz_list and intensity are loaded together for efficiency
        - Inherits all visualization and manipulation methods from SpectrumBaseModule
    """

    def __init__(self, parser: ImzMLParser, index: int, coordinates,mz_list=None,intensity=None):

        super().__init__(mz_list=mz_list, intensity=intensity, coordinates=coordinates)
        self._parser = parser
        self._index = int(index)

    @property
    def mz_list(self):
        """
        Lazily load and return m/z values for the spectrum.

        Returns:
            np.ndarray: Array of m/z values.

        Notes:
            - Triggers loading from the underlying ImzML parser on first access.
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
            - Accessing this property ensures m/z values are loaded first,
              so both arrays remain synchronized.
        """
        if self._intensity is None:
            mz, intensity = self._parser.getspectrum(self._index)
            self._intensity = intensity

            if self._mz_list is None:
                self._mz_list = mz

        return self._intensity

    @intensity.setter
    def intensity(self, value):
        self._intensity = value


class MS:
    """
    Collection class for managing multiple mass spectra with coordinate-based indexing.

    This class serves as a container and manager for multiple SpectrumBaseModule instances,
    providing efficient storage, retrieval, and manipulation of mass spectrometry data
    organized by 3D spatial coordinates. It supports both sequential and coordinate-based
    access patterns.

    The class maintains two internal data structures:
    - A queue (_queue) for sequential access and iteration
    - A nested dictionary (_coordinate_index) for fast coordinate-based lookup

    Attributes:
        _queue (List[SpectrumBaseModule]): Sequential list of all spectra
        _coordinate_index (Dict): Nested dictionary mapping coordinates to spectra
                                 Structure: {z: {x: {y: MSBaseModule}}}

    Indexing Methods:
        - ms[index]: Access by sequential index
        - ms[x, y]: Access by coordinates (z defaults to 0)
        - ms[x, y, z]: Access by full 3D coordinates
        - ms[x, y, z] = spectrum: Direct assignment

    Note:
        - Coordinates are automatically managed and indexed
        - Supports both 2D (x, y) and 3D (x, y, z) coordinate systems
        - Efficient lookup performance through coordinate indexing
        - Thread-safe for read operations
    """

    def __init__(self):
        """
        Initialize an empty MS collection.

        Creates empty internal data structures for storing and indexing mass spectra.
        No parameters are required for initialization.
        """
        self.meta : Optional[MetaDataFileBase] = None
        self._queue = []
        self._coordinate_index = {}  # Mapping from coordinates to MSBaseModule
        self._shared_mz_list = None

    @property
    def shared_mz_list(self):
        """
        Get or set the flag for using common m/z list across spectra.

        Returns:
            bool: True if a common m/z list is used; False otherwise.
        """
        return self._shared_mz_list
    @shared_mz_list.setter
    def shared_mz_list(self, value):
        if self.meta:
            if self.meta.continuous is True:
                self._shared_mz_list = value
            else:
                logger.warning("Meta data indicates non-continuous data; co_use_mz_list set to False.")

    @property
    def coordinate_index(self):
        """
        Get the nested dictionary of coordinate indices.

        Returns:
            Dict: Nested dictionary mapping coordinates to spectra
                Structure: {z: {x: {y: MSBaseModule}}}
        """
        return self._coordinate_index

    def add_spectrum(self, spectrum: SpectrumBaseModule):
        """
        Add or update a mass spectrum with coordinate indexing.

        This method ensures the spectrum is accessible via both the sequential queue
        and the coordinate index. If a spectrum already exists at the given coordinates,
        it will be updated (replaced) in both the coordinate index and the queue to
        avoid duplicates.

        Args:
            spectrum (SpectrumBaseModule): Mass spectrum to add or update.


        Notes:
            - Automatically extracts coordinates from the spectrum.
            - Creates nested dictionary structure if coordinates don't exist.
            - Updates the existing spectrum in place when coordinates already exist.
        """
        self.update_spectrum_with_coord(spectrum = spectrum)

    def get_spectrum(self, x: int, y: int, z: int = 0) -> SpectrumBaseModule:
        """
        Retrieve a mass spectrum by its 3D coordinates.

        Args:
            x (int): X coordinate
            y (int): Y coordinate
            z (int, optional): Z coordinate. Defaults to 0.

        Returns:
            SpectrumBaseModule: Mass spectrum at the specified coordinates

        Raises:
            KeyError: If no spectrum exists at the specified coordinates
        """
        # Check if the coordinates exist in the index
        if (
            z not in self._coordinate_index
            or x not in self._coordinate_index[z]
            or y not in self._coordinate_index[z][x]
        ):
            logger.error(
                f"No spectrum found at coordinates ({x}, {y}, {z})\r\n"
                f"plot mask to see the aviliable location, use ms.coordinate_index to see all the coordinates."
            )
            raise KeyError(f"No spectrum found at coordinates ({x}, {y}, {z})")
        return self._coordinate_index[z][x][y]

    def update_spectrum_with_coord(self,spectrum: SpectrumBaseModule,x=-1,y=-1,z=-1):
        """
        Add or replace a spectrum in the coordinate index and queue.

        Args:
            spectrum (SpectrumBaseModule): Spectrum to insert or update.
            x (int): Optional override for X coordinate; uses `spectrum.x` when -1.
            y (int): Optional override for Y coordinate; uses `spectrum.y` when -1.
            z (int): Optional override for Z coordinate; uses `spectrum.z` when -1.

        Returns:
            None

        Raises:
            ValueError: If updating an existing spectrum fails to synchronize in the queue.
        """
        # Sync zero_based setting if meta is available
        if self.meta :
            spectrum.coordinates.zero_based = self.meta.coordinates_zero_based

        # get locations
        x = spectrum.x if x == -1 else x
        y = spectrum.y if y == -1 else y
        z = spectrum.z if z == -1 else z

        if (z in self._coordinate_index and
            x in self._coordinate_index[z] and
            y in self._coordinate_index[z][x]
            ):
            existing = self._coordinate_index[z][x][y]
            self._coordinate_index[z][x][y] = spectrum
            try:
                idx = self._queue.index(existing)
                self._queue[idx] = spectrum
            except ValueError as exc:
                logger.error(f"Failed to update spectrum at coordinates ({x}, {y}, {z})")
                raise ValueError(f"Failed to update spectrum at coordinates ({x}, {y}, {z})") from exc
        else:
            if z not in self._coordinate_index:
                self._coordinate_index[z] = {}
            if x not in self._coordinate_index[z]:
                self._coordinate_index[z][x] = {}
            self._coordinate_index[z][x][y] = spectrum
            self._queue.append(spectrum)

    def update_spectrum_with_index(self,index:int,spectrum: SpectrumBaseModule):
        """
        Replace a spectrum at a sequential index and synchronize the coordinate index.

        Args:
            index (int): Position in the internal queue to replace.
            spectrum (SpectrumBaseModule): Spectrum to assign at the given index.

        Returns:
            None

        Raises:
            IndexError: If `index` is negative or out of range.
            ValueError: If the spectrum's coordinates do not match the existing entry at `index`.
        """
        if index < 0 or isinstance(index,int) is False :
            logger.error(f"Index {index} out of range.")
            raise IndexError(f"Index {index} out of range.")

        elif index >=0 and index < len(self._queue):
            # Sync zero_based setting if meta is available
            if self.meta:
                spectrum.coordinates.zero_based = self.meta.coordinates_zero_based

            # check if coordinates match
            if self._queue[index] == spectrum:
                #update queue
                self._queue[index] = spectrum
                self._coordinate_index[spectrum.z][spectrum.x][spectrum.y] = spectrum
            else:
                logger.error(f"Spectrum at index {index} does not match coordinates ({spectrum.x}, {spectrum.y}, {spectrum.z}).")
                raise ValueError(f"Spectrum at index {index} does not match coordinates ({spectrum.x}, {spectrum.y}, {spectrum.z}).")
        else:
            logger.warning(f"Index {index} out of range {len(self._queue)}. Update spectrum with coordinates.")
            self.update_spectrum_with_coord(spectrum = spectrum)

    def __getitem__(self,
                    key: Union[int,Tuple[int, int, int], Tuple[int, int], slice]
    ) -> Union[SpectrumBaseModule, List[SpectrumBaseModule]]:
        """
        Retrieve mass spectrum using flexible indexing methods.

        Supports multiple indexing patterns for convenient access to mass spectra:
        - Sequential indexing: ms[index]
        - 2D coordinates: ms[x, y] (z defaults to 0)
        - 3D coordinates: ms[x, y, z]
        - Slice: ms[a:b:c] returns a list of spectra from the internal queue

        Args:
            key (Union[int, Tuple[int, int], Tuple[int, int, int], slice]):
                Index, coordinates, or slice for spectrum retrieval

        Returns:
            Union[SpectrumBaseModule, List[SpectrumBaseModule]]: Single spectrum for index/coordinates,
            or a list of spectra for slice access

        Raises:
            TypeError: If key format is not supported
            KeyError: If coordinates don't exist in the collection
            IndexError: If sequential index is out of range
        """
        if isinstance(key, int):
            # Return the item from the queue by index
            return self._queue[key]
        elif isinstance(key, tuple):
            if len(key) == 3:
                x, y, z = key
                return self._coordinate_index[z][x][y]
            elif len(key) == 2:
                x, y = key
                z = 0 if 0 in self._coordinate_index else next(iter(self._coordinate_index))
                return self._coordinate_index[z][x][y]
        elif isinstance(key, slice):
            # Support slice access on the internal sequential queue
            return self._queue[key]
        else:
            logger.error("Index must be in tuple format, like [x, y, z] or [x, y]")
            raise TypeError(
                "Index must be in tuple format, like [x, y, z] or [x, y]"
            )

    def __setitem__(
        self,
        key: Union[int, Tuple[int, int, int], Tuple[int, int]],
        spectrum: SpectrumBaseModule,
    ):
        """
        Assign or replace a spectrum by index or coordinates with bounds checking and index sync.

        Args:
            key (Union[int, Tuple[int, int], Tuple[int, int, int]]): Sequential index or target coordinates.
            spectrum (SpectrumBaseModule): Spectrum to assign.

        Returns:
            None

        Raises:
            IndexError: If coordinate tuple length is not 2 or 3, or sequential index is out of range.
            ValueError: If assigning by index but the spectrum's coordinates do not match the existing entry.
        """
        if isinstance(key,int):
            self.update_spectrum_with_index(key,spectrum)
        elif len(key) == 3:
            x, y, z = key
            self.update_spectrum_with_coord(spectrum,x,y,z)
        elif len(key) == 2:
            x, y = key
            z = next(iter(self._coordinate_index))
            self.update_spectrum_with_coord(spectrum,x,y,z)
        else:
            logger.error("Coordinates must be 2 or 3 integers")
            raise IndexError("Coordinates must be 2 or 3 integers")

    def __len__(self):
        """
        Return the number of spectra in the collection.

        Returns:
            int: Total number of mass spectra in the collection
        """
        return len(self._queue)

    def __iter__(self):
        """
        Return an iterator over all spectra in the collection.

        Allows iteration through all mass spectra in the order they were added.

        Returns:
            Iterator[SpectrumBaseModule]: Iterator over mass spectra
        """
        return iter(self._queue)

    def plot_ms_mask(
        self,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (8, 8),
        dpi: int = 300,
        origin: Literal["upper", "lower"] = "lower",
        cmap: str = "Greys",
    ):
        """
        Plot the occupancy mask stored in metadata.

        Args:
            save_path (Optional[str]): Path to save the figure; show if None.
            figsize (Tuple[int, int]): Matplotlib figure size (width, height). Defaults to (8, 8).
            dpi (int): Dots-per-inch for saving. Defaults to 300.
            origin (str): Image origin for `imshow`, either 'upper' or 'lower'. Defaults to 'upper'.
            cmap (str): Colormap for rendering the mask. Defaults to 'Greys'.

        Returns:
            None

        Raises:
            ValueError: If metadata or mask is missing.

        Notes:
            - Directly visualizes `self.meta.mask` without recomputing.
            - Ensure `self.meta.mask` shape matches `(max_count_of_pixels_y, max_count_of_pixels_x)`.
        """
        # Validate metadata and mask
        if self.meta is None:
            logger.error("MS meta data is required to plot mask.")
            raise ValueError("MS meta data is required to plot mask.")

        mask = getattr(self.meta, "mask", None)
        if mask is None:
            logger.error("Meta mask is None. Create mask before plotting.")
            raise ValueError("Meta mask is None. Create mask before plotting.")

        # Plot the mask
        plt.figure(figsize=figsize)
        plt.imshow(mask, cmap=cmap, origin=origin)
        plt.title("MS Occupancy Mask")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=dpi)
        else:
            plt.show()

    def set_coordinates_zero_based(self, zero_based: bool):
        """
        Synchronize zero-based coordinate setting across metadata and all spectra, rebuilding index.

        Args:
            zero_based (bool): True to use zero-based coordinates; False for one-based.

        Raises:
            ValueError: If metadata is missing on the MS instance.

        Notes:
            - Updates `self.meta.coordinates_zero_based` and toggles each spectrum's
              `PixelCoordinates.zero_based` property.
            - Rebuilds the coordinate index to reflect new x/y/z values under the flag.
        """
        if self.meta is None:
            logger.error("MS meta data is required to synchronize coordinates.")
            raise ValueError("MS meta data is required to synchronize coordinates.")

        prev_flag = self.meta.coordinates_zero_based
        new_flag = zero_based

        # If flag changed, update all coordinates and rebuild index
        if prev_flag != new_flag:
            # Update metadata flag
            self.meta.coordinates_zero_based = new_flag
            self.meta.update_meta()
            # Rebuild coordinate index (clear in-place to avoid stale external references)
            self._coordinate_index.clear()
            # update min/max pixel values in metadata

            for spectrum in self._queue:
                spectrum.coordinates.zero_based = new_flag
                x, y, z = spectrum.x, spectrum.y, spectrum.z
                if z not in self._coordinate_index:
                    self._coordinate_index[z] = {}
                if x not in self._coordinate_index[z]:
                    self._coordinate_index[z][x] = {}
                self._coordinate_index[z][x][y] = spectrum
