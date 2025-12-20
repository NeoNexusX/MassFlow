import numpy as np
from typing import Union,Optional, Sequence
from massflow.tools.stream_imzml_writer import ImzMLWriter
from massflow.logger import get_logger

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

    def get_tuple(self) -> Optional[tuple]:
        """
        Get coordinates as a tuple.

        Returns:
            tuple: (x, y, z) coordinates, adjusted by `-1` if `zero_based` is True.
        """
        if self.x is not None and self.y is not None and self.z is not None:
            return (self.x, self.y, self.z)
        else:
            logger.error("One or more coordinates are not set.")
            return None

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
            raise TypeError("other must be a PixelCoordinates instance")
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


class Spectrum:
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
        sort_by_mz (bool): Flag indicating if data is sorted by m/z values

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
        sort_by_mz: bool = True,
    ):
        """
        Initialize a mass spectrum with m/z values, intensities, and coordinates.

        Args:
            mz_list (Optional[np.ndarray]): Array of m/z values. Can be None for lazy loading.
            intensity (Optional[np.ndarray]): Array of intensity values. Can be None for lazy loading.
            coordinates (Union[PixelCoordinates, List[int], Tuple[int,int,int]]): Spatial coordinates;
                if a list/tuple is provided, it will be converted to a PixelCoordinates object.
            sort_by_mz (bool, optional): Whether the data is already sorted by m/z. Defaults to True.

        Raises:
            AssertionError: If `coordinates` does not have length 3.
            TypeError: If `coordinates` is neither `PixelCoordinates` nor a list/tuple of three ints.
        """

        assert len(coordinates) == 3, "Coordinates must be a list of three integers."

        # Lazy loading
        self._mz_list = mz_list
        self._intensity = intensity

        # Initialize coordinates to PixelCoordinates instance
        # need to update
        if isinstance(coordinates, PixelCoordinates):
            self.coordinates = coordinates
        elif isinstance(coordinates, (list, tuple)) and len(coordinates) == 3:
            x, y, z = coordinates
            self.coordinates = PixelCoordinates(int(x), int(y), int(z))
        else:
            raise TypeError(
                "coordinates must be PixelCoordinates or a list/tuple of three ints."
            )


        # sort part , do not use _sort_by_mz
        self._sort_by_mz = sort_by_mz

    def _resolve_data(self):
        """
        Internal hook to resolve/load data.
        Base implementation handles loading from swap file if available.
        """
        pass

    # lazy load properties
    @property
    def mz_list(self):
        """
        Get the mz values of the MSI data.

        Returns:
            np.ndarray: An array of mz values.
        """
        if self._mz_list is None:
            self._resolve_data()

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
        if self._intensity is None:
            self._resolve_data()

        return self._intensity

    @intensity.setter
    def intensity(self, value):
        self._intensity = value

    @property
    def sort_by_mz(self):
        return self._sort_by_mz

    @sort_by_mz.setter
    def sort_by_mz(self, value: bool = True):
        if value :
            # make sure both mz_list and intensity are not None
            if self._sort_by_mz is False and self.mz_list is not None and self.intensity is not None:
                sorted_indices = np.argsort(self.mz_list)
                self._mz_list = self.mz_list[sorted_indices]
                self._intensity = self.intensity[sorted_indices]
                self._sort_by_mz = value

            elif self._sort_by_mz is True:
                logger.info("Data is already sorted by mz.")

            elif self.mz_list is None or self.intensity is None:
                logger.warning("mz_list or intensity is None, can not sort by mz.")
        else:
            self._sort_by_mz = value

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
        else:
            return 0

    def __eq__(self, other):
        """
        Check equality based on coordinates.

        Two Spectrum instances are considered equal if they have the same coordinates.

        Args:
            other: Object to compare with

        Returns:
            bool: True if coordinates are equal, False otherwise
        """
        if not isinstance(other, Spectrum):
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
        #make sure mz_list and intensity are not None
        if self.mz_list is not None and self.intensity is not None:
            #make sure xr is valid
            if xr is not None and len(xr) == 2:
                mz_c = None
                inten_c = None
                
                # without sort and dont want to sort
                if not self._sort_by_mz and not sort_by_mz:
                    mask_xr = np.ones_like(self.mz_list, dtype=bool)

                    if self.mz_list is not None:
                        mask_xr &= (self.mz_list >= xr[0]) & (self.mz_list <= xr[1])
                    else:
                        logger.error("xr or mz_list is None, can not crop by mz.")
                        raise ValueError("xr or mz_list is None, can not crop by mz.")

                    mz_c = self.mz_list[mask_xr]
                    inten_c = self.intensity[mask_xr]

                # with sort and want to sort  or want sort but no sort
                elif sort_by_mz:

                    if not self.sort_by_mz:
                        self.sort_by_mz = True

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

    def is_sorted(self):
        # Check if the mz_list is sorted in ascending order.
        if self.mz_list is not None and self.intensity is not None:
            return np.all(self.mz_list[:-1] <= self.mz_list[1:])
        return False

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
    
    def swap_out2disk(self, writer: ImzMLWriter):
        pass