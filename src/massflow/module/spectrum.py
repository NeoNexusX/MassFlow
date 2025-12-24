import numpy as np
from typing import Union,Optional, Sequence
from massflow.module.pixel_coordinates import PixelCoordinates
from massflow.tools.stream_imzml_writer import ImzMLWriter
from massflow.logger import get_logger

logger = get_logger("spectrum")

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
        self.coordinates = PixelCoordinates(coordinates)

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
            _ = self.intensity  # trigger data load

        return self._mz_list

    @mz_list.setter
    def mz_list(self, value):
        if value is None or isinstance(value, np.ndarray):
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
                # sort the data by mz
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

    def get_coordinates(self):
        """
        Get the coordinates of the MSI data.

        Returns:
            PixelCoordinates: Object containing x, y, and z coordinates.
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

    def crop_range(self,
                   x_range: Sequence[float],
                   sort_by_mz: bool = True,
                   mode: str = "new",):
        """
        Crop spectrum within a given m/z range and return a shallow copy.

        Args:
            xr (Optional[Tuple[float, float] | List[float]]): Lower and upper m/z bounds.
            sort_by_mz (bool): Whether to sort by m/z prior to range search.
        """

        #make sure mz_list and intensity are not None
        if self.mz_list is not None and self.intensity is not None:

            # valid test
            min_mz = np.min(self.mz_list)
            max_mz = np.max(self.mz_list)
            if x_range[0] < min_mz or x_range[1] > max_mz:

                #make sure x_range is valid
                if x_range is not None and len(x_range) == 2:
                    mz_c = None
                    inten_c = None
                    
                    # without sort and dont want to sort
                    if not sort_by_mz and not self._sort_by_mz:
                        mask_xr = np.ones_like(self.mz_list, dtype=bool)
                        mask_xr &= (self.mz_list >= x_range[0]) & (self.mz_list <= x_range[1])

                        mz_c = self.mz_list[mask_xr]
                        inten_c = self.intensity[mask_xr]

                    # with sort and want to sort  or want sort but no sort
                    elif sort_by_mz:
                        # always set sort to true
                        self.sort_by_mz = True
                        # use searchsorted to find the indices
                        start_index = np.searchsorted(self.mz_list, x_range[0], side="left")
                        end_index = np.searchsorted(self.mz_list, x_range[1], side="right")

                        # build the data
                        mz_c = self.mz_list[start_index:end_index]
                        inten_c = self.intensity[start_index:end_index]
                else:
                    logger.info("xr is empty, can not crop by mz. use full mz range")
                    return self
            else:
                mz_c = np.array([-1])
                inten_c = np.array([-1])
        else:
            logger.error("mz_list is None, can not crop by mz. xr must be a tuple of two float numbers.")
            raise ValueError("mz_list is None, can not crop by mz. xr must be a tuple of two float numbers.")

        # cut part
        if mode == "new":
            # ---build a new object ：same class、copy attributes ---
            new_obj = self.__class__.__new__(self.__class__)
            new_obj.__dict__.update(self.__dict__)  # Copy all attributes
            new_obj.mz_list = mz_c
            new_obj.intensity = inten_c

            return new_obj

        elif mode == "update":
            self.mz_list = mz_c
            self.intensity = inten_c
            return self

        else:
            logger.error(f"mode {mode} is not supported. use 'new' or 'update'")
            raise ValueError(f"mode {mode} is not supported. use 'new' or 'update'")

    def is_sorted(self):
        # Check if the mz_list is sorted in ascending order.
        if self.mz_list is not None and self.intensity is not None:
            return np.all(self.mz_list[:-1] <= self.mz_list[1:])

        return False

    def clear_data(self):
        """
        Clear the loaded m/z and intensity data to free memory.
        """
        self.mz_list = None
        self.intensity = None
        self.coordinates = None

    @property
    def x(self):
        """
        Get the X coordinate from the underlying PixelCoordinates.

        Returns:
            int: The X coordinate.

        Raises:
            AttributeError: If 'coordinates' is not initialized.
        """
        if self.coordinates is None:
            raise AttributeError("coordinates is not initialized.")
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
        if self.coordinates is None:
            raise AttributeError("coordinates is not initialized.")
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
        if self.coordinates is None:
            raise AttributeError("coordinates is not initialized.")
        return self.coordinates.z
    
    def swap_out2disk(self, writer: ImzMLWriter):
        pass