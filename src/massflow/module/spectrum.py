import copy
from typing import Union, Optional, Sequence
import numpy as np
from massflow.module.pixel_coordinates import PixelCoordinates
from massflow.tools.stream_imzml_writer import ImzMLWriter
from massflow.tools.logger import get_logger

logger = get_logger("massflow.module.spectrum")


class Spectrum:
    """
    Base class for mass spectrum data with spatial coordinates.

    Attributes:
        coordinate (PixelCoordinates): 3D spatial coordinates.
        mz_list (np.ndarray): Array of m/z values.
        intensity (np.ndarray): Array of intensity values.
    """

    def __init__(
        self,
        mz_list: Optional[np.ndarray],
        intensity: Optional[np.ndarray],
        coordinate: Union[PixelCoordinates, Sequence[int]],
        sort_by_mz: bool = True,
        shared_mz_list: Optional[np.ndarray] = None,
    ):
        """
        Initialize a mass spectrum.

        Args:
            mz_list: Array of m/z values (None for lazy loading).
            intensity: Array of intensity values (None for lazy loading).
            coordinate: Spatial coordinates (x, y, z).
            sort_by_mz: Whether data is sorted by m/z.
            shared_mz_list: Optional shared m/z array for continuous data.
        """

        assert len(coordinate) == 3, "Coordinates must be a list of three integers."

        # Lazy loading
        self.shared_mz_list = shared_mz_list
        self._mz_list = self.shared_mz_list if mz_list is None else mz_list
        self._intensity = intensity

        # Initialize coordinates to PixelCoordinates instance
        self.coordinate = PixelCoordinates(coordinate)

        # sort part , do not use _sort_by_mz
        self._sort_by_mz = sort_by_mz

    def resolve_data(self, **kwargs):
        """Method to swap data."""

    # lazy load properties
    @property
    def mz_list(self):
        """
        Get m/z values. Triggers data loading if necessary.

        Returns:
            np.ndarray: Array of m/z values.
        """
        if self.shared_mz_list is not None:
            self._mz_list = self.shared_mz_list
        elif self._mz_list is None:
            _ = self.intensity  # trigger data load

        return self._mz_list

    @mz_list.setter
    def mz_list(self, value):
        self._mz_list = value

    @property
    def intensity(self):
        """
        Get intensity values. Triggers data loading if necessary.

        Returns:
            np.ndarray: Array of intensity values.
        """
        if self._intensity is None:
            self.resolve_data()

        return self._intensity

    @intensity.setter
    def intensity(self, value):
        self._intensity = value

    @property
    def sort_by_mz(self):
        return self._sort_by_mz

    @sort_by_mz.setter
    def sort_by_mz(self, value: bool = True):
        if value:
            # make sure both mz_list and intensity are not None
            if (self._sort_by_mz is False
                and self.mz_list is not None
                and self.intensity is not None):
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
        Get spatial coordinates.

        Returns:
            PixelCoordinates: The (x, y, z) coordinates object.
        """
        return self.coordinate

    def __str__(self):

        """Return a string representation of the spectrum basic info."""

        if self.mz_list is not None and len(self.mz_list) > 0:
            mz_min, mz_max = min(self.mz_list), max(self.mz_list)
            intensity_max = max(self.intensity) if self.intensity is not None else 0

            return (f"  MS len: {len(self)}\r\n"
                    f"  MS range: {mz_min} - {mz_max}\r\n"
                    f"  MS coord: {self.coordinate}\r\n"
                    f"  max intensity: {intensity_max}\r\n")

        return f"  MS len: 0\r\n  MS coord: {self.coordinate}\r\n"

    def __len__(self):
        """Return number of peaks in the spectrum."""

        if self.mz_list is not None:
            return len(self.mz_list)
        else:
            return 0

    def __eq__(self, other):
        """Check equality based on spatial coordinates."""

        if not isinstance(other, Spectrum):
            return False
        return self.coordinate == other.coordinate

    def __getitem__(self, index):
        """
        Get (m/z, intensity) pair at index.
        """

        if self.mz_list is None or self.intensity is None:
            logger.error("mz_list or intensity is None, can not get item.")
            raise IndexError("mz_list or intensity is None, can not get item.")

        return self.mz_list[index], self.intensity[index]

    def crop_range(
        self,
        mz_range: Sequence[float],
        sort_by_mz: bool = True,
        mode: str = "new",
    ):
        """
        Crop spectrum to a specific m/z range.

        Args:
            mz_range: (min_mz, max_mz) tuple.
            sort_by_mz: Ensure data is sorted before cropping.
            mode: 'new' to return a copy, 'update' to modify in-place.

        Returns:
            Spectrum: The cropped spectrum.
        """

        if mz_range is None or len(mz_range) != 2:
            raise ValueError("mz_range must be a sequence of two values: (min_mz, max_mz).")

        min_mz = float(mz_range[0])
        max_mz = float(mz_range[1])
        if min_mz > max_mz:
            raise ValueError("mz_range is invalid: min_mz must be <= max_mz.")

        if self.mz_list is None or self.intensity is None:
            logger.error("mz_list or intensity is None, cannot crop by mz_range.")
            raise ValueError("mz_list or intensity is None, cannot crop by mz_range.")

        mz_arr = np.asarray(self.mz_list)
        inten_arr = np.asarray(self.intensity)
        if mz_arr.shape[0] != inten_arr.shape[0]:
            raise ValueError("m/z array and intensity array length mismatch.")

        sorted_after_crop = self._sort_by_mz
        if sort_by_mz and not self._sort_by_mz:
            sorted_indices = np.argsort(mz_arr)
            mz_arr = mz_arr[sorted_indices]
            inten_arr = inten_arr[sorted_indices]
            sorted_after_crop = True

        mask = (mz_arr >= min_mz) & (mz_arr <= max_mz)
        mz_c = mz_arr[mask]
        inten_c = inten_arr[mask]

        if mode == "new":
            new_obj = copy.copy(self)
            new_obj.mz_list = mz_c
            new_obj.intensity = inten_c
            new_obj._sort_by_mz = sorted_after_crop
            return new_obj

        if mode == "update":
            self.mz_list = mz_c
            self.intensity = inten_c
            self._sort_by_mz = sorted_after_crop
            return self

        logger.error(f"mode {mode} is not supported. use 'new' or 'update'")
        raise ValueError(f"mode {mode} is not supported. use 'new' or 'update'")


    def is_sorted(self):
        # Check if the mz_list is sorted in ascending order.
        if self.mz_list is not None and self.intensity is not None:
            return np.all(self.mz_list[:-1] <= self.mz_list[1:])

        return False


    def clear_data(self):
        """Clear loaded data from memory (for lazy loading)."""
        self.mz_list = None
        self.intensity = None

    @property
    def x(self):
        """Get X coordinate."""
        if self.coordinate is None:
            logger.error("coordinates is not initialized.")
            raise AttributeError("coordinates is not initialized.")
        return self.coordinate.x

    @property
    def y(self):
        """Get Y coordinate."""
        if self.coordinate is None:
            logger.error("coordinates is not initialized.")
            raise AttributeError("coordinates is not initialized.")
        return self.coordinate.y

    @property
    def z(self):
        """Get Z coordinate."""
        if self.coordinate is None:
            logger.error("coordinates is not initialized.")
            raise AttributeError("coordinates is not initialized.")
        return self.coordinate.z


    def swap_out2disk(self, writer: ImzMLWriter):
        """Method to swap data."""
