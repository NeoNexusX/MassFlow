from typing import Union, Optional, Sequence
import numpy as np
from massflow.module.pixel_coordinates import PixelCoordinates
from massflow.tools.stream_imzml_writer import ImzMLWriter
from massflow.logger import get_logger

logger = get_logger("spectrum")


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
        # need to update
        self.coordinate = PixelCoordinates(coordinate)

        # sort part , do not use _sort_by_mz
        self._sort_by_mz = sort_by_mz

    def _resolve_data(self):
        """Internal hook to trigger lazy loading of data."""

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
        if value is None or isinstance(value, np.ndarray):
            self._mz_list = value

    @property
    def intensity(self):
        """
        Get intensity values. Triggers data loading if necessary.

        Returns:
            np.ndarray: Array of intensity values.
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
        if value:
            # make sure both mz_list and intensity are not None
            if (
                self._sort_by_mz is False
                and self.mz_list is not None
                and self.intensity is not None
            ):
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

        Raises:
            IndexError: If data is not loaded.
        """
        if self.mz_list is None or self.intensity is None:
            logger.error("mz_list or intensity is None, can not get item.")
            raise IndexError("mz_list or intensity is None, can not get item.")

        return self.mz_list[index], self.intensity[index]

    def crop_range(
        self,
        x_range: Sequence[float],
        sort_by_mz: bool = True,
        mode: str = "new",
    ):
        """
        Crop spectrum to a specific m/z range.

        Args:
            x_range: (min_mz, max_mz) tuple.
            sort_by_mz: Ensure data is sorted before cropping.
            mode: 'new' to return a copy, 'update' to modify in-place.

        Returns:
            Spectrum: The cropped spectrum.
        """

        # make sure mz_list and intensity are not None
        if self.mz_list is not None and self.intensity is not None:

            # valid test
            min_mz = np.min(self.mz_list)
            max_mz = np.max(self.mz_list)
            if x_range[0] < min_mz or x_range[1] > max_mz:

                # make sure x_range is valid
                if x_range is not None and len(x_range) == 2:
                    mz_c = None
                    inten_c = None

                    # without sort and dont want to sort
                    if not sort_by_mz and not self._sort_by_mz:
                        mask_xr = np.ones_like(self.mz_list, dtype=bool)
                        mask_xr &= (self.mz_list >= x_range[0]) & (
                            self.mz_list <= x_range[1]
                        )

                        mz_c = self.mz_list[mask_xr]
                        inten_c = self.intensity[mask_xr]

                    # with sort and want to sort  or want sort but no sort
                    elif sort_by_mz:
                        # always set sort to true
                        self.sort_by_mz = True
                        # use searchsorted to find the indices
                        start_index = np.searchsorted(
                            self.mz_list, x_range[0], side="left"
                        )
                        end_index = np.searchsorted(
                            self.mz_list, x_range[1], side="right"
                        )

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
            logger.error(
                "mz_list is None, can not crop by mz. xr must be a tuple of two float numbers."
            )
            raise ValueError(
                "mz_list is None, can not crop by mz. xr must be a tuple of two float numbers."
            )

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
        """Clear loaded data from memory (for lazy loading)."""
        self.mz_list = None
        self.intensity = None

    @property
    def x(self):
        """Get X coordinate."""
        if self.coordinate is None:
            raise AttributeError("coordinates is not initialized.")
        return self.coordinate.x

    @property
    def y(self):
        """Get Y coordinate."""
        if self.coordinate is None:
            raise AttributeError("coordinates is not initialized.")
        return self.coordinate.y

    @property
    def z(self):
        """Get Z coordinate."""
        if self.coordinate is None:
            raise AttributeError("coordinates is not initialized.")
        return self.coordinate.z

    def swap_out2disk(self, writer: ImzMLWriter):
        pass
