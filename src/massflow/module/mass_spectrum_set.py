"""
Mass Spectrometry Module for MassFlow Framework

This module provides the collection class for managing multiple mass spectra.

Classes:
    MassSpectrumSet: Collection class for managing multiple mass spectra

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from typing import List, Tuple,Literal,overload,Union, Optional
import matplotlib.pyplot as plt
from massflow.logger import get_logger
from massflow.module.ms_meta_data import MetaDataBase,ImzMlMetaData
from massflow.module.spectrum import Spectrum


logger = get_logger("ms_module")

class MassSpectrumSet:
    """
    Collection class for managing multiple mass spectra with coordinate-based indexing.

    This class serves as a container and manager for multiple `Spectrum` instances,
    providing efficient storage, retrieval, and manipulation of mass spectrometry data
    organized by 3D spatial coordinates. It supports both sequential and coordinate-based
    access patterns.

    The class maintains two internal data structures:
    - A queue (_queue) for sequential access and iteration
    - A nested dictionary (_coordinate_index) for fast coordinate-based lookup

    Attributes:
        _queue (List[Spectrum]): Sequential list of all spectra
        _coordinate_index (dict): Nested dictionary mapping coordinates to spectra.
                                 Structure: {z: {x: {y: Spectrum}}}

    Indexing Methods:
        - ms[index]: Access by sequential index
        - ms[x, y]: Access by coordinates (z is resolved automatically; see Notes)
        - ms[x, y, z]: Access by full 3D coordinates
        - ms[x, y, z] = spectrum: Direct assignment

        Notes:
                - Coordinates are automatically managed and indexed.
                - Supports both 2D (x, y) and 3D (x, y, z) coordinate systems.
                - Efficient lookup performance through coordinate indexing.
                - No explicit thread-safety guarantees are provided.
                - For 2D access (`ms[x, y]`), `__getitem__` selects `z=0` when present; otherwise it selects
                    the first available z-key in `self._coordinate_index` (iteration order).
                - For 2D assignment (`ms[x, y] = spectrum`), `__setitem__` assigns into the first available
                    z-key in `self._coordinate_index` and will raise if the index is empty.
    """

    def __init__(self):
        """
        Initialize an empty MS collection.

        Creates empty internal data structures for storing and indexing mass spectra.
        No parameters are required for initialization.
        """
        self.meta :MetaDataBase = ImzMlMetaData(filepath=None)
        self._queue = []
        self._coordinate_index = {}  # Mapping from coordinates to MSBaseModule
        self._shared_mz_list = None


    @property
    def shared_mz_list(self):
        """
        Get or set the shared/common m/z list across spectra.

        In continuous-data mode, multiple spectra can share the same m/z axis. This property stores
        that shared m/z array/list and is typically populated by a data manager (e.g. ImzML loader).

        Returns:
                        Any: The stored shared m/z list/array. Defaults to `None`.

                Notes:
                        - The setter only updates the stored value when `self.meta` exists and
                            `self.meta.continuous is True`.
                        - If metadata is missing, or metadata indicates non-continuous data, the setter keeps the
                            previous value and logs a warning.
        """
        return self._shared_mz_list
   
    @shared_mz_list.setter
    def shared_mz_list(self, value):
        if self.meta:
            if self.meta.continuous is True:
                self._shared_mz_list = value
            else:
                logger.warning("Meta data indicates non-continuous data; shared_mz_list not updated.")

    @property
    def coordinate_index(self):
        """
        Get the nested dictionary of coordinate indices.

        Returns:
            dict: Nested dictionary mapping coordinates to spectra.
                  Structure: {z: {x: {y: Spectrum}}}
        """
        return self._coordinate_index

    def add_spectrum(self, spectrum: Spectrum):
        """
        Add or update a mass spectrum with coordinate indexing.

        This method ensures the spectrum is accessible via both the sequential queue
        and the coordinate index. If a spectrum already exists at the given coordinates,
        it will be updated (replaced) in both the coordinate index and the queue to
        avoid duplicates.

        Args:
            spectrum (Spectrum): Mass spectrum to add or update.


        Notes:
            - Automatically extracts coordinates from the spectrum.
            - Creates nested dictionary structure if coordinates don't exist.
            - Updates the existing spectrum in place when coordinates already exist.
        """
        self.update_spectrum_with_coord(spectrum = spectrum)

    def get_spectrum(self, x: int, y: int, z: int = 0) -> Spectrum:
        """
        Retrieve a mass spectrum by its 3D coordinates.

        Args:
            x (int): X coordinate
            y (int): Y coordinate
            z (int, optional): Z coordinate. Defaults to 0.

        Returns:
            Spectrum: Mass spectrum at the specified coordinates

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
                f"Plot mask to see the available locations; use ms.coordinate_index to list all coordinates."
            )
            raise KeyError(f"No spectrum found at coordinates ({x}, {y}, {z})")
        return self._coordinate_index[z][x][y]

    def update_spectrum_with_coord(self,spectrum: Spectrum,x=-1,y=-1,z=-1):
        """
        Add or replace a spectrum in the coordinate index and queue.

        Args:
            spectrum (Spectrum): Spectrum to insert or update.
            x (int): Optional override for X coordinate; uses `spectrum.x` when -1.
            y (int): Optional override for Y coordinate; uses `spectrum.y` when -1.
            z (int): Optional override for Z coordinate; uses `spectrum.z` when -1.

        Returns:
            None

        Raises:
            ValueError: If updating an existing spectrum fails to synchronize in the queue.
        """
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

    def update_spectrum_with_index(self,index:int,spectrum: Spectrum):
        """
        Replace a spectrum at a sequential index and synchronize the coordinate index.

        Args:
            index (int): Position in the internal queue to replace.
            spectrum (Spectrum): Spectrum to assign at the given index.

        Returns:
            None

        Raises:
            IndexError: If `index` is negative or out of range.
            ValueError: If the provided spectrum does not match the existing entry at `index`.

        Notes:
                        - This method checks equality against the existing queue entry
                            (`self._queue[index] == spectrum`) before replacing. This is stricter than coordinate
                            matching and depends on `Spectrum.__eq__` semantics.
                        - If the index is out of range, it falls back to `update_spectrum_with_coord`.
        """
        if index < 0 or isinstance(index,int) is False :
            logger.error(f"Index {index} out of range.")
            raise IndexError(f"Index {index} out of range.")

        elif index < len(self._queue):
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

    @overload
    def __getitem__(self, key: int) -> Spectrum: ...
    @overload
    def __getitem__(self, key: Tuple[int, int]) -> Spectrum: ...
    @overload
    def __getitem__(self, key: Tuple[int, int, int]) -> Spectrum: ...
    @overload
    def __getitem__(self, key: slice) -> List[Spectrum]: ...

    def __getitem__(self,
                    key: Union[int,Tuple[int, int, int], Tuple[int, int], slice]
    ) -> Union[Spectrum, List[Spectrum]]:
        """
        Retrieve mass spectrum using flexible indexing methods.

        Supports multiple indexing patterns for convenient access to mass spectra:
        - Sequential indexing: ms[index]
        - 2D coordinates: ms[x, y] (z is resolved automatically)
        - 3D coordinates: ms[x, y, z]
        - Slice: ms[a:b:c] returns a list of spectra from the internal queue

        Args:
            key (Union[int, Tuple[int, int], Tuple[int, int, int], slice]):
                Index, coordinates, or slice for spectrum retrieval

        Returns:
            Union[Spectrum, List[Spectrum]]: Single spectrum for index/coordinates,
            or a list of spectra for slice access

        Raises:
            TypeError: If key format is not supported
            KeyError: If coordinates don't exist in the collection
            IndexError: If sequential index is out of range

        Notes:
            - For a 2-tuple `(x, y)`, this method uses `z=0` if `0` exists in the index;
                            otherwise it uses the first available z-key in `self._coordinate_index`.
                        - If `self._coordinate_index` is empty and a 2D key is used, `next(iter(...))` raises
                            `StopIteration`.
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
        spectrum: Spectrum,
    ):
        """
        Assign or replace a spectrum by index or coordinates with bounds checking and index sync.

        Args:
            key (Union[int, Tuple[int, int], Tuple[int, int, int]]): Sequential index or target coordinates.
            spectrum (Spectrum): Spectrum to assign.

        Returns:
            None

        Raises:
            IndexError: If coordinate tuple length is not 2 or 3, or sequential index is out of range.
            ValueError: If assigning by index but the spectrum's coordinates do not match the existing entry.

        Notes:
            - For a 2-tuple `(x, y)`, this method assigns into the first available z-key in
              `self._coordinate_index` (and will raise if the index is empty).
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
            Iterator[Spectrum]: Iterator over mass spectra
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
            origin (str): Image origin for `imshow`, either 'upper' or 'lower'. Defaults to 'lower'.
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


