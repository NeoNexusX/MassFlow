"""
Mass Spectrometry Module for MassFlow Framework

This module provides core classes and functionality for handling mass spectrometry (MS) data,
particularly for Mass Spectrometry Imaging (MSI) applications. It includes support for
lazy loading, efficient data management, and visualization capabilities.

Classes:
    MSBaseModule: Base class for mass spectrum data with coordinates
    MSImzML: Specialized class for handling ImzML format with lazy loading
    MS: Collection class for managing multiple mass spectra

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from typing import List, Tuple, Union, Optional
import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser
from logger import get_logger

logger = get_logger("ms_module")

class PixelCoordinates:
    """
    Initialize pixel coordinates with optional zero-based adjustment.

    Args:
        x (int): X coordinate value.
        y (int): Y coordinate value.
        z (int): Z coordinate value.
        coordinate_set_zero (bool, optional): If True, getters will return (value - 1)
            to convert 1-based coordinates to 0-based. Defaults to False.

    Returns:
        None

    Raises:
        TypeError: If any of x, y, or z cannot be converted to int.
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
        self._x = None
        self._y = None
        self._z = None

        self.x = x
        self.y = y
        self.z = z
        self.zero_based = zero_based
    @property
    def x(self) -> int:
        """
        Get the X coordinate.

        Returns:
            int: Adjusted X coordinate (minus 1 if coordinate_set_zero is True).

        Raises:
            None
        """
        return self._x-1 if self.zero_based else self._x

    @x.setter
    def x(self, value: int):
        self._x = value

    @property
    def y(self) -> int:
        """
        Get the Y coordinate.

        Returns:
            int: Adjusted Y coordinate (minus 1 if coordinate_set_zero is True).

        Raises:
            None
        """
        return self._y-1 if self.zero_based else self._y

    @y.setter
    def y(self, value: int):
        self._y = value 

    @property
    def z(self) -> int:
        """
        Get the Z coordinate.

        Returns:
            int: Adjusted Z coordinate (minus 1 if coordinate_set_zero is True).

        Raises:
            None
        """
        return self._z-1 if self.zero_based else self._z

    @z.setter
    def z(self, value: int):
        self._z = value

    @property
    def zero_based(self) -> bool:
        """
        Alias flag indicating whether coordinates are zero-based.

        Returns:
            bool: True if coordinates are zero-based; False otherwise.
        """
        return self._zero_based

    @zero_based.setter
    def zero_based(self, value: bool):
        """
        Set alias flag and synchronize internal state.

        Args:
            value (bool): True for zero-based; False for one-based.

        Returns:
            None
        """
        self._zero_based = bool(value)

    def __eq__(self, other: object) -> bool:
        """
        Compare equality based on adjusted coordinates.

        Args:
            other (PixelCoor): Another PixelCoor instance to compare.

        Returns:
            bool: True if all adjusted coordinates match; False otherwise.

        Raises:
            TypeError: If `other` is not an instance of PixelCoor.
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

    def lefter(self, other: 'PixelCoordinates') -> bool:
        """
        Determine if this pixel is to the left of another pixel.

        Args:
            other (PixelCoor): The other pixel to compare against.

        Returns:
            bool: True if this pixel's X is less than the other's X; otherwise False.

        Raises:
            TypeError: If `other` is not an instance of PixelCoor.
        """
        return self.x < other.x

    def righter(self, other: 'PixelCoordinates') -> bool:
        """
        Determine if this pixel is to the right of another pixel.

        Args:
            other (PixelCoor): The other pixel to compare against.

        Returns:
            bool: True if this pixel's X is greater than the other's X; otherwise False.

        Raises:
            TypeError: If `other` is not an instance of PixelCoor.
        """
        return self.x > other.x

    def upper(self, other: 'PixelCoordinates') -> bool:
        """
        Determine if this pixel is above another pixel.

        Args:
            other (PixelCoor): The other pixel to compare against.

        Returns:
            bool: True if this pixel's Y is greater than the other's Y; otherwise False.

        Raises:
            TypeError: If `other` is not an instance of PixelCoor.
        """
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
        Return the number of dimensions (always 3 for PixelCoordinates).

        Returns:
            int: Always returns 3.

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
        sorted_by_mz_fun (bool): Flag indicating if data is sorted by m/z values
        
    Properties:
        mz_list (np.ndarray): Array of m/z values
        intensity (np.ndarray): Array of intensity values corresponding to m/z values
        
    Note:
        - Coordinates must be a list of exactly three or two integers
        - m/z and intensity arrays must have the same length
        - Lazy loading is supported through None initialization of mz_list and intensity
    """

    def __init__(self,
                mz_list: Optional[np.ndarray],
                intensity: Optional[np.ndarray],
                coordinates: Union[PixelCoordinates, List[int], Tuple[int, int, int]],
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
            AssertionError: If coordinates is not a list of exactly three integers.
        """

        assert len(coordinates) == 3, "Coordinates must be a list of three integers."

        # Lazy loading
        self._mz_list = mz_list
        self._intensity = intensity

        # Normalize coordinates to PixelCoordinates instance
        if isinstance(coordinates, PixelCoordinates):
            self.coordinates = coordinates
        elif isinstance(coordinates, (list, tuple)) and len(coordinates) == 3:
            x, y, z = coordinates
            self.coordinates = PixelCoordinates(int(x), int(y), int(z))
        else:
            raise TypeError("coordinates must be PixelCoordinates or a list/tuple of three ints.")

        if sorted_by_mz_fun:
            self.sorted_by_mz = sorted_by_mz_fun
            self.sort_by_mz()

    # lazy load properties
    @property
    def mz_list(self) -> np.ndarray:
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
    def intensity(self) -> np.ndarray:
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
            
        Example:
            >>> spectrum = MSBaseModule(mz_array, intensity_array, [0, 0, 0])
            >>> print(len(spectrum))  # Output: number of peaks
        """
        return len(self.mz_list)

    def __eq__(self, other):
        """
        Check equality based on coordinates.
        
        Two MSBaseModule instances are considered equal if they have the same coordinates.
        
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
            Tuple[float, float]: Tuple of (m/z, intensity) values at the given index
            
        Raises:
            IndexError: If index is out of range
        """
        return self.mz_list[index], self.intensity[index]

    def sort_by_mz(self):
        """
        Sort the mz_list and intensity arrays by m/z values in ascending order.
        
        This method sorts both arrays simultaneously to maintain correspondence
        between m/z values and their intensities. The operation is performed
        in-place and updates the sorted_by_mz_fun flag.
        
        Note:
            - Only sorts if data is not already sorted (sorted_by_mz_fun is False)
            - Logs a warning if mz_list or intensity is None
            - After sorting, sorted_by_mz_fun flag is set to True
        """
        if self.sorted_by_mz and self.mz_list is not None and self.intensity is not None:

            sorted_indices = np.argsort(self.mz_list)
            self._mz_list = self.mz_list[sorted_indices]
            self._intensity = self.intensity[sorted_indices]

        elif self.mz_list is None or self.intensity is None:
            logger.warning("mz_list or intensity is None, can not sort by mz.")

    def plot(self,
            save_path=None,
            figsize=(20, 5),
            dpi: int = 300,
            color='steelblue',
            plot_mode: str = "line"):
        """
        Plot the mass spectrum.

        Args:
            save_path (str, optional): Path to save the plot. If None, display the plot.
            figsize (tuple, optional): Figure size (width, height) in inches. Defaults to (20, 5).
            dpi (int, optional): Dots per inch for image quality. Defaults to 300.
            color (str, optional): Color of the plot lines. Defaults to 'steelblue'.
            plot_mode (str, optional): Plot mode. "line" for connected line plot, "stem" for default stem plot. Defaults to "line".
        """
        intensity = self.intensity
        mz = self.mz_list

        plt.figure(figsize=figsize)
        mode = (plot_mode or "stem").lower()
        if mode == "line":
            # Connected line plot
            plt.plot(mz, intensity,color=color,linewidth=0.8, alpha=0.8)
        else:
            # Default: stem plot
            markerline, stemlines, baseline = plt.stem(mz, intensity)
            plt.setp(stemlines, linewidth=0.7, color=color, alpha=0.7)
            plt.setp(markerline, markersize=3, color=color, alpha=0.7)
            plt.setp(baseline, linewidth=0.5, color='gray', alpha=0.4)

        # Reduce whitespace by setting tight axis limits
        plt.xlim(mz.min(), mz.max())
        plt.ylim(0, intensity.max() * 1.05)  # 5% margin at top
        plt.title("Mass Spectrum")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.tight_layout()  # Minimize figure padding

        if save_path:
            plt.savefig(save_path,dpi=dpi)
        else:
            plt.show()

    @property
    def x(self) -> int:
        """
        Get the X coordinate from the underlying PixelCoordinates.

        Returns:
            int: The X coordinate.

        Raises:
            AttributeError: If 'coordinates' is not initialized.
        """
        return self.coordinates.x

    @property
    def y(self) -> int:
        """
        Get the Y coordinate from the underlying PixelCoordinates.

        Returns:
            int: The Y coordinate.

        Raises:
            AttributeError: If 'coordinates' is not initialized.
        """
        return self.coordinates.y

    @property
    def z(self) -> int:
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
    
    This class extends MSBaseModule to provide efficient handling of ImzML (Imaging Mass 
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
        - Inherits all visualization and manipulation methods from MSBaseModule
    """
    def __init__(self,
                parser:ImzMLParser,
                index: int,
                coordinates):

        super().__init__(mz_list=None, intensity=None, coordinates=coordinates)
        self._parser = parser
        self._index = int(index)

    @property
    def mz_list(self):
        if self._mz_list is None:
            mz, intensity = self._parser.getspectrum(self._index)
            self._mz_list = mz
            self._intensity = intensity
        return self._mz_list

    @mz_list.setter
    def mz_list(self, value):
        self._mz_list = value

    @property
    def intensity(self):
        if self._intensity is None:
            # Ensure mz_list triggers lazy load and sets intensity
            _ = self.mz_list
        return self._intensity

    @intensity.setter
    def intensity(self, value):
        self._intensity = value

class MS:
    """
    Collection class for managing multiple mass spectra with coordinate-based indexing.
    
    This class serves as a container and manager for multiple MSBaseModule instances,
    providing efficient storage, retrieval, and manipulation of mass spectrometry data
    organized by 3D spatial coordinates. It supports both sequential and coordinate-based
    access patterns.
    
    The class maintains two internal data structures:
    - A queue (_queue) for sequential access and iteration
    - A nested dictionary (_coordinate_index) for fast coordinate-based lookup
    
    Attributes:
        _queue (List[MSBaseModule]): Sequential list of all spectra
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
        self.meta = None
        self._queue = []
        self._coordinate_index = {}  # Mapping from coordinates to MSBaseModule

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
        Add a mass spectrum to the collection with coordinate indexing.
        
        This method adds a spectrum to both the sequential queue and the coordinate
        index for efficient access. The spectrum's coordinates are used to create
        a nested dictionary structure for fast coordinate-based lookup.
        
        Args:
            spectrum (MSBaseModule): Mass spectrum to add to the collection
            
        Note:
            - Automatically extracts coordinates from the spectrum
            - Creates nested dictionary structure if coordinates don't exist
            - Spectrum is added to both queue and coordinate index

        """
        spectrum.coordinates.zero_based = self.meta.coordinates_zero_based
        self._queue.append(spectrum)
        x, y, z = spectrum.x, spectrum.y, spectrum.z

        if z not in self._coordinate_index:
            self._coordinate_index[z] = {}
        if x not in self._coordinate_index[z]:
            self._coordinate_index[z][x] = {}
        self._coordinate_index[z][x][y] = spectrum

    def get_spectrum(self, x: int, y: int, z: int =0 ) -> SpectrumBaseModule:
        """
        Retrieve a mass spectrum by its 3D coordinates.
        
        Args:
            x (int): X coordinate
            y (int): Y coordinate  
            z (int, optional): Z coordinate. Defaults to 0.
            
        Returns:
            MSBaseModule: Mass spectrum at the specified coordinates
            
        Raises:
            KeyError: If no spectrum exists at the specified coordinates
        """
        # Check if the coordinates exist in the index
        if z not in self._coordinate_index or x not in self._coordinate_index[z] or y not in self._coordinate_index[z][x]:
            logger.error(f"No spectrum found at coordinates ({x}, {y}, {z})\r\n"
                            f"min_pixel_x: {self.meta.min_pixel_x}\r\n"
                            f"max_pixel_x: {self.meta.pixel_size_x}\r\n"
                            f"min_pixel_y: {self.meta.min_pixel_y}\r\n"
                            f"max_pixel_y: {self.meta.pixel_size_y}\r\n")
            raise KeyError(f"No spectrum found at coordinates ({x}, {y}, {z})")
        return self._coordinate_index[z][x][y]

    def __getitem__(self,
                    key: Union[Tuple[int, int, int],
                    Tuple[int, int], slice]) -> Union[SpectrumBaseModule, List[SpectrumBaseModule]]:
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
            Union[MSBaseModule, List[MSBaseModule]]: Single spectrum for index/coordinates, 
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
                return self._coordinate_index[0][x][y]  # z defaults to 0
        elif isinstance(key, slice):
            # Support slice access on the internal sequential queue
            return self._queue[key]
        else:
            logger.error("Index must be in tuple format, like [x, y, z] or [x, y]")
            raise TypeError("Index must be in tuple format, like [x, y, z] or [x, y]")

    def __setitem__(self, key: Union[Tuple[int, int, int], Tuple[int, int]], spectrum: SpectrumBaseModule):
        """
        Assign mass spectrum to specific coordinates with automatic indexing.

        Args:
            key (Union[Tuple[int, int], Tuple[int, int, int]]): Target coordinates for spectrum placement.
            spectrum (SpectrumBaseModule): Mass spectrum to assign.

        Returns:
            None

        Raises:
            TypeError: If key is not a tuple.
            IndexError: If tuple length is not 2 or 3.
        """
        if len(key) == 3:
            x, y, z = key
        elif len(key) == 2:
            x, y = key
            z = 0
        else:
            raise IndexError("Coordinates must be 2 or 3 integers")

        # Update spectrum coordinates using PixelCoordinates, respecting meta zero-based setting
        spectrum.coordinates = PixelCoordinates(x, y, z,self.meta.coordinates_zero_based)

        # Add to index
        if z not in self._coordinate_index:
            self._coordinate_index[z] = {}
        if x not in self._coordinate_index[z]:
            self._coordinate_index[z][x] = {}
        if y not in self._coordinate_index[z][x]:
            self._coordinate_index[z][x][y] = spectrum

        # If not in queue, add to queue
        if spectrum not in self._queue:
            self._queue.append(spectrum)

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
            Iterator[MSBaseModule]: Iterator over mass spectra
        """
        return iter(self._queue)

    def plot_spectrum(self,
                x: int = 0,
                y: int = 0,
                z: int = 0,
                save_path=None,
                figsize=(20, 5),
                dpi: int = 300,
                color='steelblue',
                plot_mode: str = "line"):
        """Plot a mass spectrum at the given coordinates.

        Parameters
        - x, y, z: Coordinates of the spectrum to plot (z defaults to 0).
        - save_path: If provided, saves the figure to this path; otherwise, shows it.
        - figsize: Matplotlib figure size tuple.
        - dpi: Figure DPI when saving.
        - plot_mode: 'stem' for stem plot (default), 'line' to connect points with a line.
        """

        spectrum = self.get_spectrum(x, y, z)
        intensity = spectrum.intensity
        mz = spectrum.mz_list

        plt.figure(figsize=figsize)
        mode = (plot_mode or "stem").lower()
        if mode == "line":
            # Connected line plot
            plt.plot(mz, intensity,color=color,linewidth=0.8, alpha=0.8)
        else:
            # Default: stem plot
            markerline, stemlines, baseline = plt.stem(mz, intensity)
            plt.setp(stemlines, linewidth=0.7, color=color, alpha=0.7)
            plt.setp(markerline, markersize=3, color=color, alpha=0.7)
            plt.setp(baseline, linewidth=0.5, color='gray', alpha=0.4)

        # Reduce whitespace by setting tight axis limits
        plt.xlim(mz.min(), mz.max())
        plt.ylim(0, intensity.max() * 1.05)  # 5% margin at top

        plt.title(f"Mass Spectrum at Coordinates (x={x}, y={y}, z={z})")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.tight_layout()  # Minimize figure padding

        if save_path:
            plt.savefig(save_path,dpi=dpi)
        else:
            plt.show()

    def plot_ms_mask(
            self,
            save_path: Optional[str] = None,
            figsize: Tuple[int, int] = (8, 8),
            dpi: int = 300,
            origin: str = 'upper',
            cmap: str = 'Greys'):
        """
        Plot the occupancy mask stored in metadata.

        Parameters
        - save_path (Optional[str]): If provided, save the figure to this path; otherwise, display it.
        - figsize (Tuple[int, int]): Matplotlib figure size (width, height). Defaults to (8, 8).
        - dpi (int): Dots-per-inch when saving. Defaults to 300.
        - origin (str): Image origin for imshow, 'upper' or 'lower'. Defaults to 'upper'.
        - cmap (str): Colormap used to render the mask. Defaults to 'Greys'.

        Returns
        - None: This function produces a plot or saves an image.

        Raises
        - ValueError: If metadata or mask is missing, or mask dimensions are invalid.

        Notes
        - This method directly visualizes `self.meta.mask` without recomputing.
        - Ensure `self.meta.mask` matches `(max_count_of_pixels_y, max_count_of_pixels_x)`.
        """
        # Validate metadata and mask
        if self.meta is None:
            logger.error("MS meta data is required to plot mask.")
            raise ValueError("MS meta data is required to plot mask.")

        mask = getattr(self.meta, 'mask', None)
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
            - This updates `self.meta.coordinates_zero_based` and toggles each spectrum's
              `PixelCoordinates` view via `set_coordinate_zero_based`.
            - The coordinate index is rebuilt to reflect new x/y/z values under the flag.
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

