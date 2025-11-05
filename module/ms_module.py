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

from math import log
from typing import List, Tuple, Union, Optional
import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser
from logger import get_logger

logger = get_logger("ms_module")

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
                coordinates: List[int],
                sorted_by_mz_fun: bool = False):
        """
        Initialize a mass spectrum with m/z values, intensities, and coordinates.
        
        Args:
            mz_list (Optional[np.ndarray]): Array of m/z values. Can be None for lazy loading.
            intensity (Optional[np.ndarray]): Array of intensity values. Can be None for lazy loading.
            coordinates (List[int]): List of three integers [x, y, z] representing spatial coordinates.
            sorted_by_mz_fun (bool, optional): Whether the data is already sorted by m/z. Defaults to False.
            
        Raises:
            AssertionError: If coordinates is not a list of exactly three integers.
        """

        assert len(coordinates) == 3, "Coordinates must be a list of three integers."

        # Lazy loading
        self._mz_list = mz_list
        self._intensity = intensity
        self.coordinates = coordinates
        x, y, z = coordinates
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)
        self.sorted_by_mz_fun = sorted_by_mz_fun

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

    def get_coordinates(self) -> List[int]:
        """
        Get the coordinates of the MSI data.

        Returns:
            List[int]: A list of three integers representing the coordinates (x, y, z).
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
        if self.sorted_by_mz_fun is False and self.mz_list is not None and self.intensity is not None:
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
        """
        Initialize MSImzML with parser, index, and coordinates for lazy loading.
        
        Args:
            parser (ImzMLParser): ImzML parser instance for reading spectrum data
            index (int): Index of the spectrum within the ImzML file
            coordinates (List[int]): 3D coordinates [x, y, z] of the spectrum location
            
        Note:
            The actual m/z and intensity data are not loaded during initialization.
            They will be loaded on first access to the mz_list or intensity properties.
        """

        super().__init__(mz_list=None, intensity=None, coordinates=coordinates)
        self._parser = parser
        self._index = int(index)

    @classmethod
    def creator(cls, parser:ImzMLParser, index: int, coordinates: List[int]):
        """
        Class method factory for creating MSImzML instances.
        
        Args:
            parser (ImzMLParser): ImzML parser instance
            index (int): Spectrum index in the ImzML file
            coordinates (List[int]): 3D coordinates [x, y, z]
            
        Returns:
            MSImzML: New MSImzML instance
        """
        return cls(parser, index, coordinates)

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

    def __getitem__(self, key: Union[Tuple[int, int, int], Tuple[int, int], slice]) -> SpectrumBaseModule:
        """
        Retrieve mass spectrum using flexible indexing methods.
        
        Supports multiple indexing patterns for convenient access to mass spectra:
        - Sequential indexing: ms[index]
        - 2D coordinates: ms[x, y] (z defaults to 0)  
        - 3D coordinates: ms[x, y, z]
        
        Args:
            key (Union[int, Tuple[int, int], Tuple[int, int, int]]): 
                Index or coordinates for spectrum retrieval
                
        Returns:
            MSBaseModule: Mass spectrum at the specified location
            
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
        else:
            raise TypeError("Index must be in tuple format, like [x, y, z] or [x, y]")

    def __setitem__(self, key: Union[Tuple[int, int, int], Tuple[int, int]], spectrum: SpectrumBaseModule):
        """
        Assign mass spectrum to specific coordinates with automatic indexing.
        
        Allows direct assignment of spectra to coordinate positions. The method
        automatically updates the spectrum's coordinates and adds it to both
        the coordinate index and sequential queue if not already present.
        
        Args:
            key (Union[Tuple[int, int], Tuple[int, int, int]]): 
                Target coordinates for spectrum placement
            spectrum (MSBaseModule): Mass spectrum to assign
                
        Raises:
            TypeError: If key is not a tuple
            IndexError: If tuple length is not 2 or 3
            
        Note:
            - Automatically updates spectrum.coordinates to match key
            - Creates coordinate index structure if needed
            - Adds to queue only if spectrum not already present
            - For 2D coordinates, z defaults to 0
        """
        if isinstance(key, tuple):
            if len(key) == 3:
                x, y, z = key
            elif len(key) == 2:
                x, y = key
                z = 0
            else:
                raise IndexError("Coordinates must be 2 or 3 integers")

            # Update spectrum coordinates
            spectrum.coordinates = [x, y, z]
            spectrum.x, spectrum.y, spectrum.z = x, y, z

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
        else:
            raise TypeError("Index must be in tuple format, like [x, y, z] or [x, y]")

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

    def plot_ms(self,
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