from typing import List, Tuple, Union, Optional
from logger import get_logger
import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser

logger = get_logger("MSBaseModule")

class MSBaseModule:

    def __init__(self,
                mz_list: Optional[np.ndarray],
                intensity: Optional[np.ndarray],
                coordinates: List[int],
                sorted_by_mz_fun: bool = False):

        assert len(coordinates) == 3, "Coordinates must be a list of three integers."

        # 延迟加载
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

    @property
    def intensity(self) -> np.ndarray:
        """
        Get the intensity values of the MSI data.

        Returns:
            np.ndarray: An array of intensity values.
        """
        return self._intensity

    def get_coordinates(self) -> List[int]:
        """
        Get the coordinates of the MSI data.

        Returns:
            List[int]: A list of three integers representing the coordinates (x, y, z).
        """
        return self.coordinates

    def __len__(self):
        return len(self.mz_list)

    def __eq__(self, other):
        if not isinstance(other, MSBaseModule):
            return False
        return self.coordinates == other.coordinates

    def __getitem__(self, index):
        return self.mz_list[index], self.intensity[index]

    def sort_by_mz(self):
        """
        Sort the mz_list and intensity arrays by mz values.
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
        plt.title(f"Mass Spectrum")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.tight_layout()  # Minimize figure padding

        if save_path:
            plt.savefig(save_path,dpi=dpi)
        else:
            plt.show()

class MSImzML(MSBaseModule):
    """
    Lazy spectrum wrapper for ImzML.
    Holds parser + index; loads on-demand.
    """
    def __init__(self,
                parser:ImzMLParser,
                index: int,
                coordinates):

        super().__init__(mz_list=None, intensity=None, coordinates=coordinates)
        self._parser = parser
        self._index = int(index)

    @classmethod
    def creator(cls, parser:ImzMLParser, index: int, coordinates: List[int]):
        return cls(parser, index, coordinates)

    @property
    def mz_list(self):
        if self._mz_list is None:
            mz, intensity = self._parser.getspectrum(self._index)
            self._mz_list = mz
            self._intensity = intensity
        return self._mz_list

    @property
    def intensity(self):
        if self._intensity is None:
            # Ensure mz_list triggers lazy load and sets intensity
            _ = self.mz_list
        return self._intensity

class MS:

    def __init__(self):
        self._queue = []
        self._coordinate_index = {}  # Mapping from coordinates to MSBaseModule
        # TODO: dlq :Initialize metadata with default values
        # self.meta = MSIMetaData(name, version)

    def add_spectrum(self, spectrum: MSBaseModule):
        self._queue.append(spectrum)
        x, y, z = spectrum.x, spectrum.y, spectrum.z

        if z not in self._coordinate_index:
            self._coordinate_index[z] = {}
        if x not in self._coordinate_index[z]:
            self._coordinate_index[z][x] = {}
        self._coordinate_index[z][x][y] = spectrum

    def get_spectrum(self, x: int, y: int, z: int =0 ) -> MSBaseModule:
        return self._coordinate_index[z][x][y]

    def __getitem__(self, key: Union[Tuple[int, int, int], Tuple[int, int], slice]) -> MSBaseModule:
        """
        Support multiple indexing methods:
        - module[x, y, z]     # Complete coordinates
        - module[x, y]        # z defaults to 0
        - module[(x, y, z)]   # Tuple format
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

    def __setitem__(self, key: Union[Tuple[int, int, int], Tuple[int, int]], spectrum: MSBaseModule):
        """
        Support direct assignment:
        - module[x, y, z] = spectrum
        - module[x, y] = spectrum  # z defaults to 0
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
        return len(self._queue)

    def __iter__(self):
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
