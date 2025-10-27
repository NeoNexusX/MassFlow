from typing import List, Tuple, Union, Optional
import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser

class MSBaseModule:

    def __init__(self,
                mz_list: Optional[np.ndarray],
                intensity: Optional[np.ndarray],
                coordinates: List[int]):

        assert len(coordinates) == 3, "Coordinates must be a list of three integers."

        # 延迟加载友好：存私有字段 + 属性访问
        self._mz_list = mz_list
        self._intensity = intensity
        self.coordinates = coordinates
        x, y, z = coordinates
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)

    # 延迟加载兼容的属性接口
    @property
    def mz_list(self) -> np.ndarray:
        return self._mz_list

    @property
    def intensity(self) -> np.ndarray:
        return self._intensity

    def get_coordinates(self) -> List[int]:
        return self.coordinates

    def __len__(self):
        return len(self.mz_list)

    def __eq__(self, other):
        if not isinstance(other, MSBaseModule):
            return False
        return self.coordinates == other.coordinates

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
        if isinstance(key, tuple):
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

    def plot_ms(self, x: int=0, y: int=0, z: int =0):

        spectrum = self.get_spectrum(x, y, z)
        plt.figure(figsize=(10, 6))
        plt.stem(spectrum.mz_list, spectrum.intensity, use_line_collection=True)
        plt.title(f"Mass Spectrum at Coordinates (x={x}, y={y}, z={z})")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.show()