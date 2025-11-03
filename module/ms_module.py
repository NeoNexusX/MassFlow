from typing import List, Tuple, Union, Optional
import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser
from logger import get_logger
from meta_data import MetaDataBase
import os
from pyimzml.metadata import ParamGroup
logger = get_logger("ms_module")

class MSBaseModule:

    def __init__(self,
                mz_list: Optional[np.ndarray],
                intensity: Optional[np.ndarray],
                coordinates: List[int]):

        assert len(coordinates) == 3, "Coordinates must be a list of three integers."

        # Lazy loading
        self._mz_list = mz_list
        self._intensity = intensity
        self.coordinates = coordinates
        x, y, z = coordinates
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)

    # lazy load properties
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

    def plot(self,
            save_path=None,
            figsize=(20, 5),
            dpi: int = 300,
            color='steelblue',
            plot_mode: str = "line"):

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

class MetaDataImzMl(MetaDataBase):
    """ImzML metadata wrapper that loads and caches frequently used fields."""
    def __init__(self,
                 name="MSI",
                 version=1.0,
                 mz_num=None,
                 storage_mode='split',
                 parser: ImzMLParser = None,
                 filepath: str = None
                 ):
        """Initialize the metadata object with either a parser or a file path."""
        super().__init__(name, version, mz_num, storage_mode)

        self._filepath = None
        self._parser = None
        self._spectrum_count_num = None
        self._max_count_of_pixels_x = None
        self._max_count_of_pixels_y = None
        self._max_dimension_x = None
        self._max_dimension_y = None
        self._pixel_size_x = None
        self._pixel_size_y = None
        self._absolute_position_offset_x = None
        self._absolute_position_offset_y = None
        self._processed = None
        self._instrument_model = None
        self._centroid_spectrum = None
        self._profile_spectrum = None
        self._ms1_spectrum = None
        self._msn_spectrum = None

        # Set actual value through property
        if parser is not None:
            self.parser = parser
        elif filepath is not None:
            self.filepath = filepath
        else:
            raise ValueError("Either parser or filepath must be provided")

        if self.parser is not None:
            self.spectrum_count_num = len(self.parser.coordinates)
            self.extract_metadata()  # Use pyimzml.ImzMLParser to retrieve metadata

    @property
    def filepath(self):
        """Return the associated imzML file path."""
        return self._filepath

    @filepath.setter
    def filepath(self, filepath: str):
        """Set the imzML file path and initialize the parser if needed."""
        if not filepath or not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        self._filepath = filepath
        if self.parser is None:
            self.parser = ImzMLParser(filepath)

    @property
    def parser(self):
        """Return the bound pyimzML parser instance."""
        return self._parser

    @parser.setter
    def parser(self, parser: ImzMLParser):
        """Set the pyimzML parser instance."""
        if parser is not None and not isinstance(parser, ImzMLParser):
            raise TypeError("parser must be an instance of pyimzML.ImzMLParser or None")
        self._parser = parser

    @property
    def spectrum_count_num(self):
        """Return the number of spectra."""
        return self._spectrum_count_num

    @spectrum_count_num.setter
    def spectrum_count_num(self, spectrum_count_num: int):
        """Persist the number of spectra and sync it to the base storage."""
        self._spectrum_count_num = spectrum_count_num
        self._set('spectrum_count_num', spectrum_count_num)

    @property
    def max_count_of_pixels_x(self):
        """Return the pixel count along the X axis."""
        return self._max_count_of_pixels_x

    @max_count_of_pixels_x.setter
    def max_count_of_pixels_x(self, max_count_of_pixels_x):
        """Set the pixel count along the X axis."""
        if max_count_of_pixels_x is not None:
            self._max_count_of_pixels_x = max_count_of_pixels_x
            self._set('max_count_of_pixels_x', max_count_of_pixels_x)

    @property
    def max_count_of_pixels_y(self):
        """Return the pixel count along the Y axis."""
        return self._max_count_of_pixels_y

    @max_count_of_pixels_y.setter
    def max_count_of_pixels_y(self, max_count_of_pixels_y):
        """Set the pixel count along the Y axis."""
        if max_count_of_pixels_y is not None:
            self._max_count_of_pixels_y = max_count_of_pixels_y
            self._set('max_count_of_pixels_y', max_count_of_pixels_y)

    @property
    def max_dimension_x(self):
        """Return the physical dimension along the X axis."""
        return self._max_dimension_x

    @max_dimension_x.setter
    def max_dimension_x(self, max_dimension_x):
        """Set the physical dimension along the X axis."""
        self._max_dimension_x = max_dimension_x
        self._set('max_dimension_x', max_dimension_x)

    @property
    def max_dimension_y(self):
        """Return the physical dimension along the Y axis."""
        return self._max_dimension_y

    @max_dimension_y.setter
    def max_dimension_y(self, max_dimension_y):
        """Set the physical dimension along the Y axis."""
        self._max_dimension_y = max_dimension_y
        self._set('max_dimension_y', max_dimension_y)

    @property
    def pixel_size_x(self):
        """Return the pixel size along the X axis."""
        return self._pixel_size_x

    @pixel_size_x.setter
    def pixel_size_x(self, pixel_size_x):
        """Set the pixel size along the X axis."""
        self._pixel_size_x = pixel_size_x
        self._set('pixel_size_x', pixel_size_x)

    @property
    def pixel_size_y(self):
        """Return the pixel size along the Y axis."""
        return self._pixel_size_y

    @pixel_size_y.setter
    def pixel_size_y(self, pixel_size_y):
        """Set the pixel size along the Y axis."""
        self._pixel_size_y = pixel_size_y
        self._set('pixel_size_y', pixel_size_y)

    @property
    def absolute_position_offset_x(self):
        """Return the absolute position offset on the X axis."""
        return self._absolute_position_offset_x

    @absolute_position_offset_x.setter
    def absolute_position_offset_x(self, absolute_position_offset_x):
        """Set the absolute position offset on the X axis."""
        self._absolute_position_offset_x = absolute_position_offset_x
        self._set('absolute_position_offset_x', absolute_position_offset_x)

    @property
    def absolute_position_offset_y(self):
        """Return the absolute position offset on the Y axis."""
        return self._absolute_position_offset_y

    @absolute_position_offset_y.setter
    def absolute_position_offset_y(self, absolute_position_offset_y):
        """Set the absolute position offset on the Y axis."""
        self._absolute_position_offset_y = absolute_position_offset_y
        self._set('absolute_position_offset_y', absolute_position_offset_y)

    @property
    def processed(self):
        """Return whether the data has been processed."""
        return self._processed
    @processed.setter
    def processed(self, processed):
        """Set whether the data has been processed."""
        self._processed = processed
        self._set('processed', processed)

    @property
    def instrument_model(self):
        """Return the mass spectrometer model."""
        return self._instrument_model

    @instrument_model.setter
    def instrument_model(self, instrument_model):
        """Set the mass spectrometer model."""
        self._instrument_model = instrument_model
        self._set('instrument_model', instrument_model)

    @property
    def centroid_spectrum(self):
        """Return whether centroid spectra are present."""
        return self._centroid_spectrum

    @centroid_spectrum.setter
    def centroid_spectrum(self, centroid_spectrum):
        """Set whether centroid spectra are present."""
        self._centroid_spectrum = centroid_spectrum
        self._set('centroid_spectrum', centroid_spectrum)

    @property
    def profile_spectrum(self):
        """Return whether profile spectra are present."""
        return self._profile_spectrum

    @profile_spectrum.setter
    def profile_spectrum(self, profile_spectrum):
        """Set whether profile spectra are present."""
        self._profile_spectrum = profile_spectrum
        self._set('profile_spectrum', profile_spectrum)

    @property
    def ms1_spectrum(self):
        """Return whether MS1 spectra are present."""
        return self._ms1_spectrum

    @ms1_spectrum.setter
    def ms1_spectrum(self, ms1_spectrum):
        """Set whether MS1 spectra are present."""
        self._ms1_spectrum = ms1_spectrum
        self._set('ms1_spectrum', ms1_spectrum)

    @property
    def msn_spectrum(self):
        """Return whether MSn spectra are present."""
        return self._msn_spectrum

    @msn_spectrum.setter
    def msn_spectrum(self, msn_spectrum):
        """Set whether MSn spectra are present."""
        self._msn_spectrum = msn_spectrum
        self._set('msn_spectrum', msn_spectrum)

    # Approach 1: load metadata through pyimzML
    def extract_metadata(self):
        """Iterate _meta_index and populate matching attributes from the parser."""

        logger.info("Extracting metadata...")

        if self._parser is None:
            logger.error("Parser is not initialized. Please set parser or filepath first.")

        for accession_id, prop_name in self._meta_index.items():
            param_value = self.find_param_by_accession_id(accession_id)
            if param_value is not None:
                setattr(self, prop_name, param_value)


    def find_param_by_accession_id(self, accession_id: str): # Use pyimzML.ImzMLParser to fetch metadata
        """Search the predefined metadata areas for the given accession identifier."""

        search_areas = [
            self.parser.metadata.file_description,  # File description (data type, creation time, etc.)
            self.parser.metadata.scan_settings,  # Scan settings (scan mode, m/z range, etc.)
            self.parser.metadata.instrument_configurations,  # Instrument configuration (model, ion source, etc.)
            self.parser.metadata.samples,  # Sample information (sample name, preparation, etc.)
            self.parser.metadata.softwares,  # Software information (parser, version, etc.)
            self.parser.metadata.data_processings,  # Data processing (peak picking, normalization, etc.)
            self.parser.metadata.referenceable_param_groups,  # Referenceable parameter groups (shared metadata)
        ]

        for area in search_areas:
            if area is None:
                continue

            result = self._search_in_area(area, accession_id)
            if result is not None:
                return result

        return None

    def _search_in_area(self, area, accession_id):
        """Search a single parameter area for the given accession identifier."""
        if isinstance(area, ParamGroup):
            if accession_id in area:
                return area[accession_id]

        elif isinstance(area, dict):
            for param_group in area.values():
                if accession_id in param_group:
                    return param_group[accession_id]

        return None
