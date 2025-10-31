# module/msi_meta_data.py
meta_index = {
    "IMS:1000042": "max_count_of_pixels_x", # Image width (pixel count)
    "IMS:1000043": "max_count_of_pixels_y",# Image height (pixel count)
    "IMS:1000044": "max_dimension_x",# Image width (physical size, µm)
    "IMS:1000045": "max_dimension_y",# Image height (physical size, µm)
    "IMS:1000046": "pixel_size_x",# Pixel width (µm)
    "IMS:1000047": "pixel_size_y",# Pixel height (µm)
    "IMS:1000053": "absolute_position_offset_x",# X-axis position offset
    "IMS:1000054": "absolute_position_offset_y",# Y-axis position offset
    "IMS:1000031": "processed",# Whether the data is processed
    "MS:1000031": "instrument_model",# Instrument model
    "MS:1000127": "centroid_spectrum",# Mass spectrum in centroid mode
    "MS:1000128": "profile_spectrum",# Mass spectrum in profile mode
    "MS:1000579": "ms1_spectrum",# MS1 spectrum
    "MS:1000580": "msn_spectrum"# MSn spectrum
}

class MetaDataBase:
    """
    Abstract base class for MSI data models.

    Manages all common metadata fields, properties, and setters
    that are shared between different data representations.

    All metadata fields starting with _meta_ are automatically synchronized
    to the metadata dictionary via property setters.
    """

    def __init__(self, name, version=None, mz_num=None, storage_mode='split'):

        self._meta = {}
        # Initialize all metadata fields via properties to trigger auto-sync
        self._name = None
        self._version = None
        self._mz_num = None
        self._storage_mode = None
        self._meta_index = None

        #set actual values through properties
        self.name = name
        self.version = version
        self.mz_num = mz_num
        self.storage_mode = storage_mode
        self.meta_index = meta_index

    def _set(self, key, value):
        self._meta[key] = value

    #Metadata properties with auto-sync
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if name is not None:
            self._name = name
            self._set('name', name)

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, version):
        if version is not None:
            assert version > 0, "Version must be positive"
            self._version = version
            self._set('version', version)

    @property
    def mz_num(self):
        return self._mz_num

    @mz_num.setter
    def mz_num(self, mz_num):
        if mz_num is not None:
            self._mz_num = int(mz_num)
            self._set('mz_num', self._mz_num)

    @property
    def storage_mode(self):
        return self._storage_mode

    @storage_mode.setter
    def storage_mode(self, mode):
        if mode is not None:
            self._storage_mode = mode
            self._set('storage_mode', self._storage_mode)

    def __getitem__(self, key):
        return self._meta[key]

    def __iter__(self):
        return iter(self._meta)

    def __len__(self):
        return len(self._meta)

    def keys(self):
        return self._meta.keys()

    def items(self):
        return self._meta.items()

    def values(self):
        return self._meta.values()

    def get(self, key, default=None):
        return self._meta.get(key, default)

    def to_dict(self):
        return dict(self._meta)

    @property
    def meta_index(self):
        return self._meta_index

    @meta_index.setter
    def meta_index(self, meta_index):
        if meta_index is not None:
            if not isinstance(meta_index, dict):
                raise TypeError("meta_index must be a dict")
            self._meta_index = meta_index