# module/msi_meta_data.py

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

        #set actual values through properties
        self.name = name
        self.version = version
        self.mz_num = mz_num
        self.storage_mode = storage_mode

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