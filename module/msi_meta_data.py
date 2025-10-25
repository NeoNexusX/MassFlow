# module/msi_meta_data.py

class MSIMetadataBase:
    """
    Abstract base class for MSI data models.

    Manages all common metadata fields, properties, and setters
    that are shared between different data representations.

    All metadata fields starting with _meta_ are automatically synchronized
    to the metadata dictionary via property setters.
    """

    def __init__(self, name, version=1.0, mask=None, mz_num=None,
                 storage_mode='split', need_base_mask: bool = False):

        self._metadata = {}

        # Initialize all metadata fields via properties to trigger auto-sync
        self._meta_name = None
        self._meta_version = None
        self._meta_mask = None
        self._meta_mz_num = None
        self._meta_storage_mode = None
        self._meta_need_base_mask = None

        #set actual values through properties
        self.meta_name = name
        self.meta_version = version
        self.meta_mask = mask
        self.meta_mz_num = mz_num
        self.meta_storage_mode = storage_mode
        self.meta_need_base_mask = need_base_mask

    #Metadata dictionary accessors

    @property
    def metadata(self):
        """Returns the internal dictionary of all metadata."""
        return self._metadata

    def update_metadata(self):
        """
        Manually update the metadata dictionary with all _meta_ attributes.

        Note: With the current implementation, this is automatically called
        by all property setters, so manual calls are usually not needed.
        This method is kept for backward compatibility.
        """
        for msi_attr in dir(self):
            if msi_attr.startswith('_meta_'):
                self._metadata[msi_attr] = getattr(self, msi_attr)

    #Metadata properties with auto-sync

    @property
    def meta_name(self):
        return self._meta_name

    @meta_name.setter
    def meta_name(self, name):
        self._meta_name = name
        self._metadata['_meta_name'] = name

    @property
    def meta_version(self):
        return self._meta_version

    @meta_version.setter
    def meta_version(self, version):
        self._meta_version = version
        self._metadata['_meta_version'] = version

    @property
    def meta_mask(self):
        return self._meta_mask

    @meta_mask.setter
    def meta_mask(self, mask):
        self._meta_mask = mask
        self._metadata['_meta_mask'] = mask

    @property
    def meta_need_base_mask(self):
        return self._meta_need_base_mask

    @meta_need_base_mask.setter
    def meta_need_base_mask(self, need_base_mask):
        self._meta_need_base_mask = bool(need_base_mask)
        self._metadata['_meta_need_base_mask'] = self._meta_need_base_mask

    @property
    def meta_mz_num(self):
        return self._meta_mz_num if self._meta_mz_num is not None else 0

    @meta_mz_num.setter
    def meta_mz_num(self, mz_num):
        if mz_num is None:
            self._meta_mz_num = None
        else:
            self._meta_mz_num = int(mz_num)
        self._metadata['_meta_mz_num'] = self._meta_mz_num

    @property
    def meta_storage_mode(self):
        return self._meta_storage_mode

    @meta_storage_mode.setter
    def meta_storage_mode(self, mode):
        self._meta_storage_mode = mode
        self._metadata['_meta_storage_mode'] = mode