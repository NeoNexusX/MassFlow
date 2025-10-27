import numpy as np
from matplotlib import pyplot as plt
from .meta_data import MetaDataBase


class MSIBaseModule:
    """
    Basic MSI data slice containing m/z value, intensity matrix, and optional mask.
    """

    def __init__(self, mz, msroi, base_mask=None):
        self.mz = mz
        self.msroi = msroi
        self.base_mask = base_mask


class MSIMetaData(MetaDataBase):
    """
    Metadata class for MSI data models.
    Inherits all metadata handling from MetaDataBase.
    Adds an additional field for need_base_mask.
    """
    def __init__(self,
                 name,
                 version=1.0,
                 mz_num=None,
                 storage_mode='split',
                 mask=None,
                 need_base_mask: bool = False):

        super().__init__(name, version, mz_num, storage_mode)

        self._mask = None
        self._need_base_mask = None

        # Set actual value through property
        self.mask = mask
        self.need_base_mask = need_base_mask

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, mask):
        if mask is not None:
            self._mask = mask
            self._set('mask', mask)

    @property
    def need_base_mask(self):
        return self._need_base_mask

    @need_base_mask.setter
    def need_base_mask(self, need_base_mask: bool):
        if need_base_mask is not None:
            self._need_base_mask = need_base_mask
            self._set('need_base_mask', need_base_mask)

class MSI:
    """
    Domain model for MSI data (Image Matrix format).
    Adds logic for the slice queue and the 3D data matrix.

    All metadata properties are inherited
    from MSIMetadataBase and automatically synchronized to the metadata dict.
    """

    def __init__(self,
                 name,
                 version=1.0,
                 mask=None,
                 mz_num=None,
                 storage_mode='split',
                 need_base_mask: bool = False):

        # Call parent class __init__ to initialize all metadata
        self.meta = MSIMetaData(name, version, mz_num, storage_mode, mask, need_base_mask)
        # msi.meta[name]
        # Initialize MSI-specific private fields
        self._queue = []
        self._data = None

    # ---- Data matrix property ----
    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data_matrix):
        self._data = data_matrix

    @property
    def queue(self):
        return self._queue

    # ---- Queue accessors ----
    def add_msi_img(self, msi):
        if isinstance(msi, MSIBaseModule):
            self._queue.append(msi)
        else:
            raise ValueError("Only MSIBaseModule instances can be added to the queue.")

    def __getitem__(self, index):
        return self._queue[index]

    def __len__(self):
        return len(self._queue)

    def __iter__(self):
        return iter(self._queue)

    # ---- Business methods ----
    def get_msi_by_mz(self, mz_value_min: float, mz_value_max: float = 0, tol=1e-3):
        """
        Return MSI slices within the given m/z range.
        """
        mz_value_max = mz_value_min if mz_value_max == 0 else mz_value_max
        buffer = []
        for msi_data in self._queue:
            if (mz_value_min - tol) <= msi_data.mz <= (mz_value_max + tol):
                buffer.append(msi_data)
        return buffer

    def get_image_by_mz(self, mz_value_min: float, mz_value_max: float, tol=1e-3):
        """
        Return image matrices within the given m/z range.
        """
        images = []
        for msi_data in self._queue:
            if (mz_value_min - tol) <= msi_data.mz <= (mz_value_max + tol):
                images.append(msi_data.msroi)
        return images

    def plot_msi(self, target_mz_range=None, display_threshold_percent=95,
                 figure_size=(12, 8), cmap='inferno', output_path=None):
        """
        Plot MSI images within the specified m/z range.
        """
        target_mz_range = [0, 1000] if target_mz_range is None else target_mz_range
        assert len(target_mz_range) == 2, "target_mz_range should be a list with two elements"
        assert target_mz_range[0] < target_mz_range[1], "target_mz_range[0] should be less than target_mz_range[1]"

        for msi in self._queue:
            if target_mz_range[0] <= msi.mz <= target_mz_range[1]:
                display_threshold = np.percentile(msi.msroi, display_threshold_percent)
                plt.figure(figsize=figure_size)
                plt.imshow(msi.msroi, aspect='auto', cmap=cmap, vmax=display_threshold)
                plt.colorbar(label='Intensity')
                plt.title(f'MSI Image at m/z {msi.mz:.4f}')
                if output_path:
                    plt.savefig(output_path.format(msi.mz), dpi=300)
                else:
                    plt.show()

    def allocate_data_from_meta(self, dtype=np.float32):

        assert self.meta.mask is not None, "meta_mask cannot be None"
        assert self.meta.mz_num > 0, "meta_mz_num must be greater than 0"

        self.data = np.zeros(
            (int(self.meta.mz_num), int(self.meta.mask.shape[0]), int(self.meta.mask.shape[1])),
            dtype=dtype
        )
        return self.data
