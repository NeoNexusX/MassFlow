import numpy as np
from matplotlib import pyplot as plt

class MSIBaseModule:

    def __init__(self, mz, msroi, base_mask=None):
        #TODO ：mz 要转变成set
        self.mz = mz
        self.msroi = msroi
        self.base_mask = base_mask


class MSI:
    """
    Domain model for MSI data: metadata, slice queue, and optional data matrix.
    Provides API via methods/properties only to avoid external direct field access.
    """

    def __init__(self, name, version=1.0, mask=None, mz_num=None, storage_mode='split', need_base_mask: bool = False):
        # Private fields
        self._meta_name = name
        self._queue = []
        self._meta_mask = mask
        self._meta_storage_mode = storage_mode
        self._meta_version = version
        self._meta_need_base_mask = need_base_mask
        self._meta_mz_num = mz_num
        self._metadata = {}
        self._data = None
        self.update_metadata()

    # ---- Metadata and data accessors ----
    @property
    def metadata(self):
        return self._metadata

    def update_metadata(self):
        for msi_attr in dir(self):
            if msi_attr.startswith('_meta_'):
                self.metadata[msi_attr] = getattr(self, msi_attr)

    @property
    def meta_name(self):
        return self._meta_name

    @meta_name.setter
    def meta_name(self, name):
        self._meta_name = name

    @property
    def meta_version(self):
        return self._meta_version

    @meta_version.setter
    def meta_version(self, version):
        self._meta_version = version

    @property
    def meta_mask(self):
        return self._meta_mask

    @meta_mask.setter
    def meta_mask(self, mask):
        self._meta_mask = mask

    @property
    def meta_need_base_mask(self):
        return self._meta_need_base_mask

    @meta_need_base_mask.setter
    def meta_need_base_mask(self, need_base_mask):
        self._meta_need_base_mask = bool(need_base_mask)

    @property
    def meta_mz_num(self):
        return self._meta_mz_num if self._meta_mz_num is not None else 0

    @meta_mz_num.setter
    def meta_mz_num(self, mz_num):
        if mz_num is None:
            self._meta_mz_num = None
        else:
            self._meta_mz_num = int(mz_num)

    @property
    def meta_storage_mode(self):
        return self._meta_storage_mode

    @meta_storage_mode.setter
    def meta_storage_mode(self, mode):
        self._meta_storage_mode = mode

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data_matrix):
        self._data = data_matrix

    # ---- Queue accessors ----
    def add_msi_slice(self, msi):
        if isinstance(msi, MSIBaseModule):
            self._queue.append(msi)
        else:
            raise ValueError("Only MSIBaseModule instances can be added to the queue.")

    def get_queue(self):
        return self._queue

    def __getitem__(self, index):
        return self._queue[index]

    def __len__(self):
        return len(self._queue)

    def __iter__(self):
        return iter(self._queue)

    # ---- Business methods (query/plot) ----
    def get_msi_by_mz(self, mz_value_min: float, mz_value_max: float = 0, tol=1e-3):
        """Return MSI slices within the given m/z range."""
        mz_value_max = mz_value_min if mz_value_max == 0 else mz_value_max
        buffer = []
        for msi_data in self._queue:
            if (mz_value_min - tol) <= msi_data.mz <= (mz_value_max + tol):
                buffer.append(msi_data)
        return buffer

    def get_image_by_mz(self, mz_value_min: float, mz_value_max: float, tol=1e-3):
        """Return image matrices within the given m/z range."""
        images = []
        for msi_data in self._queue:
            if (mz_value_min - tol) <= msi_data.mz <= (mz_value_max + tol):
                images.append(msi_data.msroi)
        return images

    def plot_msi(self, target_mz_range=None, display_threshold_percent=95,
                 figure_size=(12, 8), cmap='inferno', output_path=None):
        """Plot MSI images within the specified m/z range."""

        target_mz_range = [0,1000] if target_mz_range is None else target_mz_range
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
                    plt.savefig(output_path.format(msi.mz))
                else:
                    plt.show()

    def allocate_data_from_meta(self, dtype=np.float32):
        """Allocate data matrix from current metadata and bind to `data`.

        Requirements:
        - `meta_mask` is set and is a 2D matrix
        - `meta_mz_num` > 0
        """
        assert self.meta_mask is not None, "meta_mask is None"
        assert hasattr(self.meta_mask, 'shape') and len(self.meta_mask.shape) == 2, "meta_mask must be 2D"
        assert self.meta_mz_num > 0, "meta_mz_num must be greater than 0"
        self.data = np.zeros(
            (int(self.meta_mz_num), int(self.meta_mask.shape[0]), int(self.meta_mask.shape[1])),
            dtype=dtype
        )
        return self.data
