"""
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
import numpy as np
from matplotlib import pyplot as plt
from logger import get_logger
from .meta_data import MSIMetaData

logger = get_logger("msi_module")


class MSIBaseModule:
    """
    Basic MSI data slice containing m/z value, intensity matrix, and optional mask.
    """

    def __init__(self, mz, msroi, base_mask=None):
        self.mz = mz
        self.msroi = msroi
        self.base_mask = base_mask

# MSIMetaData has been relocated to module/meta_data.py
class MSI:
    """
    Domain model for MSI data (Image Matrix format).
    Adds logic for the slice queue and the 3D data matrix.
    All metadata properties are inherited from MSIMetaData
    and automatically synchronized to the metadata dict.
    """

    def __init__(self,
                 meta:MSIMetaData=None):

        # Call parent class __init__ to initialize all metadata
        # self.meta = MSIMetaData(name, version, mz_num, storage_mode, mask, need_base_mask)
        self.meta = meta if meta is not None else MSIMetaData()
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
        """Return the MSI slice queue."""
        return self._queue

    # ---- Queue accessors ----
    def add_msi_img(self, msi):
        """Add a MSIBaseModule instance to the queue."""
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

        Parameters:
            target_mz_range (list|tuple|None): Two-element range [min_mz, max_mz].
                                               If None, defaults to [0, 1000].
            display_threshold_percent (int|float): Percentile for display clipping.
            figure_size (tuple): Matplotlib figure size (width, height).
            cmap (str): Colormap name used for imshow.
            output_path (str|None): Template path to save each image, should
                                    contain a placeholder for mz, e.g. "out_{:.4f}.png".

        Returns:
            None

        Raises:
            ValueError: If `target_mz_range` is not a 2-element sequence or
                        min value is not less than max value.

        Notes:
            This method validates `target_mz_range` with explicit checks,
            logging errors via `logger.error` before raising exceptions.
        """
        target_mz_range = [0, 1000] if target_mz_range is None else target_mz_range

        if not len(target_mz_range) == 2:
            logger.error("target_mz_range should be a list with two elements")
            raise ValueError("target_mz_range should be a list with two elements")

        if not (target_mz_range[0] < target_mz_range[1]):
            logger.error("target_mz_range[0] should be less than target_mz_range[1]")
            raise ValueError("target_mz_range[0] should be less than target_mz_range[1]")

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
        """
        Allocate the 3D data matrix according to metadata.

        Parameters:
            dtype (numpy.dtype): The dtype used for the allocated data matrix.

        Returns:
            numpy.ndarray: The allocated 3D data matrix with shape
                           (mz_num, mask_height, mask_width).

        Raises:
            ValueError: If `meta.mask` is None or `meta.mz_num` is not greater than 0.

        Notes:
            This method logs detailed error messages when preconditions are not met,
            instead of relying on Python `assert`, so it works even under
            optimization mode (`python -O`).
        """

        if self.meta is None or self.meta.mask is None or self.meta.mz_num is None:
            logger.error(" should load meta first ,meta or meta element cannot be None")
            raise ValueError("should load meta first ,meta or meta element cannot be None")

        if not self.meta.mz_num > 0:
            logger.error("meta_mz_num must be greater than 0")
            raise ValueError("meta_mz_num must be greater than 0")

        self.data = np.zeros(
            (int(self.meta.mz_num), int(self.meta.mask.shape[0]), int(self.meta.mask.shape[1])),
            dtype=dtype
        )
        return self.data