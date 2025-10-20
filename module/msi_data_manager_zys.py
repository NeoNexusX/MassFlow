from typing import Optional
import h5py
import numpy as np
from .msi_data_manager import MSIDataManager
from .msi_module import MSIBaseModule


class MSIDataManagerZYS(MSIDataManager):

    def __init__(self, msi, target_mz_range=None,threshold=0.3,filepath: str = None):
        super().__init__( msi, target_mz_range,filepath)
        self.threshold = threshold
        self.h5_data_zys: Optional[h5py.File] = None

    def load_data_from_zys_mat(self):
        assert self.filepath.endswith('.mat'), "Error: filepath is not a .mat file."
        self.h5_data_zys = h5py.File(self.filepath, 'r')

    def get_dataset_from_zys(self, dataset_path: str):
        return self.h5_data_zys[dataset_path]  # type: ignore

    def get_dataset_numpy_from_zys(self, dataset_path: str):
        dataset = self.get_dataset_from_zys(dataset_path)
        return dataset[()]

    def rebuild_hdf5_file_from_zys(self):
        """Rebuild MSI images from an HDF5 file; each image object contains the m/z value and corresponding image matrix.

        Returns:
            list: List of dictionaries, each containing 'mz' and 'image' keys
        """
        # Read required datasets
        _mask = self.get_dataset_numpy_from_zys('datamsi/mask')
        _mz_values = self.get_dataset_numpy_from_zys('datamsi/mzroi')  # m/z array
        _msroi = self.get_dataset_numpy_from_zys('datamsi/MSroi')  # spectral data

        # Get valid pixel coordinates (assuming mask is a 2D)
        coords = np.argwhere(_mask)
        num_pixels = len(coords)

        # Validate data dimensions
        if _msroi.ndim != 2:
            raise ValueError("MSroi should be a 2D array")

        # Auto-detect data layout: (num_pixels, num_mz) or (num_mz, num_pixels)
        if _msroi.shape[0] == num_pixels:
            _msroi = _msroi.T
        # assume now (num_mz, num_pixels)
        # First select valid channels, then allocate the final matrix to avoid size mismatch
        selected_channels = []  # List[Tuple[mz, channel_data]]

        # Iterate over each m/z channel and collect valid ones to a temporary list
        for i, (mz, sparsity, _) in enumerate(_mz_values):
            if self.target_mz_range is not None and not self.target_mz_range[0] <= mz <= self.target_mz_range[1]:
                continue

            # Extract valid pixel data for the i-th m/z channel
            channel_data = _msroi[i, :]

            # Per-channel normalization
            ch_min = np.min(channel_data)
            ch_max = np.max(channel_data)
            if ch_max - ch_min > 1e-8 and sparsity >= self.threshold:  # threshold
                channel_data = (channel_data - ch_min) / (ch_max - ch_min)
            else:
                continue

            # Do not write to the final matrix yet; record the filtered channel
            selected_channels.append((mz, channel_data))

        # If valid channels exist, allocate the final matrix and add to the queue
        if len(selected_channels) > 0:
            # Set mask first, then allocate the final matrix based on selection count
            self._msi.meta_mask = _mask
            self._msi.meta_mz_num = len(selected_channels)
            self._msi.allocate_data_from_meta(dtype=np.float32)

            # Fill valid pixels into the managed data matrix and bind slices
            for i, (mz, channel_data) in enumerate(selected_channels):
                self._msi.data[i, coords[:, 0], coords[:, 1]] = channel_data
                base_mask = np.where(self._msi.data[i] > 0, 1, 0) if self._msi.meta_need_base_mask else None
                self._msi.add_msi_slice(
                    MSIBaseModule(
                        mz=mz,
                        msroi=self._msi.data[i],
                        base_mask=base_mask
                    )
                )

            # If msroi shape differs from _mask, fix mask orientation by transposing
            mask_to_set = _mask.T if self._msi.get_queue()[0].msroi.shape != _mask.shape else _mask
            self._msi.meta_mask = mask_to_set
        else:
            # No valid channels; clear m/z count in metadata and keep queue empty
            self._msi.meta_mask = _mask
            self._msi.meta_mz_num = 0
