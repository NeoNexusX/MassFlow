"""
MSI Data Manager for ZYS Format

This module provides specialized functionality for handling MSI data in ZYS format,
extending the base MSIDataManager with specific methods for reading and processing
.mat files containing MSI data with threshold-based filtering and normalization.
"""
from typing import Optional
import h5py
import numpy as np
from .msi_data_manager import MSIDataManager
from .msi_module import MSIBaseModule


class MSIDataManagerZYS(MSIDataManager):
    """
    MSI Data Manager specialized for ZYS format files.
    
    This class extends MSIDataManager to handle MSI data stored in ZYS format (.mat files).
    It provides functionality for loading, processing, and normalizing MSI data with
    threshold-based filtering and sparsity control.
    
    Attributes:
        threshold (float): Sparsity threshold for filtering m/z channels (default: 0.3)
        h5_data_zys (Optional[h5py.File]): Handle to the opened HDF5/MAT file
    """

    def __init__(self, msi, target_mz_range=None,threshold=0.0,filepath: str = None):
        """
        Initialize MSIDataManagerZYS instance.
        
        Args:
            msi: MSI object to manage data
            target_mz_range (tuple, optional): Target m/z range as (min, max) for filtering
            threshold (float, optional): Sparsity threshold for m/z channel filtering (default: 0.3)
            filepath (str, optional): Path to the input .mat file
        """
        super().__init__( msi, target_mz_range,filepath)
        self.threshold = threshold
        self.h5_data_zys: Optional[h5py.File] = None

    def load_full_data_from_file(self):
        """
        Load ZYS format data from a .mat file.
        
        Opens the .mat file specified in self.filepath and stores the file handle
        for subsequent data access operations.
        
        Raises:
            AssertionError: If the filepath does not end with '.mat'
        """
        assert self.filepath.endswith('.mat'), "Error: filepath is not a .mat file."
        self.h5_data_zys = h5py.File(self.filepath, 'r')

    def get_dataset_from_zys(self, dataset_path: str):
        """
        Retrieve a dataset from the loaded ZYS file.
        
        Args:
            dataset_path (str): Path to the dataset within the HDF5 file structure
            
        Returns:
            h5py.Dataset: The requested dataset object
        """
        return self.h5_data_zys[dataset_path]  # type: ignore

    def get_dataset_numpy_from_zys(self, dataset_path: str):
        """
        Retrieve a dataset from the ZYS file and convert it to a NumPy array.
        
        Args:
            dataset_path (str): Path to the dataset within the HDF5 file structure
            
        Returns:
            numpy.ndarray: The dataset converted to a NumPy array
        """
        dataset = self.get_dataset_from_zys(dataset_path)
        return dataset[()]

    def rebuild_hdf5_file_from_zys(self):
        """
        Rebuild MSI images from a ZYS format HDF5 file with threshold-based filtering.
        
        This method processes the ZYS format data by:
        1. Reading mask, m/z values, and spectral data from the file
        2. Applying sparsity threshold filtering to select valid m/z channels
        3. Normalizing channel data using min-max normalization
        4. Allocating data matrices and populating MSI slices
        5. Handling mask orientation adjustments if needed
        
        The method filters m/z channels based on:
        - Target m/z range (if specified)
        - Sparsity threshold (self.threshold)
        - Data normalization requirements
        
        Note:
            This method modifies the internal MSI object state by setting metadata
            and adding MSI slices to the processing queue.
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
            self._msi.meta.mask = _mask
            self._msi.meta.mz_num = len(selected_channels)
            self._msi.allocate_data_from_meta(dtype=np.float32)

            # Fill valid pixels into the managed data matrix and bind slices
            for i, (mz, channel_data) in enumerate(selected_channels):
                self._msi.data[i, coords[:, 0], coords[:, 1]] = channel_data
                base_mask = np.where(self._msi.data[i] > 0, 1, 0) if self._msi.meta.need_base_mask else None
                self._msi.add_msi_img(
                    MSIBaseModule(
                        mz=mz,
                        msroi=self._msi.data[i],
                        base_mask=base_mask
                    )
                )

            # If msroi shape differs from _mask, fix mask orientation by transposing
            mask_to_set = _mask.T if self._msi.queue[0].msroi.shape != _mask.shape else _mask
            self._msi.meta.mask = mask_to_set
        else:
            # No valid channels; clear m/z count in metadata and keep queue empty
            self._msi.meta.mask = _mask
            self._msi.meta.mz_num = 0
