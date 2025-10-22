"""
MSI Data Management Module

Provides functions for reading/writing MSI data, memory statistics, and visualization.
Supports .h5/.msi files and batch import from directories, filters by m/z range,
and generates merged or split outputs.
"""

import os
from abc import ABC, abstractmethod
import h5py
import numpy as np
from .msi_module import MSI

class MSIDataManager(ABC):
    """
    Abstract base class for MSI Data Manager.

    This abstract class defines the interface and common functionalities for MSI data management.
    All concrete implementations must inherit from this class and implement the required abstract methods.

    Abstract Methods:
        load_full_data_from_file(filepath): Load complete MSI data from file
    
    Common Methods:
        write2local(filepath): Write MSI data to local file
        inspect_data(): Inspect and analyze MSI data structure

    Attributes:
        _msi (MSI): Internal MSI object for data management
        target_mz_range (tuple, optional): Target m/z range for filtering
        filepath (str, optional): Path to the data file
    """
    def __init__(self,
                 msi: MSI,
                 target_mz_range=None,
                 filepath=None):
        self._msi = msi
        self.target_mz_range = target_mz_range
        self.filepath = filepath
        self.current_image_num = 0

    def get_msi(self) -> MSI:
        """
        Get the MSI object.

        Returns:
            MSI: The MSI object.
        """
        return self._msi

    @abstractmethod
    def load_full_data_from_file(self):
        """
        Load metadata from a file.

        Args:
            filepath (str): Path to the input file.
        """

    def inspect_data(self):
        """
        Inspect the data structure of the MSI object.

        Prints metadata shapes and queue information, including max/min m/z values,
        queue length, and count of non-empty base masks.
        """

        print("MSI meta data:")
        for msi_attr, meta_data in self._msi.metadata.items():
            shape = getattr(meta_data, 'shape', None)
            if shape is not None:
                print(f"    MSI_{msi_attr}: {shape}")
            else:
                print(f"    MSI_{msi_attr}: {meta_data}")

        print("MSI  information:")
        if self._msi.get_queue():
            mz_values = [module.mz for module in self._msi.get_queue()]
            max_mz = max(mz_values)
            min_mz = min(mz_values)
            non_empty_count = sum(1 for module in self._msi.get_queue() if module.base_mask is not None)
            print(f"    MSI max mz: {max_mz}")
            print(f"    MSI min mz: {min_mz}")
            print(f"    MSI len : {len(self._msi.get_queue())}")
            print(f"    base_mask not empty is {non_empty_count}")
            # print(f"MSI queue mz values: {[mz.item() for mz in mz_values]}")
        else:
            print("MSI queue is empty.")

    def _write_meta_data(self, output_path):

        with h5py.File(output_path, 'a') as file_handle:
            self._msi.update_metadata()
            # 先写入不带下划线前缀的属性名（如 meta_version），优先使用 getter 的值
            for private_key, _ in self._msi.metadata.items():
                public_key = private_key[1:]  # 去掉前导下划线：'_meta_xxx' -> 'meta_xxx'
                value = getattr(self._msi, public_key)
                if file_handle.get(public_key) is None:
                    if value is None:
                        print("leak meta data with {value}")
                        continue
                    elif 'num' in public_key or 'version' in public_key:
                        dtype = np.float32
                    else:
                        dtype = h5py.string_dtype(encoding='utf-8') if isinstance(value, str) else None
                    file_handle.create_dataset(public_key, data=value, dtype=dtype)

    def write2local(self, mode="merge", prefix="MSI", output_fold=None, compression_opts=9):
        """
        Write the MSI data to local disk.

        Parameters:
        - mode (str): Writing mode, either "split" or "merge". Default is "merge".
        - prefix (str): Prefix for output file names. Default is "MSI".
        - output_fold (str): Path to the output folder. Default is None.
        - compression_opts (int): Compression level for HDF5 datasets. Default is 9.

        Notes:
        - Writes m/z images and metadata to disk in the specified format.
        - Creates output folder if it does not exist.
        """

        if len(self._msi) == 0:
            print("No MSI images to write. Please rebuild or load the HDF5 file first.")

        self._msi.meta_storage_mode = mode

        output_fold = (
            f'./{prefix}_{self._msi.meta_name}_{self._msi.meta_version}'
            if output_fold is None else output_fold
        )
        os.makedirs(output_fold, exist_ok=True)

        for msi_base in self._msi:
            mz_data = msi_base.mz
            if mode == "split":
                file_name = f"{prefix}_{mz_data:.4f}.msi"
            elif mode == "merge":
                file_name = f"{prefix}_{self._msi.meta_name}_merge_{self._msi.meta_version}.msi"
            else:
                assert False, f"Error: {mode} is not a valid mode. Please use 'split' or 'merge'."

            # update file name
            output_path = os.path.join(output_fold, file_name)

            # filter the mask data
            self._create_datasets(output_path, mz_data, msi_base.msroi, compression_opts,
                                            group_name=f"mz_{mz_data:.4f}")

        if mode == "split" and len(self._msi) > 0:
            output_path = os.path.join(output_fold, f"{prefix}_metadata.msi")
        self._write_meta_data(output_path)

    @staticmethod
    def _create_datasets(output_path, mz_data, msroi, compression_opts, compression='gzip',
                         group_name=None):

        with h5py.File(output_path, 'a') as file_handle:

            if group_name is None:
                group_name = f"mz_{mz_data:.4f}" if group_name is None else 'default'

            if group_name in file_handle:
                group = file_handle[group_name]
            else:
                group = file_handle.create_group(group_name)

            if 'mz' not in group and mz_data is not None:
                group.create_dataset('mz', data=mz_data)

            if 'msroi' not in group and msroi is not None:
                group.create_dataset('msroi',
                                     data=msroi,
                                     compression=compression,
                                     compression_opts=compression_opts)
