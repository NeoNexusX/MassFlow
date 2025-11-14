"""
MSI Data Management Module

Provides functions for reading/writing MSI data, memory statistics, and visualization.
Supports .h5/.msi files and batch import from directories, filters by m/z range,
and generates merged or split outputs.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
import os
from abc import ABC, abstractmethod
import h5py
import numpy as np
from logger import get_logger
from .msi_module import MSI

logger = get_logger("msi_data_manager")


class MSIDataManagerBase(ABC):
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

        meta_data = "MSI meta data:\r\n"
        for attr, value in self._msi.meta.items():
            shape = getattr(value, 'shape', None)
            if shape is not None and len(shape) > 0:
                meta_data+=(f"    meta_{attr}: {shape}\r\n")
            else:
                meta_data+=(f"    meta_{attr}: {value}\r\n")
        logger.info(meta_data)

        if self._msi.queue:
            mz_values = [module.mz for module in self._msi.queue]
            max_mz = max(mz_values)
            min_mz = min(mz_values)
            non_empty_count = sum(1 for module in self._msi.queue if module.base_mask is not None)
            logger.info(f"MSI  information:\r\n"
                        f"    MSI max mz: {max_mz}\r\n"
                        f"    MSI min mz: {min_mz}\r\n"
                        f"    MSI len : {len(self._msi)}\r\n"
                        f"    base_mask not empty is {non_empty_count}\r\n")
            # print(f"MSI queue mz values: {[mz.item() for mz in mz_values]}")
        else:
            logger.info("MSI queue is empty.")

    def _write_meta_data(self, output_path):

        with h5py.File(output_path, 'a') as file_handle:

            for attr, value in self._msi.meta.items():
                ds_name = f"meta_{attr}"
                if value is None:
                    continue
                if isinstance(value, str):
                    dtype = h5py.string_dtype(encoding='utf-8')
                elif ('num' in attr) or ('version' in attr):
                    dtype = np.float32
                else:
                    dtype = None

                # 始终 upsert（覆盖或重建），避免旧值残留
                self._upsert_dataset(file_handle, ds_name, value, dtype=dtype)

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
            logger.info("No MSI images to write. Please rebuild or load the HDF5 file first.")

        self._msi.meta.storage_mode = mode

        meta_name = self._msi.meta.get('name')
        meta_version = self._msi.meta.get('version')

        output_fold = (
            f'./{prefix}_{meta_name}_{meta_version}'
            if output_fold is None else output_fold
        )
        os.makedirs(output_fold, exist_ok=True)

        for msi_base in self._msi:
            mz_data = msi_base.mz
            if mode == "split":
                file_name = f"{prefix}_{mz_data:.4f}.msi"
            elif mode == "merge":
                file_name = f"{prefix}_{meta_name}_merge_{meta_version}.msi"
            else:
                logger.error(f"Error: {mode} is not a valid mode. Please use 'split' or 'merge'.")
                raise ValueError(f"Error: {mode} is not a valid mode. Please use 'split' or 'merge'.")

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
                MSIDataManagerBase._upsert_dataset(group, 'mz', data=mz_data)

            if 'msroi' not in group and msroi is not None:
                MSIDataManagerBase._upsert_dataset(group,
                                               'msroi',
                                                data=msroi,
                                                compression=compression,
                                                compression_opts=compression_opts)

    @staticmethod
    def _upsert_dataset(group, name, data, dtype=None, compression=None, compression_opts=None):
        """
        Update dataset in-place when possible; otherwise delete and recreate.
        """
        if data is None:
            return
        if name in group:
            ds = group[name]
            try:
                ds[...] = data
                return
            except (TypeError, ValueError):
                # dtype/shape not compatible -> drop and recreate
                del group[name]
        kwargs = {}
        if dtype is not None:
            kwargs['dtype'] = dtype
        if compression is not None:
            kwargs['compression'] = compression
            kwargs['compression_opts'] = compression_opts
        group.create_dataset(name, data=data, **kwargs)
