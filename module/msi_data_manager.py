"""
MSI Data Management Module

Provides functions for reading/writing MSI data, memory statistics, and visualization.
Supports .h5/.msi files and batch import from directories, filters by m/z range,
and generates merged or split outputs.
"""
import sys
import os
import glob
import h5py
from pympler import tracker
import numpy as np

from .msi_module import MSI, MSIBaseModule

class MSIDataManager:
    """
    MSI Data Manager

    Handles batch read/write of .h5/.msi files, filters data by m/z range,
    and provides memory statistics, visualization, and split/merge output.
    """
    def __init__(self,
                 msi: MSI,
                 target_mz_range=None,
                 filepath=None):

        self.filepath = filepath
        self.current_num = 0
        self.me_tr = None
        # self.memory_tracker()
        self._msi = msi
        self.target_mz_range = target_mz_range

    def get_msi(self) -> MSI:
        return self._msi

    #file action
    def load_full_data_from_file(self):
        """
        Load MSI data completely from files (metadata + all m/z images).

        Steps:
        1. Call load_data_helper to read all metadata.
        2. Preallocate a 3D data matrix using meta_mz_num and mask dimensions.
        3. Call load_data_helper again to read all image data and fill the matrix.

        Notes:
        - Requires meta_mz_num > 0, otherwise an assertion error is raised.
        - If target_mz_range is set, only load m/z images within that range.
        """

        # Read metadata
        self.load_data_helper(fn_name='meta')
        assert self.get_msi().meta_mz_num > 0, "meta_mz_num must be greater than 0"
        # Initialize data storage matrix via MSI API
        self.get_msi().allocate_data_from_meta(dtype=np.float32)
        self.load_data_helper(fn_name='data')
        # self.me_tr.store_summary('after_load_data')

    def load_data_helper(self, fn_name = 'meta'):
        # Find all .h5 files in the specified directory
        fn = self.load_meta_from_file if fn_name == 'meta' else self.load_data_from_file

        if (self.filepath.endswith('.h5') or self.filepath.endswith('.msi') )and os.path.isfile(self.filepath):
            fn(self.filepath)
        elif os.path.isdir(self.filepath):
            msi_files = glob.glob(os.path.join(self.filepath, "*.msi"))
            for msi_file in msi_files:
                fn(msi_file)
        else:
            assert False, f"Error: {self.filepath} is not a valid .h5 file or .msi file or directory."

        if fn_name != 'meta':
            # Check if current_num exceeds meta_mz_num after loading data
            assert self.current_num <= self._msi.meta_mz_num, (
                f"current_num {self.current_num} != meta_mz_num "
                f"{self._msi.meta_mz_num}"
            )
            #update meta_mz_num if current_num is smaller
            if self.target_mz_range is not None and self.current_num <= self._msi.meta_mz_num:
                self._msi.meta_mz_num = self.current_num

    def load_meta_from_file(self, file):
        with h5py.File(file, 'r') as h5_file:
            for key, group in h5_file.items():
                if key.startswith('meta_') and self._msi.metadata.get(key) is None:
                    # Read raw dataset value
                    dataset_value = group[()]
                    # Decode bytes to string if necessary
                    if isinstance(dataset_value, bytes):
                        dataset_value = dataset_value.decode('utf-8')
                    setattr(self._msi, key, dataset_value)

        self._msi.update_metadata()

    def load_data_from_file(self, file):
        with h5py.File(file, 'r') as h5_file:
            for key, group in h5_file.items():
                if isinstance(group, h5py.Group) and not key.startswith('meta_'):
                    mz= group['mz'][()]
                    if self.target_mz_range is not None and not self.target_mz_range[0] <= mz <= self.target_mz_range[1]:
                        continue
                    # Compute base_mask based on metadata if needed
                    msi_image = group['msroi'][()]
                    base_mask = np.where(msi_image > 0, 1, 0) if self._msi.meta_need_base_mask else None

                    self._msi.data[self.current_num, :, :] = msi_image
                    self._msi.add_msi_slice(
                        MSIBaseModule(mz=mz,
                                      msroi=(self._msi.data[self.current_num] if self._msi.data is not None else msi_image),
                                      base_mask=base_mask
                                      )
                    )
                    print(f'loading {key}')
                    self.current_num += 1
            print(f"finish loading {file}")

    def write2local(self, mode="merge", prefix="MSI", output_fold=None, compression_opts=9):
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
            MSIDataManager._create_datasets(output_path, mz_data, msi_base.msroi, compression_opts,
                                            group_name=f"mz_{mz_data:.4f}")

        if mode == "split" and len(self._msi) > 0:
            output_path = os.path.join(output_fold, f"{prefix}_metadata.msi")
        self._write_meta_data(output_path)

    def _inspect_hdf5_structure(self, group, indent=0, max_depth=2):

        if indent > max_depth:
            return

        for key in group.keys():
            # Skip the '#refs#' group (do not traverse references)
            item = group[key]
            # Compute indentation for hierarchical display
            prefix = '    ' * indent

            # If it is a subgroup, recurse
            if isinstance(item, h5py.Group):
                print(f"{prefix}Group: {key}")
                self._inspect_hdf5_structure(item, indent + 1, max_depth)

            # If it is a dataset, print name and shape (skip reference type)
            elif isinstance(item, h5py.Dataset):
                # Skip reference-type datasets (dtype.kind == 'O')
                print(f"{prefix}Dataset: {key}  Shape: {item.shape}  type: {item.dtype}")

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

    def inspect_data(self):

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

    def calculate_memory_usage(self):

        print("=== memory usage ===")

        # Metadata section
        print("\n--- Metadata part ---")
        metadata_total_size = 0

        # Update metadata to ensure latest values
        self._msi.update_metadata()

        for key, meta_item in self._msi.metadata.items():
            item_size = sys.getsizeof(meta_item)
            metadata_total_size += item_size
            print(f"Metadata[{key}]: {item_size} bytes")

        print(f"Metadata fullsize: {metadata_total_size} bytes ({metadata_total_size / 1024:.2f} KB)")

        # Data matrix section
        data_matrix_size = 0
        data = self._msi.data
        if data is not None:
            data_matrix_size = data.nbytes
            print("\n--- Data Matrix part ---")
            print(f"Data matrix: {data_matrix_size} bytes ({data_matrix_size / (1024 * 1024):.2f} MB)")

        # Queue section
        print("\n--- Queue part ---")
        queue_total_size = 0

        for module in self._msi.get_queue():
            # Calculate the actual memory usage of each module
            queue_total_size += sys.getsizeof(module)

        print(f"Queue size: {queue_total_size} bytes ({queue_total_size / (1024 * 1024):.2f} MB)")

        # Total
        total_size = metadata_total_size + queue_total_size + data_matrix_size
        print("\n================ Sum ================")
        print(f"sum usage: ({total_size / 1024:.2f} KB, {total_size / (1024 * 1024):.2f} MB)")

        print("================ tracker ==================")
        self.me_tr.print_diff(self.me_tr.summaries['start'], self.me_tr.summaries['after_load_data'])

    def memory_tracker_init(self):
        self.me_tr = tracker.SummaryTracker(ignore_self=True)
        self.me_tr.store_summary('start')
