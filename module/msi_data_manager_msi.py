import glob
import sys
import os
from pympler import tracker
import h5py
import numpy as np
from .msi_module import MSI, MSIBaseModule
from .msi_data_manager import MSIDataManager
from logger import get_logger
logger = get_logger("msi_data_manager_msi")


class MSIDataManagerMSI(MSIDataManager):
    """
    MSI Data Manager for .msi Data

    Handles batch read/write of .h5/.msi files, filters data by m/z range,
    and provides memory statistics, visualization, and split/merge output.
    """
    def __init__(self,
                 msi: MSI,
                 target_mz_range=None,
                 filepath=None):

        super().__init__(msi, target_mz_range, filepath)
        self.me_tr = None
        # self.memory_tracker()

    def get_msi(self) -> MSI:
        """
        Get the MSI object.

        Returns:
            MSI: The MSI object.
        """
        return self._msi

    #file action
    def load_full_data_from_file(self):
        """
        Load MSI data completely from files (metadata + all m/z images).

        Steps:
        1. Call load_data_helper to read all metadata.
        2. Preallocate a 3D data matrix using meta.mz_num and mask dimensions.
        3. Call load_data_helper again to read all image data and fill the matrix.

        Notes:
        - Requires meta.mz_num > 0, otherwise an assertion error is raised.
        - If target_mz_range is set, only load m/z images within that range.
        """

        # Read metadata
        self.__load_data_helper(fn_name='meta')
        logger.info(f"MSI meta data: {self.get_msi().meta}")

        if self.get_msi().meta.mz_num <1 :
            logger.error("meta.mz_num must be greater than 0")
            raise ValueError("meta.mz_num must be greater than 0")

        # Initialize data storage matrix via MSI API
        self.get_msi().allocate_data_from_meta(dtype=np.float32)
        self.__load_data_helper(fn_name='data')
        # self.me_tr.store_summary('after_load_data')

    def __load_data_helper(self, fn_name = 'meta'):
        # This method should not be used by son class its only design for this class
        # Find all .h5 files in the specified directory
        fn = self.__load_meta_from_file if fn_name == 'meta' else self.__load_data_from_file

        if (self.filepath.endswith('.h5') or self.filepath.endswith('.msi') )and os.path.isfile(self.filepath):
            fn(self.filepath)
        elif os.path.isdir(self.filepath):
            msi_files = glob.glob(os.path.join(self.filepath, "*.msi"))
            for msi_file in msi_files:
                fn(msi_file)
        else:
            logger.error(f"Error: {self.filepath} is not a valid .h5 file or .msi file or directory.")
            raise ValueError(f"Error: {self.filepath} is not a valid .h5 file or .msi file or directory.")

        if fn_name != 'meta':
            # Check if current_num exceeds meta.mz_num after loading data
            if self.current_image_num < self._msi.meta.mz_num:
                logger.info(f"current_num {self.current_image_num} < meta.mz_num {self._msi.meta.mz_num}, "
                               f"some m/z images may be missing.")
            #update meta.mz_num if current_num is smaller
            if self.target_mz_range is not None and self.current_image_num <= self._msi.meta.mz_num:
                self._msi.meta.mz_num = self.current_image_num

    def __load_meta_from_file(self, file):
        """
        Load metadata from an .msi file.

        Parameters:
        - file (str): Path to the HDF5 file.

        Notes:
        - Reads metadata groups starting with 'meta_' and stores them in the MSI object.
        - Decodes byte strings to UTF-8 if necessary.
        - The function that actually performs the meta reading, called with the help of __load_data_helper.
        """
        with h5py.File(file, 'r') as h5_file:
            for key, group in h5_file.items():
                if key.startswith('meta_') and self._msi.meta.get(key) is None:
                    attr = key.replace('meta_', '', 1)
                    if self._msi.meta.get(attr) is None:
                        # Read raw dataset value
                        dataset_value = group[()]
                        # Decode bytes to string if necessary
                        if isinstance(dataset_value, bytes):
                            dataset_value = dataset_value.decode('utf-8')
                        setattr(self._msi.meta, attr, dataset_value)

    def __load_data_from_file(self, file):
        """
        Load MSI data images from an .msi file.

        Parameters:
        - file (str): Path to the HDF5 file.

        Notes:
        - Reads m/z images groups starting with 'mz_' and stores them in the MSI object.
        - Decodes byte strings to UTF-8 if necessary.
        - The function that actually performs the data reading, called with the help of __load_data_helper.
        """
        with h5py.File(file, 'r') as h5_file:
            for key, group in h5_file.items():
                if isinstance(group, h5py.Group) and not key.startswith('meta_'):
                    mz= group['mz'][()]
                    if self.target_mz_range is not None and not self.target_mz_range[0] <= mz <= self.target_mz_range[1]:
                        continue
                    # Compute base_mask based on metadata if needed
                    msi_image = group['msroi'][()]
                    base_mask = np.where(msi_image > 0, 1, 0) if self._msi.meta.need_base_mask else None

                    self._msi.data[self.current_image_num, :, :] = msi_image
                    self._msi.add_msi_img(
                        MSIBaseModule(mz=mz,
                                      msroi=(self._msi.data[self.current_image_num] if self._msi.data is not None else msi_image),
                                      base_mask=base_mask
                                      )
                    )
                    logger.info(f'loading {key}')
                    self.current_image_num += 1
            logger.info(f"finish loading {file}")

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
                logger.info(f"{prefix}Group: {key}")
                self._inspect_hdf5_structure(item, indent + 1, max_depth)

            # If it is a dataset, print name and shape (skip reference type)
            elif isinstance(item, h5py.Dataset):
                # Skip reference-type datasets (dtype.kind == 'O')
                logger.info(f"{prefix}Dataset: {key}  Shape: {item.shape}  type: {item.dtype}")

    def calculate_memory_usage(self):
        """
        Calculate and print the memory usage of the current MSI object.

        Separately counts memory usage of metadata, data matrices, and modules in the queue,
        then summarizes the total usage and outputs in KB/MB units.
        """

        logger.info("=== memory usage ===")

        # Data matrix section
        data_matrix_size = 0
        data = self._msi.data
        if data is not None:
            data_matrix_size = data.nbytes
            logger.info("\n--- Data Matrix part ---")
            logger.info(f"Data matrix: {data_matrix_size} bytes ({data_matrix_size / (1024 * 1024):.2f} MB)")

        # Queue section
        logger.info("\n--- Queue part ---")
        queue_total_size = 0

        for module in self._msi.get_queue():
            # Calculate the actual memory usage of each module
            queue_total_size += sys.getsizeof(module)

        logger.info(f"Queue size: {queue_total_size} bytes ({queue_total_size / (1024 * 1024):.2f} MB)")

        # Total
        total_size = queue_total_size + data_matrix_size
        logger.info("\n================ Sum ================")
        logger.info(f"sum usage: ({total_size / 1024:.2f} KB, {total_size / (1024 * 1024):.2f} MB)")

        logger.info("================ tracker ==================")
        self.me_tr.print_diff(self.me_tr.summaries['start'], self.me_tr.summaries['after_load_data'])

    def memory_tracker_init(self):
        """
        Initialize the memory tracker for the MSI object.

        Sets up a SummaryTracker to monitor memory usage,
        storing a summary at the start of loading data.
        """
        self.me_tr = tracker.SummaryTracker(ignore_self=True)
        self.me_tr.store_summary('start')
