"""
MS Data Management Module

Provides abstract interfaces and common utilities for managing MS/MSI data,
including metadata inspection, batch spectrum loading, and occupancy mask
creation over spatial coordinates.
Provides abstract interfaces and common utilities for managing MS/MSI data,
including metadata inspection, batch spectrum loading, and occupancy mask
creation over spatial coordinates.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""

import os
import uuid
from abc import ABC, abstractmethod
import random
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from massflow.logger import get_logger
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.tools.stream_imzml_writer import ImzMLWriter

logger = get_logger("ms_data_manager")



class MSDataManager(ABC):
    """
    Abstract base class for MS Data Manager.
    """

    def __init__(
        self,
        ms=None,
        target_mz_range=None,
        target_locs=None,
        filepath=None,
        temp_dir=None,
        max_threads: int = 8,
        mz_dtype=np.float64,
        intensity_dtype=np.float32,
    ):
        """Initialize the MS data manager."""
        # create MassSpectrumSet object if not provided
        self._ms = MassSpectrumSet() if ms is None else ms
        self.target_mz_range = target_mz_range

        # Multi-threading control part
        self.max_threads = max_threads
        self._threads_executor = None

        # data type part
        self.mz_dtype = mz_dtype
        self.intensity_dtype = intensity_dtype

        # swap part
        if filepath is None:
            base_dir = temp_dir if temp_dir is not None else tempfile.gettempdir()
            if not os.path.exists(base_dir):
                os.makedirs(base_dir)
            self.filepath = os.path.join(base_dir, f"temp_{uuid.uuid4().hex}.imzML")
            self.swap_filepath = self.filepath
        else:
            self.filepath = filepath

        self._writer = None

        # target_locs input verification
        if target_locs is not None:
            # target locs input
            if len(target_locs) <= 1 or (
                target_locs[0][0] > target_locs[1][0]
                or target_locs[0][1] > target_locs[1][1]
            ):
                logger.error("target_locs must be non-empty or locs must x1<x2,y1<y2")
                raise ValueError(
                    "target_locs must be non-empty or locs must x1<x2,y1<y2"
                )

        # update target_locs
        self.target_locs = target_locs
        self.current_spectrum_num = 0

    @property
    def threads_executor(self) -> ThreadPoolExecutor:
        """lazy load for thread pool executor"""
        if self._threads_executor is None:
            logger.info(
                f"Creating ThreadPoolExecutor with max_workers={self.max_threads}"
            )
            self._threads_executor = ThreadPoolExecutor(max_workers=self.max_threads)
        return self._threads_executor

    @threads_executor.setter
    def threads_executor(self, executor: ThreadPoolExecutor):
        """set thread pool executor"""
        self._threads_executor = executor

    @property
    def writer(self):
        """return writer."""

        if self.ms.meta is not None:

            if self._writer is None:
                spec_type = "profile" if self.ms.meta.profile_spectrum is not None else "centroid"
                scan_direction = self.ms.meta.scan_direction if self.ms.meta.scan_direction is not None else "top_down"
                line_scan_direction = self.ms.meta.line_scan_direction if self.ms.meta.line_scan_direction is not None else "line_left_right"
                scan_pattern = self.ms.meta.scan_pattern if self.ms.meta.scan_pattern is not None else "one_way"
                scan_type = self.ms.meta.scan_type if self.ms.meta.scan_type is not None else "horizontal_line"

                if self.ms.meta.continuous is not None:
                    mode = "continuous"
                elif self.ms.meta.processed is not None:
                    mode = "processed"
                else:
                    mode = "auto"
                logger.info(f"Creating ImzMLWriter with mode={mode}, spec_type={spec_type}")
                self._writer = ImzMLWriter(str(self.filepath),
                                            mz_dtype=self.mz_dtype, #type: ignore
                                            intensity_dtype=self.intensity_dtype,#type: ignore
                                            mode=mode,
                                            spec_type=spec_type,
                                            scan_direction=scan_direction,
                                            line_scan_direction=line_scan_direction,
                                            scan_pattern=scan_pattern,
                                            scan_type=scan_type)
            return self._writer

        else:
            logger.error(
                "MS meta data is None. Please load or copy[use copy_meta method] meta data first."
            )
            raise ValueError(
                "MS meta data is None. Please load or copy meta data first."
            )

    def copy_meta(self, source_dm: "MSDataManager"):
        """copy metadata from source MS to self MS object."""

        if source_dm.ms.meta is None or self._ms.meta is None:
            logger.warning("Source MS has no metadata to copy or meta is none.")
            return

        exclude_keys = {
            "mask",
            "filepath",
            "spectrum_count_num",
            "continuous",
            "processed",
            "centroid_spectrum",
            "profile_spectrum",
        }

        for key, value in source_dm.ms.meta.items():
            if key in exclude_keys:
                continue

            if hasattr(self._ms.meta, key):
                try:
                    setattr(self._ms.meta, key, value)
                except (AttributeError, TypeError, ValueError) as e:
                    logger.warning(f"Failed to copy meta field {key}: {e}")

        logger.info("Metadata copy completed.")

    @property
    def ms(self) -> MassSpectrumSet:
        """Get the MS object."""
        return self._ms

    @ms.setter
    def ms(self, ms: MassSpectrumSet):
        """Set the MS object."""
        self._ms = ms

    @property
    def swap_writer(self):
        """Return the writer."""
        return self.writer

    @abstractmethod
    def load_full_data_from_file(self):
        """Load data from the configured input file path."""

    def inspect_data(self, inspect_num=2):
        """Inspect the data structure of the MSI object."""
        meta_info = "MS meta data:\r\n"
        meta_info += f"  target_mz_range: {self.target_mz_range}\r\n"
        meta_info += f"  target_locs: {self.target_locs}\r\n"
        meta_info += f"  filepath: {self.filepath}\r\n"
        meta_info += f"  current_spectrum_num: {self.current_spectrum_num}\r\n"

        # get meta data info
        if self.ms.meta is not None:
            for attr, value in self.ms.meta.items():
                shape = getattr(value, "shape", None)
                if shape is not None and len(shape) > 0:
                    meta_info += f"  meta_{attr}: {shape}\r\n"
                else:
                    meta_info += f"  meta_{attr}: {value}\r\n"
            logger.info(meta_info)

        # get base info
        base_info = "MS  information:\r\n"
        pointer4num = 0
        for spectrum in self._ms:
            if pointer4num >= inspect_num:
                break
            base_info += (
                f"  MS len: {len(spectrum)}\r\n"
                f"  MS range: {min(spectrum.mz_list)} - {max(spectrum.mz_list)}\r\n"
                f"  MS coord: {spectrum.coordinate}\r\n"
                f"  max and min mz_list: {max(spectrum.mz_list)} - {min(spectrum.mz_list)}\r\n"
                f"  max intensity: {max(spectrum.intensity)}\r\n\r\n"
            )
            pointer4num += 1
        logger.info(base_info)

    def get_batch_generator_st(self, batch_size: int = 256):
        """Get a single-threaded batch generator."""
        total_length = len(self.ms)
        for i in range(0, total_length, batch_size):
            batch = self.ms[i : i + batch_size]  # Extract the current batch
            yield batch

    def get_batch_generator(self, batch_size: int = 256, max_threads: int = -1):
        """Get a multi-threaded batch generator."""
        # update max_threads if needed
        if max_threads != -1 and max_threads != self.max_threads:
            self.max_threads = max_threads
            # reload thread pool executor
            if self.threads_executor is not None:
                self.threads_executor.shutdown(wait=True)
            self.threads_executor = ThreadPoolExecutor(max_workers=self.max_threads)

        total_length = len(self.ms)

        # Nested function for processing a mini-batch
        def process_mini_batch(mini_batch):
            """Process a mini-batch."""
            for sp in mini_batch:
                _ = sp.intensity  # Trigger the loading of spectrum intensity data

        # Batch processing logic
        for i in range(0, total_length, batch_size):
            batch = self.ms[i : i + batch_size]  # Extract the current batch

            if self.threads_executor is not None:
                # Split the batch into mini-batches based on the number of threads
                mini_batch_size = max(1, len(batch) // self.max_threads)
                mini_batches = [
                    batch[j : j + mini_batch_size]
                    for j in range(0, len(batch), mini_batch_size)
                ]

                # Multi-threaded processing of each mini-batch
                futures = []
                for mini_batch in mini_batches:
                    future = self.threads_executor.submit(
                        process_mini_batch, mini_batch
                    )
                    futures.append(future)

                # Wait for all mini-batch tasks to complete
                for future in as_completed(futures):
                    future.result()  # Ensure all data in the mini-batch has been processed

            # Yield the fully processed batch
            yield batch

    def detect_sorted(self):
        """Check if the m/z arrays in the MS data are sorted."""

        total_length = len(self.ms)
        check_num = random.randint(0, total_length)
        check_rsult = self.ms[check_num].is_sorted()
        if not check_rsult:
            logger.info("Detected unsorted m/z data in the MS object.")
            for spectrum in self.ms:
                spectrum.sort_by_mz = False
        else:
            logger.info("random test m/z data in the MS object are sorted.")

    def create_ms_meta_mask(self):
        """Create a binary occupancy mask for all available spectra coordinates."""
        if self.ms.meta is None:
            logger.error("MS meta data is None. Please load meta data first.")
            raise ValueError("MS meta data is None. Please load meta data first.")

        if (
            self.ms.meta.max_count_of_pixels_x is None
            or self.ms.meta.max_count_of_pixels_y is None
        ):
            logger.error("Image dimensions missing in meta data.")
            raise ValueError("Image dimensions missing in meta data.")

        # Get image dimensions
        # Get image dimensions
        width = int(self.ms.meta.max_count_of_pixels_x)
        height = int(self.ms.meta.max_count_of_pixels_y)

        # Prepare an empty occupancy mask
        mask = np.zeros((height, width), dtype=np.int8)

        # Fill the occupancy mask across all z-planes
        for z_dict in self.ms.coordinate_index.values():
            for x, y_dict in z_dict.items():
                for y in y_dict.keys():
                    ix = int(x)
                    iy = int(y)
                    if 0 <= ix < width and 0 <= iy < height:
                        mask[iy, ix] = 1

        # Cache mask in metadata and return
        self.ms.meta.mask = mask
    
    @abstractmethod
    def preload_meta(self, *args, **kwargs):
        """Hook for pre-loading or preparing metadata before full data loading."""

    @abstractmethod
    def loading_meta(self, *args, **kwargs):
        """Hook invoked while spectra are being loaded."""

    def loaded_meta(self, *args, **kwargs):
        """Finalize metadata after loading spectra."""
        logger.info("creating ms mask.")
        self.create_ms_meta_mask()
        self.detect_sorted()

    def close(self):
        """
        Close any resources held by the data manager.
        """
        # thread pool delete
        if self._threads_executor is not None:
            self._threads_executor.shutdown(wait=True)
            self._threads_executor = None

        # temporary directory delete
        if getattr(self, "swap_filepath", None) is not None:
            if os.path.exists(self.swap_filepath):
                ibd_filepath = self.swap_filepath.replace('.imzML', '.ibd')
                fdata_filepath = self.swap_filepath.replace('.imzML', '.fdata')
                log_filepath = self.swap_filepath.replace('.imzML', '.log')
                metadata_filepath = self.swap_filepath.replace('.imzML', '.metadata')
                os.remove(self.swap_filepath)
                os.remove(ibd_filepath)
                if os.path.exists(fdata_filepath): os.remove(fdata_filepath)
                if os.path.exists(log_filepath): os.remove(log_filepath)
                if os.path.exists(metadata_filepath): os.remove(metadata_filepath)
                logger.info(f"Temporary swap file {self.swap_filepath}, {ibd_filepath} removed.")

    def close_writer(self):
        """close ImzMLWriter."""
        if hasattr(self, "writer") and self.writer is not None:
            self.writer.close()

    def clear_batch_data_memory(self, batch):
        """Clear data in a batch of spectra to free memory."""
        for spec in batch:
            spec.clear_data()

    def clear_all_data_memory(self):
        """
        Clear data in all spectra to free memory.
        """
        for spec in self.ms:
            spec.clear_data()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()
