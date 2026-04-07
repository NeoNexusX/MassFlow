"""
MSI Data Manager for imzML format.

Handles reading .imzML files and loading data into the MSI domain model.
Supports coordinate-based sub-region loading and a shared m/z list (shared m/z axis) for
continuous data.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
import os
from contextlib import nullcontext
from concurrent.futures import as_completed
from typing import Generator,Any
import numpy as np
from pyimzml.ImzMLParser import ImzMLParser
from massflow.data_manager.ms_data_manager import MSDataManager
from massflow.module.spectrum import Spectrum
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.ms_meta_data_imzml import MetaDataImzMl
from massflow.tools.funs import is_valid_file
import massflow.tools.read_matrix as matrix_tools
from massflow.tools.logger import get_logger

logger = get_logger("massflow.datamanager")


class MSDataManagerImzML(MSDataManager):
    """
    MSI Data Manager for .imzML files.

    Handles reading .imzML files, building MS objects with lazily-loaded
    spectra, optional spatial sub-region selection via `target_locs`, and
    metadata extraction into `ImzMlMetaData`.
    """

    def __init__(
        self,
        ms=None,
        target_locs=None,
        filepath=None,
        max_threads: int = 8,
        temp_dir = None,
        mz_dtype = np.float64,
        intensity_dtype =  np.float32,
    ):
        """
        Initialize the ImzML data manager.


        Args:
            ms (MassSpectrumSet): Mass-spectrometry domain model instance.
            target_locs (optional): Spatial window for sub-region loading.
                Expected shape: ([x1, y1], [x2, y2]) (inclusive bounds). When None, load all.
            filepath (str, optional): Path to the .imzML file (XML header; paired .ibd is
                resolved automatically).
            max_threads (int): Maximum worker threads for batch loading utilities.
        """
        super().__init__(ms=ms,
                         target_mz_range=None,
                         target_locs=target_locs,
                         filepath=filepath,
                         max_threads=max_threads,
                         mz_dtype=mz_dtype,
                         intensity_dtype=intensity_dtype,
                         temp_dir=temp_dir)

        # Initialize ImzML parser and reader
        self._parser = None
        self._reader = None
        self.ms.meta  = MetaDataImzMl()

    @property
    def reader(self):
        """Get the ImzML spectrum reader."""
        if self._reader is None:
            if self.parser is not None:
                self._reader = self.parser.portable_spectrum_reader()
            else:
                logger.error("ImzML parser is not initialized. Cannot create reader.")
                raise ValueError("ImzML parser is not initialized. Cannot create reader.")
        return self._reader

    @property
    def parser(self):

        if self._parser is not None:
            return self._parser

        if is_valid_file(self.imzml_filepath) and is_valid_file(self.ibd_filepath):
            self._parser = ImzMLParser(self.imzml_filepath)
            logger.info(f"Initialized ImzML parser for file: {self.imzml_filepath}")
            return self._parser
        else:
            logger.error(f"ImzML file {self.imzml_filepath} or {self.ibd_filepath} does not exist or is empty.")
            raise FileNotFoundError(f"ImzML file {self.imzml_filepath} or {self.ibd_filepath} does not exist or is empty.")

    def load_head_data(self):
        """
        Stream spectra from .imzML and add lazy SpectrumImzML placeholders into `ms`.

        Behavior:
                - Iterates over `parser.coordinates` and creates a `SpectrumImzML` placeholder per pixel.
                - If `target_locs` is provided, only pixels within the inclusive window are added.
                - For continuous data, `ms.shared_mz_list` is populated once (on first seen spectrum)
                    and each placeholder gets `spectrum.mz_list = ms.shared_mz_list`.
        """
        # Start loading process
        logger.info(f"Loading data from ImzML file: {self.imzml_filepath}")


        # pre loading part : meta data load part
        self.preload_hook()

        # spectrum data load and create part  [find by coordinate]
        for i, c in enumerate(self.parser.coordinates):
            x, y, z = c
            if self.target_locs[0][0] <= x <= self.target_locs[1][0] and self.target_locs[0][1] <= y <= self.target_locs[1][1]:

                # Build (x,y,z)->index mapping and add SpectrumImzML placeholders
                spectrum = SpectrumImzML(reader=self.reader,
                                            ibd_path=self.ibd_filepath,
                                            index=i,
                                            coordinates=[x, y, z],
                                            shared_mz_list=self.ms.shared_mz_list)

                self.ms.add_spectrum(spectrum)
                self.current_spectrum_num += 1

        # loaded part : finalize meta data load part
        self.loaded_hook()

    def extract_metadata(self):
        """Extract parser metadata and batch-populate `MetaDataImzMl` fields."""

        logger.info("Extracting metadata...")

        if self.parser is None:
            logger.error("Parser is not initialized. Please set parser or filepath first.")
            raise ValueError("Parser is not initialized. Please set parser or filepath first.")

        if not isinstance(self.ms.meta, MetaDataImzMl):
            logger.error("ms.meta is not MetaDataImzMl; cannot populate from ImzML parser.")
            raise TypeError("ms.meta is not MetaDataImzMl; cannot populate from ImzML parser.")

        self.ms.meta.populate_from_parser(self.parser)

        logger.info("Metadata extraction completed.")

    def close(self):
        """
        Safely close underlying resources associated with the ImzML parser.
        """
        if self.parser is not None:
            # Try closing memory-mapped handle if present
            m = getattr(self.parser, "m", None)
            if m is not None and hasattr(m, "close"):
                m.close()

        super().close()

    def copy_meta(self, source_dm: MSDataManager):
        # 检查 source_dm.ms.meta 是否也是 MetaDataImzMl 类型（或其子类）
        if isinstance(source_dm.ms.meta, type(self.ms.meta)):
            self.ms.meta.copy(source_dm.ms.meta)
        else:
            logger.warning("Source and target metadata types do not match.")

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def batch_generator(self, batch_size: int = 256, max_threads: int = 0) -> Generator[list[Spectrum], None, None]:
        """Get a multi-threaded batch generator that yields lists of Spectrum objects."""
        # update max_threads if needed
        if max_threads != 0 and max_threads != self.max_threads:
            self.max_threads = max_threads

        total_length = len(self.ms)
        has_ibd_file = hasattr(self, "ibd_filepath") and os.path.exists(self.ibd_filepath)

        # Nested function for processing a mini-batch
        def process_mini_batch(mini_batch, file_handle=None):
            """Process a mini-batch with optional resolve kwargs."""
            with  nullcontext(file_handle) if file_handle is  None else file_handle as f:
                for sp in mini_batch:
                    sp.resolve_data(file_handle=f)

        # Batch processing logic
        for i in range(0, total_length, batch_size):
            batch = self.ms[i : i + batch_size]  # Extract the current batch

            if self.threads_executor is not None:
                # Split the batch into mini-batches based on the number of threads
                mini_batch_size = max(1, len(batch) // self.max_threads)
                mini_batches = [batch[j : j + mini_batch_size] for j in range(0, len(batch), mini_batch_size)]

                # Multi-threaded processing of each mini-batch
                futures = []
                for mini_batch in mini_batches:
                    file_handle= open(self.ibd_filepath, "rb")  if has_ibd_file else None

                    future = self.threads_executor.submit(
                        process_mini_batch,
                        mini_batch,
                        file_handle)
                    futures.append(future)

                # Wait for all mini-batch tasks to complete
                for future in as_completed(futures):
                    future.result()  # Ensure all data in the mini-batch has been processed

            # Yield the fully processed batch
            yield batch

    def matrix_generator(
        self,
        batch_size: int = 256,
        include_mz: bool = True,
        max_threads: int = 0,
    ) -> Generator[tuple, Any, None]:
        """
        Read data directly from .ibd into pre-allocated NumPy matrices.

        Args:
            include_mz: For Processed data only, True returns mz_matrix, False returns None.

        Yields:
            (mz_data, intensity_matrix, lengths, coordinates)

            mz_data varies by mode:
                - Continuous → shared_mz (1-D)
                - Processed + include_mz → mz_matrix (2-D)
        """
        # Thread pool management
        if max_threads != 0 and max_threads != self.max_threads:
            self.max_threads = max_threads

        # Extract file layout in one pass
        layout = matrix_tools.extract_ibd_layout(
            self.parser,
            self.ibd_filepath,
            out_intensity_dtype=np.dtype(self.intensity_dtype),
            out_mz_dtype=np.dtype(self.mz_dtype),
        )

        # Determine data mode
        is_continuous: bool = bool(self.ms.meta and self.ms.meta.continuous)
        shared_mz: np.ndarray | None = (np.asarray(self.ms.shared_mz_list, dtype=self.mz_dtype) if is_continuous else None)

        # Filter target spectrum indices
        target_indices = matrix_tools.filter_target_indices(layout.coordinates, self.target_locs)
        total = len(target_indices)
        if total == 0:
            return

        # Batch iteration
        for batch_start in range(0, total, batch_size):
            batch_end = min(batch_start + batch_size, total)
            batch_idx = target_indices[batch_start:batch_end]
            n = len(batch_idx)

            # Batch metadata
            lengths, max_len, coords = matrix_tools.batch_lengths_and_coords(layout, batch_idx)

            # Prepare mz container
            if is_continuous:
                mz_data = shared_mz
            else:
                mz_data = np.zeros((n, max_len), dtype=self.mz_dtype) if include_mz else None

            # Continuous mode + disk contiguous
            if (is_continuous and n > 1 and matrix_tools.is_disk_contiguous(layout, batch_idx, max_len)):
                intensity = matrix_tools.read_contiguous_block(layout, batch_idx, n, max_len)
                yield (mz_data, intensity, lengths, coords)
                continue

            # Fragmented multi-threaded reading
            intensity = np.zeros((n, max_len), dtype=self.intensity_dtype)
            matrix_tools.read_fragmented_block(
                layout,
                batch_idx,
                intensity,
                mz_data if (not is_continuous and include_mz) else None,
                self.threads_executor,
                self.max_threads,
            )
            yield (mz_data, intensity, lengths, coords)

####### hook part #########
    def preload_hook(self, *args, **kwargs):
        """Pre-load and populate metadata from the ImzML parser before loading spectra."""
        #load meatdata from parser and assign to ms.meta
        self.extract_metadata()

        # assign shared mz_list for continuous data
        if self.parser and self.ms.meta.continuous and self.ms.shared_mz_list is None:  # type: ignore
            self.ms.shared_mz_list, _ = self.parser.getspectrum(0)
            logger.info(f"Assigning shared mz_list for continuous data; len is {len(self.ms.shared_mz_list)} subsequent spectra will reuse it.")

    def loading_hook(self, *args, **kwargs):
        """Hook invoked while spectra are being loaded."""

    def loaded_hook(self):
        """Finalize metadata after loading spectra."""
        logger.info("creating ms mask.")
        self.create_ms_meta_mask()
        self.detect_sorted()
        self.inspect_data()
####### hook end  #########
