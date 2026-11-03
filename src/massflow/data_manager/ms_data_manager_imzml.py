"""
MSI Data Manager for imzML format.

Handles reading .imzML files and loading data into the MSI domain model.
Supports coordinate-based sub-region loading and a shared m/z list (shared m/z axis) for
continuous data.

Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Generator

import numpy as np
from pyimzml.ImzMLParser import ImzMLParser

from massflow.data_manager.ms_data_manager import MSDataManager
from massflow.module.ms_meta_data_imzml import MetaDataImzMl
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.tools.funs import is_valid_file
import massflow.tools.read_flat as flat_tools
import massflow.tools.read_matrix as matrix_tools
from massflow.tools.logger import get_logger
from massflow.tools.stream_imzml_writer import ImzMLWriter

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
        """
        super().__init__(ms=ms,
                         target_mz_range=None,
                         target_locs=target_locs,
                         filepath=filepath,
                         mz_dtype=mz_dtype,
                         intensity_dtype=intensity_dtype,
                         temp_dir=temp_dir)

        # Initialize ImzML parser and reader
        self._parser = None
        self._reader = None
        self._writer: ImzMLWriter | None = None
        self.swap_filepath = self.file_base_path if filepath is None else None
        self.ms.meta  = MetaDataImzMl()

    @property
    def imzml_filepath(self) -> str:
        """Return the XML metadata path."""
        return f"{self.file_base_path}.imzML"

    @property
    def ibd_filepath(self) -> str:
        """Return the paired binary-data path."""
        return f"{self.file_base_path}.ibd"

    def _require_writer(self) -> ImzMLWriter:
        if self.ms.meta is None:
            raise ValueError("MS meta data is None. Please load or copy meta data first.")
        if self._writer is None:
            self._writer = ImzMLWriter.from_ms_meta(
                output_filename=str(self.file_base_path),
                meta=self.ms.meta,
                mz_dtype=self.mz_dtype,
                intensity_dtype=self.intensity_dtype,
            )
        return self._writer

    def add_spectrum(self, mz, intensity, coordinates) -> None:
        """Write one spectrum to imzML/ibd."""
        self._require_writer().add_spectrum(mz, intensity, coordinates)

    def swap_matrix_data_out2disk(
        self,
        mz_data: np.ndarray | None,
        intensity_matrix: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
    ) -> None:
        """Write a matrix batch to imzML/ibd."""
        if mz_data is None:
            raise ValueError("mz_data is required for imzML output.")
        if not np.issubdtype(coordinates.dtype, np.integer):
            coordinates = coordinates.astype(np.int32, copy=False)
        is_shared_mz = mz_data.ndim == 1
        for row, (intensity_row, length, coord) in enumerate(zip(
            intensity_matrix,
            lengths,
            coordinates,
        )):
            length = int(length)
            mz = mz_data[:length] if is_shared_mz else mz_data[row, :length]
            self.add_spectrum(
                mz,
                intensity_row[:length],
                tuple(coord.tolist()),
            )

    def write_flat_batch(
        self,
        intensity_flat: np.ndarray,
        lengths: np.ndarray,
        coordinates: np.ndarray,
        mz_axis: np.ndarray | None,
        user_params_list: list[list] | None = None,
    ) -> None:
        """Write a flat batch to imzML/ibd."""
        if len(lengths) == 0:
            return
        if (
            intensity_flat.ndim != 1
            or lengths.ndim != 1
            or (mz_axis is not None and mz_axis.ndim != 1)
        ):
            raise ValueError("intensity_flat, lengths, and mz_axis must be 1-D.")
        if not np.issubdtype(lengths.dtype, np.integer):
            lengths = lengths.astype(np.int64, copy=False)
        self._require_writer().swap_flat_out2disk(
            mz_flat=mz_axis,
            intensity_flat=intensity_flat,
            lengths=lengths,
            coords=coordinates,
            user_params_list=user_params_list,
        )

    def close_writer(self) -> None:
        """Close the imzML writer."""
        writer = getattr(self, "_writer", None)
        if writer is not None:
            writer.close()
            self._writer = None

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
                - If no pixels match, a warning is logged and the finalize hooks are skipped.
                - For continuous data, `ms.shared_mz_list` is populated once (on first seen spectrum)
                    and each placeholder gets `spectrum.mz_list = ms.shared_mz_list`.
        """
        # Start loading process
        logger.info(f"Loading data from ImzML file: {self.imzml_filepath}")


        # pre loading part : meta data load part
        self.preload_hook()

        # spectrum data load and create part  [find by coordinate]
        total_pixels = len(self.parser.coordinates)
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
        if self.current_spectrum_num == 0:
            logger.warning(
                f"target_locs {self.target_locs} matched 0 of {total_pixels} pixels in file; no spectra loaded, skipping finalize hooks."
            )
            return
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
        parser = getattr(self, "_parser", None)
        if parser is not None:
            # Try closing memory-mapped handle if present
            m = getattr(parser, "m", None)
            if m is not None and hasattr(m, "close"):
                m.close()

        self.close_writer()
        super().close()
        swap_filepath = getattr(self, "swap_filepath", None)
        if swap_filepath is not None:
            swap_base = Path(swap_filepath)
            for suffix in (".imzML", ".ibd", ".fdata", ".pdata", ".log", ".metadata", ".IBD"):
                swap_base.with_suffix(suffix).unlink(missing_ok=True)

    def copy_meta(self, source_dm: MSDataManager):
        """Copy metadata from another data manager, dropping fields that may change here.

        MetaDataBase.copy matches by field name, so cross-type sources
        (e.g. MetaDataZarr from a spectra Zarr input) copy their shared
        fields; source-only fields are ignored (LOGIC-01).
        """
        self.ms.meta.copy(source_dm.ms.meta)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _prepare_direct_read_context(
        self,
    ) -> tuple[matrix_tools.IbdFileLayout, bool, np.ndarray | None, np.ndarray]:
        """Build and return common context shared by matrix and flat generators."""
        layout = matrix_tools.extract_ibd_layout(
            self.parser,
            self.ibd_filepath,
            out_intensity_dtype=np.dtype(self.intensity_dtype),
            out_mz_dtype=np.dtype(self.mz_dtype),
        )
        is_continuous: bool = bool(self.ms.meta and self.ms.meta.continuous)
        shared_mz: np.ndarray | None = (
            np.asarray(self.ms.shared_mz_list, dtype=self.mz_dtype) if is_continuous else None
        )
        target_indices = matrix_tools.filter_target_indices(layout.coordinates, self.target_locs)

        return layout, is_continuous, shared_mz, target_indices

    def matrix_generator(
        self,
        batch_size: int = 256,
        include_mz: bool = True,
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
        layout, is_continuous, shared_mz, target_indices = self._prepare_direct_read_context()
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
                self.read_workers,
            )
            yield (mz_data, intensity, lengths, coords)

    def flat_generator(
        self,
        batch_size: int = 256,
        include_mz: bool = True,
    ) -> Generator[tuple, Any, None]:
        """
        Read data directly from .ibd into a flat 1D array.

        Args:
            batch_size: The number of spectra per batch.
            include_mz: For Processed data only, True returns mz_flat, False returns None for mz.

        Yields:
            (mz_data, intensity_flat, lengths, coordinates)

            mz_data varies by mode:
                - Continuous → shared_mz (1-D, not flattened per spectrum)
                - Processed + include_mz → mz_flat (1-D flat array)
        """
        layout, is_continuous, shared_mz, target_indices = self._prepare_direct_read_context()
        total = len(target_indices)
        if total == 0:
            return

        # Batch iteration
        for batch_start in range(0, total, batch_size):
            batch_end = min(batch_start + batch_size, total)
            batch_idx = target_indices[batch_start:batch_end]

            # Batch metadata
            lengths, coords, total_elements = flat_tools.batch_lengths_and_coords_flat(layout, batch_idx)

            # Prepare mz container and check contiguous for continuous mode
            if is_continuous:
                mz_data = shared_mz
                actual_bytes = (layout.intensity_offsets[batch_idx[-1]] + layout.intensity_lengths[batch_idx[-1]] * layout.intensity_element_bytes) - layout.intensity_offsets[batch_idx[0]]
                if actual_bytes == total_elements * layout.intensity_element_bytes:
                    intensity_flat = flat_tools.read_contiguous_block_flat(layout, batch_idx, total_elements)
                    yield (mz_data, intensity_flat, lengths, coords)
                    continue
            else:
                mz_data = np.empty(total_elements, dtype=self.mz_dtype) if include_mz else None

            # Calculate offsets for fragmented reading
            offsets = np.empty(len(lengths), dtype=np.int64)
            offsets[0] = 0
            if len(lengths) > 1:
                np.cumsum(lengths[:-1], dtype=np.int64, out=offsets[1:])

            # Fragmented multi-threaded reading
            intensity_flat = np.empty(total_elements, dtype=self.intensity_dtype)
            flat_tools.read_fragmented_block_flat(
                layout,
                batch_idx,
                intensity_flat,
                mz_data if (not is_continuous and include_mz) else None,
                offsets,
                self.threads_executor,
                self.read_workers,
            )
            yield (mz_data, intensity_flat, lengths, coords)

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
