from typing import Optional, Sequence
import numpy as np
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager import MSDataManager
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.tools.logger import get_logger
from massflow.preprocess.batch_pre_fun import BatchPreprocess
from massflow.preprocess.helper.peak_align_helper import compute_reference
from massflow.r_preprocess.adapter import CardinalAdapter

logger = get_logger("dm_pre_fun")

class Preprocess:
    """
    Data Manager Preprocessing Functions
    """
    @staticmethod
    def _total_batches(data_manager: MSDataManager,
                       batch_size: int) -> int:
        """
        Calculate total number of batches for given data manager and batch size.
        """
        total_spectra = len(data_manager.ms)
        total_batches = (total_spectra + batch_size - 1) // batch_size
        logger.info(f"Starting process: {total_spectra} spectra in {total_batches} batches")
        return total_batches

    @staticmethod
    def _log_progress(total_batches: int,
                           batch_idx: int) -> None:
        """
        Calculate number of already processed batches based on total batches and batch index.
        """
        if total_batches <= 10 or batch_idx % max(1, total_batches // 10) == 0 or batch_idx == total_batches:
            progress = (batch_idx / total_batches) * 100
            logger.info(f"process progress: {batch_idx}/{total_batches} batches ({progress:.1f}%)")

    @staticmethod
    def peak_align(data_manager: MSDataManager,
                   ref: Optional[np.ndarray] = None,
                   tolerance: Optional[float] = None,
                   units: str = 'ppm',
                   binfun: str = 'median',
                   binratio: int = 2,
                   backend_method: str = "python",
                   batch_size: int = 256,
                   clear_memory: bool = True,
                   temp_dir: str = "./temp_align_data"
                   ) -> MSDataManagerImzML:
        """
        Align peaks across spectra in MSDataManager data using specified backend.

        This method provides a unified interface for peak alignment, supporting both
        Python-based implementation and R-based Cardinal implementation. It can align
        an entire dataset (MSDataManager).

        Parameters:
            data_manager (MSDataManager, optional): The data manager containing the mass spectra to align.
            ref (np.ndarray, optional): External reference m/z axis.
                If None, it will be estimated from the data (for data_manager).
            units (str): Units for tolerance and resolution ('ppm' or 'mz'). Default is 'ppm'.
            tolerance (float, optional): The tolerance window for peak matching.
                If None, it will be estimated from the data.
            binfun (str): Aggregation function for estimating resolution ('median', 'min', 'max', 'mean').
            binratio (int): Ratio to scale the estimated resolution to determine tolerance. Default is 2.
            backend_method (str, optional): The backend to use for alignment.
                - 'cardinal': Use the R Cardinal package (requires R environment).
                - 'python' (or None): Use the native Python implementation.

        Returns:
            MSDataManagerImzML: The aligned data manager.

        Raises:
            ValueError: If `data_manager` is not provided.
        """

        if data_manager is None:
            raise ValueError("data_manager must be provided for peak alignment.")

        logger.info(
            "peak_align_entry: backend=%s, binfun=%s, tolerance=%s, units=%s",
            backend_method,
            binfun,
            tolerance,
            units,
        )

        if backend_method == "cardinal" and isinstance(data_manager, MSDataManagerImzML):
            return CardinalAdapter.align(data_manager=data_manager,
                                         reference=ref,
                                         tolerance=tolerance,
                                         units=units,
                                         binfun=binfun,
                                         binratio=binratio,
                                         temp_dir=temp_dir)

        aligned_ms = MassSpectrumSet()
        aligned_data_manager = MSDataManagerImzML(aligned_ms, temp_dir=temp_dir)
        aligned_data_manager.copy_meta(data_manager)

        writer = aligned_data_manager.writer

        if ref is None or tolerance is None:
            ref, tolerance = compute_reference(data_manager=data_manager,
                                               reference=ref,
                                               binfun=binfun,
                                               binratio=binratio,
                                               tolerance=tolerance,
                                               units=units,
                                               batch_size=batch_size,
                                               clear_memory=clear_memory)

            tolerance = tolerance * 1e6 if units == "ppm" else tolerance

        total_batches = Preprocess._total_batches(data_manager, batch_size)

        for batch_idx, batch in enumerate(data_manager.get_batch_generator(batch_size=batch_size), start=1):
            aligned_batch = BatchPreprocess.peak_align_batch(batch_spectra=batch,
                                                             ref=ref,
                                                             tolerance=tolerance,
                                                             units=units)
            data_manager.clear_batch_data_memory(batch=batch)
            aligned_data_manager.swap_batch_data_out2disk(batch=aligned_batch, writer=writer)

            Preprocess._log_progress(total_batches, batch_idx)

        aligned_data_manager.close_writer()
        aligned_data_manager.load_full_data_from_file()

        return aligned_data_manager

    @staticmethod
    def peak_pick(data_manager: MSDataManager,
                  width: int | Sequence[int] = 2,
                  method: str = 'scipy',
                  relheight: float = 0.1,
                  snr: float = 3.0,
                  return_type: str = 'height',
                  backend_method: str = "python",
                  batch_size: int = 256,
                  temp_dir: str = "./temp_pick_data"
                  ) -> MSDataManagerImzML:
        """
        Perform peak picking on MSDataManager data using specified backend.

        This method provides a unified interface for peak picking, supporting both
        Python-based implementation and R-based Cardinal implementation. It can process
        an entire dataset (MSDataManager).
        Parameters:
            data_manager (MSDataManager): The data manager containing the mass spectra to process.
            method (str): The peak picking method to use ('diff', 'sd', 'mad', 'quantile', 'filter', 'cwt', 'scipy').
                Default is 'scipy'.
            snr (float): Signal-to-noise ratio threshold for peak detection. Default is 3.0.
            return_type (str): Type of peak representation to return ('height' or 'area'). Default is 'height'.
            backend_method (str, optional): The backend to use for peak picking.
                - 'cardinal': Use the R Cardinal package (requires R environment).
                - 'python' (or None): Use the native Python implementation.

        Returns:
                MSDataManagerImzML: The picked data manager.

            Raises:
                ValueError: If `data_manager` is not provided.

        """

        if data_manager is None:
            raise ValueError("data_manager must be provided for peak picking.")

        logger.info(
            "peak_pick_entry: backend=%s, method=%s, snr=%s, return_type=%s",
            backend_method,
            method,
            snr,
            return_type,
        )

        if backend_method == "cardinal" and isinstance(data_manager, MSDataManagerImzML):
            if method != "scipy":
                return CardinalAdapter.peak_pick(data_manager=data_manager,
                                                 method=method,
                                                 snr=snr,
                                                 return_type=return_type,
                                                 temp_dir=temp_dir)
            else:
                logger.warning("Cardinal backend does not support 'scipy' method. Falling back to Python implementation.")

        picked_ms = MassSpectrumSet()
        picked_data_manager = MSDataManagerImzML(picked_ms, temp_dir=temp_dir)
        picked_data_manager.copy_meta(data_manager)

        writer = picked_data_manager.writer

        total_batches = Preprocess._total_batches(data_manager, batch_size)

        for batch_idx, batch in enumerate(data_manager.get_batch_generator(batch_size=batch_size), start=1):
            picked_batch = BatchPreprocess.peak_pick_batch(batch_spectra=batch,
                                                           width=width,
                                                           method=method,
                                                           relheight=relheight,
                                                           return_type=return_type)
            data_manager.clear_batch_data_memory(batch=batch)
            picked_data_manager.swap_batch_data_out2disk(batch=picked_batch, writer=writer)

            Preprocess._log_progress(total_batches, batch_idx)

        picked_data_manager.close_writer()
        picked_data_manager.load_full_data_from_file()

        return picked_data_manager
