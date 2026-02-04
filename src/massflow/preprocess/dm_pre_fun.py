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
    def _process_in_batches(
        data_manager: MSDataManager,
        batch_size: int,
        temp_dir: str,
        batch_func,
        **batch_kwargs,
    ) -> MSDataManagerImzML:
        processed_ms = MassSpectrumSet()
        processed_data_manager = MSDataManagerImzML(processed_ms, temp_dir=temp_dir)
        processed_data_manager.copy_meta(data_manager)

        writer = processed_data_manager.writer

        total_batches = Preprocess._total_batches(data_manager, batch_size)

        for batch_idx, batch in enumerate(data_manager.get_batch_generator(batch_size=batch_size), start=1):
            processed_batch = batch_func(batch_spectra=batch, **batch_kwargs)
            data_manager.clear_batch_data_memory(batch=batch)
            processed_data_manager.swap_batch_data_out2disk(batch=processed_batch, writer=writer)

            Preprocess._log_progress(total_batches, batch_idx)

        processed_data_manager.close_writer()
        processed_data_manager.load_full_data_from_file()

        return processed_data_manager

    @staticmethod
    def peak_align(data_manager: MSDataManager,
                   ref: Optional[np.ndarray] = None,
                   tolerance: Optional[float] = None,
                   units: str = 'ppm',
                   binfun: str = 'median',
                   binratio: int = 2,
                   backend: str = "python",
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
            backend (str, optional): The backend to use for alignment.
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
            f"peak_align_entry: backend={backend}, binfun={binfun}, tolerance={tolerance}, units={units}"
        )

        if backend == "cardinal" and isinstance(data_manager, MSDataManagerImzML):
            return CardinalAdapter.align(data_manager=data_manager,
                                         reference=ref,
                                         tolerance=tolerance,
                                         units=units,
                                         binfun=binfun,
                                         binratio=binratio,
                                         temp_dir=temp_dir)

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

        return Preprocess._process_in_batches(
            data_manager=data_manager,
            batch_size=batch_size,
            temp_dir=temp_dir,
            batch_func=BatchPreprocess.peak_align_batch,
            ref=ref,
            tolerance=tolerance,
            units=units,
        )

    @staticmethod
    def peak_pick(data_manager: MSDataManager,
                  width: int | Sequence[int] = 2,
                  method: str = 'origin',
                  relheight: float = 0.012,
                  snr: float = 2.0,
                  return_type: str = 'height',
                  backend: str = "python",
                  batch_size: int = 256,
                  temp_dir: str = "./temp_pick_data",
                  use_numba: bool = True
                  ) -> MSDataManagerImzML:
        """
        Perform peak picking on MSDataManager data using specified backend.

        This method provides a unified interface for peak picking, supporting both
        Python-based implementation and R-based Cardinal implementation. It can process
        an entire dataset (MSDataManager).
        Parameters:
            data_manager (MSDataManager): The data manager containing the mass spectra to process.
            method (str): The peak picking method to use ('diff', 'sd', 'mad', 'quantile', 'filter', 'cwt', 'origin').
                Default is 'origin'.
            snr (float): Signal-to-noise ratio threshold for peak detection. Default is 2.0.
            return_type (str): Type of peak representation to return ('height' or 'area'). Default is 'height'.
            backend (str, optional): The backend to use for peak picking.
                - 'cardinal': Use the R Cardinal package (requires R environment).
                - 'python' (or None): Use the native Python implementation.

        Returns:
                MSDataManagerImzML: The picked data manager.

            Raises:
                ValueError: If `data_manager` is not provided.

        """
        method = "diff" if backend == "cardinal" and method == "origin" else method

        if data_manager is None:
            raise ValueError("data_manager must be provided for peak picking.")

        logger.info(
            f"peak_pick_entry: backend={backend}, method={method}, snr={snr}, return_type={return_type}"
        )

        if backend == "cardinal" and isinstance(data_manager, MSDataManagerImzML):
            return CardinalAdapter.peak_pick(data_manager=data_manager,
                                                method=method,
                                                snr=snr,
                                                return_type=return_type,
                                                temp_dir=temp_dir)

        return Preprocess._process_in_batches(
            data_manager=data_manager,
            batch_size=batch_size,
            temp_dir=temp_dir,
            batch_func=BatchPreprocess.peak_pick_batch,
            width=width,
            method=method,
            relheight=relheight,
            snr=snr,
            return_type=return_type,
            use_numba=use_numba,
        )

    @staticmethod
    def noise_reduction(
        data_manager: MSDataManager,
        method: str = "ma",
        window: int = 5,
        sd: float | None = None,
        sd_intensity: float | None = None,
        p: int = 2,
        coef: np.ndarray | None = None,
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0,
        wavelet: str = "db4",
        threshold_mode: str = "soft",
        batch_size: int = 256,
        temp_dir: str = "./temp_noise_data",
    ) -> MSDataManagerImzML:
        """Perform noise reduction on MSDataManager data.

        This method follows the same batching pattern as :meth:`peak_pick`, but forwards
        all denoising parameters to :meth:`SpectrumPreprocess.noise_reduction_spectrum` via
        :meth:`BatchPreprocess.noise_reduction_batch`.

        Parameters:
            data_manager: Data manager containing spectra to be denoised.
            method: One of {'ma','gaussian','savgol','savgol_numba','wavelet','ma_ns','ma_ns_numba','gaussian_ns','gaussian_ns_numba','bi_ns','bi_ns_numba'}.
            window: Window size or neighbor count depending on method.
            sd: Gaussian scale parameter.
            sd_intensity: Intensity scale for bilateral method.
            p: Minkowski metric for NS queries.
            coef: Custom kernel for 'ma'.
            polyorder: Polynomial order for Savitzky-Golay.
            deriv: Derivative order for Savitzky-Golay.
            delta: Sample spacing for Savitzky-Golay.
            wavelet: Wavelet family for wavelet denoising.
            threshold_mode: 'soft' or 'hard' thresholding.
            batch_size: Number of spectra per batch.
            temp_dir: Temporary directory for writing denoised data.

        Returns:
            MSDataManagerImzML containing denoised spectra.
        """

        if data_manager is None:
            raise ValueError("data_manager must be provided for noise reduction.")

        method_norm = (method or "ma").strip().lower()
        if method_norm in {"ma", "gaussian", "savgol", "savgol_numba"}:
            logger.info(
                f"noise_reduction_entry: method={method}, window={window}, polyorder={polyorder}, sd={sd}, deriv={deriv}, delta={delta}"
            )
        elif method_norm in {"wavelet"}:
            logger.info(
                f"noise_reduction_entry: method={method}, wavelet={wavelet}, threshold_mode={threshold_mode}"
            )
        else:
            logger.info(
                f"noise_reduction_entry: method={method}, window(k)={window}, p={p}, sd={sd}, sd_intensity={sd_intensity}"
            )

        return Preprocess._process_in_batches(
            data_manager=data_manager,
            batch_size=batch_size,
            temp_dir=temp_dir,
            batch_func=BatchPreprocess.noise_reduction_batch,
            method=method,
            window=window,
            sd=sd,
            sd_intensity=sd_intensity,
            p=p,
            coef=coef,
            polyorder=polyorder,
            deriv=deriv,
            delta=delta,
            wavelet=wavelet,
            threshold_mode=threshold_mode,
        )

    @staticmethod
    def normalization(
        data_manager: MSDataManager,
        scale_method: str = "none",
        method: str = "tic",
        scale: float = 1.0,
        batch_size: int = 256,
        temp_dir: str = "./temp_normalization_data",
    ) -> MSDataManagerImzML:
        """Perform intensity normalization on MSDataManager data.

        This method mirrors :meth:`noise_reduction` but forwards normalization
        parameters to :meth:`SpectrumPreprocess.normalization_spectrum` via
        :meth:`BatchPreprocess.normalization_batch`.

        Parameters:
            data_manager: Data manager containing spectra to be normalized.
            scale_method: 'none' or 'unit' min-max scaling.
            method: One of {'tic', 'rms', 'median'}.
            scale: Cardinal-like amplitude scaling factor applied after normalization.
            batch_size: Number of spectra per batch.
            temp_dir: Temporary directory for writing normalized data.

        Returns:
            MSDataManagerImzML containing normalized spectra.
        """
        if data_manager is None:
            raise ValueError("data_manager must be provided for normalization.")

        logger.info(
            f"normalization_entry: method={method}, scale_method={scale_method}, scale={scale}"
        )

        return Preprocess._process_in_batches(
            data_manager=data_manager,
            batch_size=batch_size,
            temp_dir=temp_dir,
            batch_func=BatchPreprocess.normalization_batch,
            scale_method=scale_method,
            method=method,
            scale=scale,
        )

    @staticmethod
    def baseline_correction(
        data_manager: MSDataManager,
        method: str = "asls",
        smooth: str = "none",
        span: float = 0.1,
        s: float | None = 0.0,
        upper: bool = False,
        width: int = 5,
        lam: float = 1e7,
        p: float = 0.01,
        niter: int = 15,
        baseline_scale: float = 1.0,
        m: int | None = None,
        decreasing: bool = True,
        batch_size: int = 256,
        temp_dir: str = "./temp_baseline_data",
    ) -> MSDataManagerImzML:
        """Perform baseline correction on MSDataManager data.

        This method mirrors :meth:`noise_reduction` but forwards baseline parameters to
        :meth:`SpectrumPreprocess.baseline_correction_spectrum` via
        :meth:`BatchPreprocess.baseline_correction_batch`.

        Parameters:
            data_manager: Data manager containing spectra to be baseline-corrected.
            method: One of {'locmin', 'snip', 'asls'}.
            smooth: LocMin smoothing method; one of {'none', 'loess', 'spline'}.
            span: Loess span proportion in (0, 1]; used when method='locmin'.
            s: Spline smoothing target RSS (>= 0); used when method='locmin'.
            upper: If True, use local maxima as anchors; used when method='locmin'.
            width: Neighborhood width (>= 3) for extrema detection; used when method='locmin'.
            lam: ASLS smoothness parameter; used when method='asls'.
            p: ASLS asymmetry parameter in (0, 1); used when method='asls'.
            niter: ASLS iteration count (> 0); used when method='asls'.
            baseline_scale: Scale factor in (0, 1] applied to the estimated baseline.
            m: SNIP window half-size (>= 1); used when method='snip'.
            decreasing: SNIP decreasing rule; used when method='snip'.
            batch_size: Number of spectra per batch.
            temp_dir: Temporary directory for writing baseline-corrected data.

        Returns:
            MSDataManagerImzML containing baseline-corrected spectra.
        """
        if data_manager is None:
            raise ValueError("data_manager must be provided for baseline correction.")

        method_norm = (method or "asls").strip().lower()
        if method_norm == "asls":
            logger.info(
                f"baseline_correction_entry: method={method}, lam={lam}, p={p}, niter={niter}"
            )
        elif method_norm in {"snip", "snip_numba"}:
            logger.info(
                f"baseline_correction_entry: method={method}, m={m}, decreasing={decreasing}, baseline_scale={baseline_scale}"
            )
        else:
            logger.info(
                f"baseline_correction_entry: method={method}, smooth={smooth}, span={span}, s={s}, upper={upper}, width={width}"
            )

        return Preprocess._process_in_batches(
            data_manager=data_manager,
            batch_size=batch_size,
            temp_dir=temp_dir,
            batch_func=BatchPreprocess.baseline_correction_batch,
            method=method,
            smooth=smooth,
            span=span,
            s=s,
            upper=upper,
            width=width,
            lam=lam,
            p=p,
            niter=niter,
            baseline_scale=baseline_scale,
            m=m,
            decreasing=decreasing,
        )
