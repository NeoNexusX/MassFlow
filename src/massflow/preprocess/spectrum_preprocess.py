"""
Mass Spectrometry Imaging (MSI) Preprocessing Module

This module defines an abstract base class for preprocessing Mass Spectrometry Imaging (MSI) data.
It provides a framework for implementing essential preprocessing techniques such as:

- Peak picking: Identifying significant peaks in mass spectra.
- TIC normalization: Normalizing data to account for variations in total ion intensity.
- Peak alignment: Aligning peaks across spectra to correct for mass calibration drift.
- Baseline correction: Removing baseline drift and background signals.
- Noise reduction: Reducing noise while preserving spectral features.

The `MSIPreprocessor` class is designed to be extended by concrete implementations that define
specific algorithms for each preprocessing step. It supports both traditional input formats
(e.g., mz, msroi arrays) and MSI object integration for seamless compatibility with the MSI framework.

Classes:
- MSIPreprocessor: Abstract base class for defining preprocessing workflows.

Methods:
- peak_pick: Abstract method for peak picking.
- peak_pick_pixel: Abstract method for pixel-level peak picking.
- tic_normalization: Abstract method for TIC-based normalization.
- peak_alignment: Abstract method for aligning peaks across spectra.
- baseline_correction: Abstract method for correcting baseline drift.
- noise_reduction: Abstract method for reducing noise in mass spectra.
- preprocess_pipeline: Abstract method for executing a complete preprocessing pipeline.

This module is intended for developers extending the MSI framework with custom preprocessing
algorithms.

Author: MassFlow Development Team Bionet/NeoNexus lyk
License: See LICENSE file in project root
"""
from typing import Union, Optional, Sequence
import numpy as np
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.spectrum import Spectrum
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.logger import get_logger
from massflow.preprocess.peak_align_helper import align_massdata, align_spectrum
from massflow.r_preprocess.adapter import CardinalAdapter
from massflow.preprocess.filter_helper import smoother
from massflow.preprocess.normalizer_helper import normalizer
from massflow.preprocess.baseline_correction_helper import baseline_corrector
from massflow.preprocess.est_noise_helper import estimator
from massflow.preprocess.peak_pick_helper import peak_pick_fun
import os
import time

logger = get_logger("ms_preprocess")

class SpectrumPreprocess:
    """
    Abstract base class for MSI data preprocessing.

    This class provides a framework for implementing various preprocessing techniques
    for Mass Spectrometry Imaging (MSI) data, including peak picking, normalization,
    alignment, baseline correction, and noise reduction.

    Supports both traditional input (mz, msroi arrays) and MSI object input for
    better integration with the MSI framework.

    Attributes:
        msi_object (MSI, optional): MSI object containing data and metadata
        preprocessing_params (dict): Parameters for preprocessing operations
        processed_data (np.ndarray, optional): Processed MSI data
    """

    def __init__(self):
        """
        Initialize MSI preprocessor.

        Args:
            msi_object (MSI, optional): MSI object containing data and metadata
            preprocessing_params (dict, optional): Parameters for preprocessing operations

        Raises:
            ValueError: If neither msi_object nor (mz, msroi) are provided
        """


    @staticmethod
    def base_input_check(data):
        """
        Validate input data.

        Raises:
            ValueError: If input data is invalid or missing required attributes
        """
        if data.intensity is None or data.mz_list is None:
            raise ValueError("Input spectrum must have both intensity and mz_list data.")

    @staticmethod
    def peak_pick_spectrum( data: Spectrum,
                            width: int | Sequence[int] = 2,
                            method: str = 'scipy',
                            relheight: float = 0.1,
                            return_type: str = 'height') -> Spectrum:
        """
        Perform peak picking on a single spectrum and return a reduced spectrum.

        Parameters:
            data (SpectrumBaseModule): Input spectrum.
            width (int): Required peak width for underlying detector.
            method (str): Peak pick backend; currently supports 'scipy'.
            relheight (float): Relative height threshold for candidate peaks.
            return_type (str): 'height' for peak heights or 'area' for integrated areas.

        Returns:
            SpectrumBaseModule: New spectrum containing picked peaks and original coordinates.

        Raises:
            ValueError: If `method` or `return_type` is unsupported.
        """
        SpectrumPreprocess.base_input_check(data=data)

        intensity = data.intensity
        index = data.mz_list
        peak_intensity,peak_index = peak_pick_fun(  intensity, # type: ignore
                                                    index, # type: ignore
                                                    width=width,
                                                    method=method,
                                                    relheight=relheight,
                                                    return_type=return_type)

        return Spectrum(
            mz_list=peak_index,
            intensity=peak_intensity,
            coordinate=data.coordinate,
        )

    @staticmethod
    def normalization_spectrum(
        data: Union[Spectrum, SpectrumImzML],
        scale_method: str = "none",
        method: str = "tic",
        scale: float = 1.0
    ) -> Union[Spectrum, SpectrumImzML]:
        """
        Normalize a single spectrum using TIC, RMS, or Median, with optional scaling.

        Parameters:
            data (SpectrumBaseModule | SpectrumImzML): Spectrum to normalize.
            scale_method (str): 'none' or 'unit' min-max scaling.
            method (str): One of {'tic', 'rms', 'median'}.
            scale (float): Cardinal-like amplitude scaling factor applied after normalization.

        Returns:
            SpectrumBaseModule | SpectrumImzML: Spectrum with normalized intensity.

        Raises:
            ValueError: If `method` is unsupported or if TIC/RMS/Median is invalid (≤ 0) in helper functions.
        """
        intensity = data.intensity
        norm_intensity = normalizer(
            intensity, # type: ignore
            intensity, # type: ignore
            scale_method=scale_method,
            method=method,
            scale=scale
        )

        return Spectrum(
            mz_list=data.mz_list,
            intensity=norm_intensity,
            coordinate=data.coordinate,
        )

    @staticmethod
    def peak_alignment(data_manager: Optional[MSDataManagerImzML] = None,
                       spectrum: Optional[SpectrumImzML] = None,
                       ref: Optional[np.ndarray] = None,
                       units: str = 'ppm',
                       tolerance: Optional[float] = None,
                       binfun: str = 'median',
                       binratio: int = 2,
                       backend_method: Optional[str] = "python",
                       batch_size: int = 256,
                       clear_memory: bool = True,
                       temp_dir: str = "./temp_align_data"
                       ):
        """
        Align peaks across spectra in MS data using specified backend.

        This method provides a unified interface for peak alignment, supporting both
        Python-based implementation and R-based Cardinal implementation. It can align
        an entire dataset (MSDataManagerImzML) or a single spectrum (SpectrumImzML).

        Parameters:
            data_manager (MSDataManagerImzML, optional): The data manager containing the mass spectra to align.
                Required if `spectrum` is None.
            spectrum (SpectrumImzML, optional): A single spectrum to align.
                Required if `data_manager` is None.
            ref (np.ndarray, optional): External reference m/z axis.
                If None, it will be estimated from the data (for data_manager) or must be provided (for single spectrum).
            units (str): Units for tolerance and resolution ('ppm' or 'absolute'). Default is 'ppm'.
            tolerance (float, optional): The tolerance window for peak matching.
                If None, it will be estimated from the data.
            binfun (str): Aggregation function for estimating resolution ('median', 'min', 'max', 'mean').
                Default is 'median'.
            binratio (int): Ratio to scale the estimated resolution to determine tolerance. Default is 2.
            backend_method (str, optional): The backend to use for alignment.
                - 'cardinal': Use the R Cardinal package (requires R environment).
                - 'python' (or None): Use the native Python implementation.

        Returns:
            MSDataManagerImzML | SpectrumImzML: The aligned data manager or spectrum.

        Raises:
            ValueError: If neither `data_manager` nor `spectrum` is provided, or if required parameters
                (ref, tolerance) are missing for single spectrum alignment.
        """

        logger.info(
            "peak_alignment_entry: backend=%s, binfun=%s, tolerance=%s, units=%s",
            backend_method,
            binfun,
            tolerance,
            units,
        )

        if backend_method == "cardinal" and data_manager is not None:
            aligned_data_manager = CardinalAdapter.align_massdata(
                                                                  dm_data=data_manager,
                                                                  reference=ref,
                                                                  tolerance=tolerance,
                                                                  units=units,
                                                                  binfun=binfun,
                                                                  binratio=binratio,
                                                                  temp_dir=temp_dir)
            return aligned_data_manager

        tolerance = tolerance * 1e-6 if tolerance is not None and units == "ppm" else tolerance

        if spectrum is not None:
            if ref is None or tolerance is None:
                raise ValueError("Reference and tolerance must be provided for single spectrum alignment.")
            aligned_spectrum = align_spectrum(spectrum=spectrum,
                                              reference=ref,
                                              tolerance=tolerance,
                                              units=units)
            return aligned_spectrum

        if data_manager is not None and spectrum is None:
            aligned_data_manager = align_massdata(data_manager=data_manager,
                                                  reference=ref,
                                                  tolerance=tolerance,
                                                  units=units,
                                                  binfun=binfun,
                                                  binratio=binratio,
                                                  batch_size=batch_size,
                                                  clear_memory=clear_memory,
                                                  temp_dir=temp_dir)
            return aligned_data_manager

        raise ValueError("Either 'data_manager' or 'spectrum' must be provided for alignment.")

    @staticmethod
    def baseline_correction_spectrum(
        data: Union[Spectrum, SpectrumImzML],
        method: str = "asls",
        smooth: str = "none",
        span: float = 0.1,
        s: Optional[float] = 0.0,
        upper: bool = False,
        width: int = 5,
        lam: float = 1e7,
        p: float = 0.01,
        niter: int = 15,
        baseline_scale: float = 1.0,
        m: Optional[int] = None,
        decreasing: bool = True,
    ) -> tuple[Spectrum, np.ndarray]:
        """
        Baseline correction using LocMin, SNIP, or ASLS with optional baseline scaling.

        Parameters:
            data (SpectrumBaseModule | SpectrumImzML): Input spectrum containing 1D `intensity` and `mz_list`.
                The corrected spectrum preserves `mz_list` and `coordinates`.
            method (str): One of {'locmin', 'snip', 'asls'}.
            smooth (str): LocMin smoothing method; one of {'none', 'loess', 'spline'}. Used when `method='locmin'`.
            span (float): Loess span proportion in (0, 1]; typical 0.05–0.2. Used when `method='locmin'`.
            s (float, optional): Spline smoothing target residual sum of squares (>= 0); `0.0` is interpolation. Used when `method='locmin'`.
            upper (bool): If True, use local maxima as anchors (upper envelope); otherwise minima. Used when `method='locmin'`.
            width (int): Neighborhood width (>= 3) for local extrema detection. Used when `method='locmin'`.
            lam (float): ASLS smoothness parameter. Used when `method='asls'`.
            p (float): ASLS asymmetry parameter in (0, 1). Used when `method='asls'`.
            niter (int): ASLS iteration count (> 0). Used when `method='asls'`.
            baseline_scale (float): Scale factor in (0, 1] applied to the estimated baseline prior to subtraction.
            m (int, optional): SNIP window half-size (>= 1). Used when `method='snip'`.
            decreasing (bool): SNIP decreasing rule; iterate from large window to small when True.

        Returns:
            Tuple[SpectrumBaseModule, np.ndarray]:
                - corrected spectrum (retains original `mz_list` and `coordinates`)
                - scaled baseline vector

        Raises:
            ValueError: Unsupported `method`.
            TypeError: Invalid input type or non-1D arrays (propagated from helpers).
            ValueError/ImportError: Smoothing dependency/parameter errors when `smooth in {'loess','spline'}` (propagated).

        Notes:
            - `index` is passed to the helper for validation only; current algorithms do not directly use it.
            - Loess smoothing uses `preprocess.peak_alignment._smooth1d` and spline smoothing uses
              `scipy.interpolate.UnivariateSpline`; ensure dependencies are installed and parameters valid.
        """
        
        intensity = data.intensity
        index = data.mz_list
        
        corrected_intensity, baseline = baseline_corrector(
            intensity, # type: ignore
            intensity, # type: ignore
            index=index,
            method=method,
            lam=lam,
            p=p,
            niter=niter,
            baseline_scale=baseline_scale,
            m=m,
            decreasing=decreasing,
            smooth=smooth,
            span=span,
            s=s,
            upper=upper,
            width=width,
        )
        corrected_spectrum = Spectrum(
            mz_list=data.mz_list,
            intensity=corrected_intensity,
            coordinate=data.coordinate,
        )
        
        return corrected_spectrum, baseline

    @staticmethod
    def noise_reduction_spectrum(
        data: Union[Spectrum, SpectrumImzML],
        method: str = "ma",
        window: int = 5,
        sd: Optional[float] = None,
        sd_intensity: Optional[float] = None,
        p: int = 2,
        coef: Optional[np.ndarray] = None,
        polyorder: int = 3,
        derive: int = 0,
        delta: float = 1.0,
        wavelet: str = "db4",
        threshold_mode: str = "soft",
    ) -> Union[Spectrum, SpectrumImzML]:
        """
        Reduce spectral noise while preserving features using multiple algorithms.

        Parameters:
            data (SpectrumBaseModule | SpectrumImzML): Spectrum to denoise.
            method (str): One of {'ma','gaussian','savgol','wavelet','ma_ns','gaussian_ns','bi_ns'}.
            window (int): Window size or neighbor count depending on method.
            sd (float, optional): Gaussian scale parameter.
            coef (np.ndarray, optional): Custom kernel for 'ma'.
            polyorder (int): Polynomial order for Savitzky-Golay.
            wavelet (str): Wavelet family for wavelet denoising.
            threshold_mode (str): 'soft' or 'hard' thresholding.
            sd_intensity (float, optional): Intensity scale for bilateral method.
            p (int): Minkowski metric for NS queries.

        Returns:
            SpectrumBaseModule: New spectrum with smoothed intensity and original coordinates.

        Raises:
            ValueError: If `method` is unsupported.
        """
        SpectrumPreprocess.base_input_check(data=data)

        intensity = data.intensity
        index = data.mz_list

        smoothed_intensity = smoother(
            intensity,# type: ignore
            index=index,
            method=method,
            window=window,
            sd=sd,# type: ignore
            sd_intensity=sd_intensity,# type: ignore
            p=p,
            coef=coef,# type: ignore
            polyorder=polyorder,
            derive=derive,
            delta=delta,
            wavelet=wavelet,
            threshold_mode=threshold_mode,
        )

        return Spectrum(
            mz_list=data.mz_list,
            intensity=smoothed_intensity,
            coordinate=data.coordinate,
        )

    @staticmethod
    def noise_estimation_spectrum(
        x: Spectrum,
        nbins: int = 1,
        overlap: float = 0.2,
        method: str = "sd",
        denoise_method: str = "gaussian",
        dynamic: bool = False,
    ):
        """
        Estimate noise level for a spectrum using binning and SSE metrics.

        Parameters:
            x (SpectrumBaseModule): Input spectrum.
            nbins (int): Initial number of bins for segmentation.
            overlap (float): Overlap ratio between adjacent bins.
            method (str): Estimator identifier (e.g., 'sd').
            denoise_method (str): Denoising method for pre-estimation.
            dynamic (bool): Whether to adapt bins dynamically.

        Returns:
            float | np.ndarray: Estimated noise scalar or per-bin array depending on method.
        """
        SpectrumPreprocess.base_input_check(data=x)

        intensity = x.intensity
        index = x.mz_list

        return estimator(
            intensity,# type: ignore
            index,# type: ignore
            nbins=nbins,
            overlap=overlap,
            method=method,
            dynamic=dynamic,
            denoise_method=denoise_method,
        )

    @staticmethod
    def calculate_snr_spectrum(
        spectrum: Spectrum, method="sd"
    ) -> float:
        """
        Calculate the signal-to-noise ratio (SNR) for a spectrum.

        Parameters:
            spectrum (SpectrumBaseModule): Input spectrum.
            method (str): Noise estimation method forwarded to `noise_estimation_spectrum`.

        Returns:
            float: SNR computed as quantile-based signal level divided by estimated noise.
        """
        SpectrumPreprocess.base_input_check(data=spectrum)

        signal_level = np.percentile(spectrum.intensity, 95)# type: ignore

        noise = SpectrumPreprocess.noise_estimation_spectrum(spectrum, method=method)

        logger.info(f"SNR: signal_level:{signal_level}, noise:{noise}")
        return signal_level / noise


