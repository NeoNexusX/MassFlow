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
from typing import Union,Optional
import numpy as np
from module.ms_module import SpectrumBaseModule, SpectrumImzML, MS
from logger import get_logger
from .filter import (smooth_signal_ma, smooth_signal_gaussian, smooth_ns_signal_ma,
                     smooth_ns_signal_gaussian, smooth_ns_signal_bi,smooth_signal_savgol,
                     smooth_signal_wavelet,smooth_preprocess)
from .baseline_correction import asls_baseline, snip_baseline
from .est_noise_helper import _findbins,estimation_fun
from scipy.interpolate import InterpolatedUnivariateSpline
from .peak_alignment import peak_matching, get_reference_mz_axis, get_reference_mz_axis2

logger = get_logger("ms_preprocess")

class MSIPreprocessor():
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
    def peak_pick(data:SpectrumBaseModule,method: str) -> SpectrumBaseModule:
        pass

    @staticmethod
    def peak_pick_spectrum(data:SpectrumBaseModule,method: str) -> SpectrumBaseModule:
        pass

    @staticmethod
    def tic_normalization(data:SpectrumBaseModule,method: str = 'total_ion_current') -> SpectrumBaseModule:
        """
        Perform Total Ion Current (TIC) normalization.
        
        Abstract method for normalizing MSI data using various TIC-based approaches
        to account for variations in total ion intensity across pixels.
        
        Args:
            method (str): Normalization method ('total_ion_current', 'median', 'max')
            
        Returns:
            processed_data (MSI): Processed MSI object with peak-picked data
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """

    @staticmethod
    def peak_alignment(
            # reference_mz_axis：'mean' 或 'histogram'
            ref_method: str = 'mean',
            # Computes reference axis parameters (mean spectrum + local extrema) using the mean method.
            agg: str = 'mean',
            round_digits: int = 6,
            half_window: int = 3,
            findlimits: bool = False,
            snr_threshold: Optional[float] = None,
            min_height: Optional[float] = None,
            min_distance_da: Optional[float] = None,
            min_distance_factor: float = 1.0,
            noise_method: str = 'mad',
            # Histogram/Loess reference_axis_parameters
            mz_res: Optional[float] = None,
            px_perc: float = 0.01,
            N_sample: Optional[int] = None,
            smoothing_window: int = 11,
            # parameters for peak alignment (intended for use in 'peak_matching')
            n_sample: int = 2000,
            gap_to_tol_factor: float = 0.5,
            combiner: str = 'max',
            unique: bool = False,
            match_method: str = 'diff',
            # Input and common parameters
            ms_data: MS = None,
            reference_mz_axis: Optional[np.ndarray] = None,
            tolerance: Optional[float] = None,
            units: str = 'ppm',
    ) -> MS:
        """
        Unified peak alignment API: selects reference axis generation method based on ref_method,
        then aligns spectra using diff or dp method.

        Parameters:
        - ref_method: 'mean' for mean spectrum + local extrema; 'histogram' for histogram + Loess.
        - mean parameters: agg, round_digits, half_window, findlimits, snr_threshold, min_height, min_distance_da, min_distance_factor, noise_method.
        - Histogram/Loess reference_axis_parameters: mz_res, px_perc, N_sample, smoothing_window.
        - peak alignment parameters: n_sample, gap_to_tol_factor, combiner('max'/'sum'/'mean'), unique, match_method('diff'/'dp').
        - common parameters: ms_data, reference_mz_axis（若提供则跳过参考轴生成）, tolerance, units('ppm'/'da').

        Returns:
        - MS: MS object with aligned spectra, all m/z unified to reference axis.
        """
        # Parameter validation and method selection
        if ms_data is None:
            raise ValueError("ms_data cannot be empty")
        ref_method = (ref_method or 'mean').lower().strip()
        if ref_method not in ('mean', 'histogram'):
            raise ValueError("`ref_method` must be either 'mean' or 'histogram'")
        match_method = (match_method or 'diff').lower().strip()
        if match_method not in ('diff', 'dp'):
            raise ValueError("`match_method` must be either 'diff' or 'dp'")

        # Reference axis generation: if not provided, generate it based on the specified method.
        if reference_mz_axis is None:
            if ref_method == 'mean':
                reference_mz_axis = get_reference_mz_axis(
                    ms_data=ms_data,
                    agg=agg,
                    round_digits=round_digits,
                    half_window=half_window,
                    findlimits=findlimits,
                    snr_threshold=snr_threshold,
                    noise_method=noise_method,
                    min_height=min_height,
                    n_sample=n_sample,
                    min_distance_da=min_distance_da,
                    min_distance_factor=min_distance_factor,
                )
            else:  # 'histogram'
                reference_mz_axis = get_reference_mz_axis2(
                    ms_data=ms_data,
                    mz_res=mz_res,
                    px_perc=px_perc,
                    N_sample=N_sample,
                    smoothing_window=smoothing_window,
                )

        # Peak alignment: map spectra to reference axis and aggregate using combiner
        ms_aligned = peak_matching(
            ms_data=ms_data,
            reference_mz_axis=reference_mz_axis,
            tolerance=tolerance,
            units=units,
            n_sample=n_sample,
            gap_to_tol_factor=gap_to_tol_factor,
            combiner=combiner,
            unique=unique,
            method=match_method,
        )

        return ms_aligned

    @staticmethod
    def baseline_correction(
        data: Union[np.ndarray, SpectrumBaseModule],
        method: str = "asls",
        lam: float = 1e7,
        p: float = 0.01,
        niter: int = 15,
        baseline_scale: float = 0.8,
        m: int = 5,
        decreasing: bool = True,
        epsilon: float = 1e-3
    ) -> tuple[Union[np.ndarray, SpectrumBaseModule], np.ndarray]:
        """
        Remove baseline drift and background signals from mass spectra.

        Supports:
        - ASLS (asymmetric least squares): robust baseline estimation with peak preservation.
        - SNIP (statistics-sensitive non-linear iterative peak-clipping): progressive clipping with adaptive early-stop.

        Args:
            data: intensity array (np.ndarray) or SpectrumBaseModule
            method: 'asls' (default) or 'snip'
            lam: ASLS smoothness parameter (1e4-1e8, higher = smoother baseline)
            p: ASLS asymmetry parameter (0.001-0.1, lower = more peak preservation)
            niter: ASLS iterations (5-20)
            baseline_scale: scale factor (0-1) applied to the estimated baseline
            m: SNIP max half-window; None -> auto(min(50, n//10))
            decreasing: SNIP iteration order; True: p=m..1 (coarse->fine), False: p=1..m
            epsilon: SNIP adaptive early-stop threshold on relative change

        Returns:
            A tuple containing:
            - np.ndarray or SpectrumBaseModule with baseline-corrected intensity.
            - np.ndarray: The estimated baseline.
        """

        def apply_single(intensity: np.ndarray):
            xi = np.array(intensity, dtype=np.float64, copy=True)
            xi = np.ascontiguousarray(xi)
            # Estimate baseline via module-level functions
            if method == "asls":
                baseline = asls_baseline(xi, lam=lam, p=p, niter=niter)
            elif method == "snip":
                baseline = snip_baseline(xi, m=m, decreasing=decreasing, epsilon=epsilon)
            else:
                raise ValueError("Unsupported baseline method: use 'asls' or 'snip'")
            scale = float(np.clip(baseline_scale, 0.0, 1.0))
            scaled_baseline = scale * baseline
            corrected = xi - scaled_baseline
            corrected = np.maximum(corrected, 0.0)
            return corrected, scaled_baseline

        if isinstance(data, np.ndarray):
            return apply_single(data)
        elif isinstance(data, SpectrumBaseModule):
            corrected_intensity, baseline = apply_single(np.array(data.intensity, dtype=np.float64))
            corrected_spectrum = SpectrumBaseModule(
                mz_list=np.array(data.mz_list, dtype=np.float64, copy=True) if data.mz_list is not None else None,
                intensity=np.ascontiguousarray(corrected_intensity),
                coordinates=data.coordinates
            )
            return corrected_spectrum, baseline
        else:
            raise TypeError("data must be np.ndarray or SpectrumBaseModule")

    @staticmethod
    def noise_reduction(
        data: Union[SpectrumBaseModule, SpectrumImzML],
        method: str = "ma",
        window: int = 5,
        sd: float = None,
        sd_intensity: float = None,
        p: int = 2,
        coef: np.ndarray = None,
        polyorder: int = 2,
        wavelet: str = 'db4',
        threshold_mode: str = 'soft'
    ) -> Union[SpectrumBaseModule, SpectrumImzML]:
        """
        Perform noise reduction on MSI data.
        
        Reduce noise in spectra while preserving important features using
        several algorithms.
        
        Args:
            data: MSBaseModule or MS object containing spectral data
            method (str): Denoising method ('ma', 'gaussian', 'savgol', 'wavelet')
            window (int): Window size for denoising (must be a positive integer)
            sd (float): Standard deviation for Gaussian filter (used when method='gaussian')
            coef (np.ndarray, optional): Custom convolution kernel coefficients for 'ma' method
            polyorder (int): Polynomial order for Savitzky-Golay filter (must be less than window)
            wavelet (str): Wavelet type for wavelet denoising ('db4', 'db8', 'haar', 'coif2', etc.)
            threshold_mode (str): Thresholding mode for wavelet denoising ('soft' or 'hard')

            
        Returns:
            MSBaseModule: New spectrum with the same coordinates and mz_list, but smoothed intensity.
            
        Raises:
            TypeError: If data is not MSBaseModule or MS
            ValueError: If method is not supported
        """
        # Dispatch based on input type without nested helper function
        # Apply smoothing based on method by passing MSBaseModule
        data = smooth_preprocess(data)  # Preprocess to ensure no negative intensities
        if method == "ma":
            smoothed_intensity = smooth_signal_ma(data, coef=coef, window=window)
        elif method == "gaussian":
            smoothed_intensity = smooth_signal_gaussian(data, sd=sd, window=window)
        elif method == "ma_ns":
            smoothed_intensity = smooth_ns_signal_ma(data, p=p, k=window)
        elif method == "savgol":
            smoothed_intensity = smooth_signal_savgol(data, window=window, polyorder=polyorder)
        elif method == "wavelet":
            smoothed_intensity = smooth_signal_wavelet(data, wavelet=wavelet, threshold_mode=threshold_mode)
        elif method == "gaussian_ns":
            smoothed_intensity = smooth_ns_signal_gaussian(data, sd=sd,p=p,k=window)
        elif method == "bi_ns":
            smoothed_intensity = smooth_ns_signal_bi(data, sd_dist=sd, sd_intensity=sd_intensity, p=p,k=window)
        else:
            supported = "ma, gaussian, ma_ns, gaussian_ns, bi_ns"
            logger.error(f"Unsupported smoothing method: {method}. Use one of: {supported}.")
            raise ValueError(f"Unsupported smoothing method: {method}. Use one of: {supported}.")

        return SpectrumBaseModule(
            mz_list=data.mz_list,
            intensity=smoothed_intensity,
            coordinates=data.coordinates,
        )

    @staticmethod
    def preprocess_pipeline(data:SpectrumBaseModule) -> SpectrumBaseModule:
        """
        Execute a complete preprocessing pipeline.
        
        Applies multiple preprocessing steps in sequence to the MSI data.
        
        Args:
            steps (List[str]): List of preprocessing steps to apply
            **kwargs: Additional parameters for preprocessing methods
            
        Returns:
            np.ndarray: Fully preprocessed MSI data
            
        Raises:
            ValueError: If invalid preprocessing step is specified
        """

    @staticmethod
    def est_noise(x: Optional[SpectrumBaseModule] = None,
                nbins: int = 1,
                overlap: float = 0.5,
                k: int = 10,
                method: str = "sd",
                denoise_method: str = "wavelet") -> Union[np.ndarray,float]:
        """
        Estimate noise level in the MSI data.
        
        Args:
            data (MSBaseModule): Input MSI data
            
        Returns:
            np.ndarray: A 1D array of length equal to `x.intensity.size` representing the
                estimated noise baseline. When `nbins <= 1`, the array is constant with the
                global noise estimate value.

        Raises:
            TypeError: If `x.intensity` is not a 1D numpy array.
            ValueError: If `nbins < 1` or `k <= 0`.
        """
        if x.intensity is None or x.intensity.ndim != 1:
            logger.error("x.intensity must be a 1D numpy array")
            raise TypeError("x.intensity must be a 1D numpy array")
        if nbins < 1:
            logger.error("nbins must be >= 1")
            raise ValueError("nbins must be >= 1")
        if k <= 0:
            logger.error("k must be a positive integer")
            raise ValueError("k must be a positive integer")

        # Smooth signal (neighborhood-search Gaussian) and compute absolute residuals
        smoothed = MSIPreprocessor.noise_reduction(x, window=k, method=denoise_method)

        residuals = np.abs(smoothed.intensity - x.intensity)

        #noise_estimation part
        noise_estimation = estimation_fun(method)(residuals)

        #find bins
        if nbins > 1:
            bins, meta = _findbins(residuals, nbins=nbins, overlap=overlap)
            lower = meta["lower"].astype(int)
            upper = meta["upper"].astype(int)

            #core location
            midpoints = (lower + upper) // 2

            #noise_estimation for all bins
            noise_estimations = []
            for bin_data in bins:
                noise_estimations.append(estimation_fun(method)(bin_data))

            #spline for noise line
            rank_spline = int(max(1, min(3, len(midpoints) - 1)))
            spline_fn = InterpolatedUnivariateSpline(midpoints, noise_estimations, k=rank_spline)
            noise_estimation = spline_fn(np.arange(len(residuals), dtype=float))

        return noise_estimation

    @staticmethod
    def calculate_snr(spectrum: SpectrumBaseModule, method="mad") -> float:
        signal_level = np.percentile(spectrum.intensity, 95)
        noise = MSIPreprocessor.est_noise(spectrum, method=method)
        return signal_level / noise
