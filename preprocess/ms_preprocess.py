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
from math import log
from typing import Union, Optional
import numpy as np
import inspect
from module.ms_module import SpectrumBaseModule, SpectrumImzML
from logger import get_logger
from .filter_helper import (
    smoother,
)
from .baseline_correction_helper import asls_baseline, snip_baseline
from .est_noise_helper import estimator

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
    def peak_pick_spectrum(data:SpectrumBaseModule,method: str) -> SpectrumBaseModule:
        pass

    @staticmethod
    def normalization_spectrum(data:SpectrumBaseModule,method: str = 'total_ion_current') -> SpectrumBaseModule:
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
    def peak_alignment_spectrum(data:SpectrumBaseModule) -> SpectrumBaseModule:
        """
        Perform peak alignment across spectra.
        
        Abstract method for aligning peaks across different spectra to correct
        for mass calibration drift and improve data consistency.
        
        Args:
            fill this
            
        Returns:
            processed_data (MSI): Processed MSI object with peak_aligned data
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """
        pass

    @staticmethod
    def baseline_correction_spectrum(
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
    def noise_reduction_spectrum(
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
        intensity = data.intensity
        index = data.mz_list

        smoothed_intensity = smoother(intensity, 
                                      index=index, 
                                      method=method, 
                                      window=window, 
                                      sd=sd,
                                      sd_intensity=sd_intensity,
                                      p=p,
                                      coef=coef,
                                      polyorder=polyorder,
                                      wavelet=wavelet,
                                      threshold_mode=threshold_mode)

        return SpectrumBaseModule(
            mz_list=data.mz_list,
            intensity=smoothed_intensity,
            coordinates=data.coordinates,
        )

    @staticmethod
    def noise_estimation_spectrum(
        x: SpectrumBaseModule,
        nbins: int = 1,
        overlap: float = 0.2,
        method: str = "sd",
        denoise_method: str = "bi_ns",
        dynamic: bool = False
    ):
        """
        Estimate noise level in the MSI data.
        """
        intensity = x.intensity
        index = x.mz_list

        return estimator(intensity,
                         index,
                         nbins=nbins,
                         overlap=overlap,
                         method=method,
                         dynamic=dynamic,
                         denoise_method=denoise_method)


    @staticmethod
    def calculate_snr_spectrum(spectrum: SpectrumBaseModule, 
                               method="quantile") -> float:
        """
        pass
        """
        signal_level = np.percentile(spectrum.intensity, 95)

        noise = MSIPreprocessor.noise_estimation_spectrum(spectrum, 
                                                          method=method)
        
        logger.info(f"SNR: signal_level:{signal_level}, noise:{noise}")
        return signal_level / noise

    @staticmethod
    def preprocess_pipeline(data:SpectrumBaseModule) -> SpectrumBaseModule:
        """
        pass
        """