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
from typing import Union
import numpy as np
from module.ms_module import MSBaseModule, MSImzMLBase
from logger import get_logger
from .filter import (smooth_signal_ma, smooth_signal_gaussian, smooth_ns_signal_ma,
                     smooth_ns_signal_gaussian, smooth_ns_signal_bi,smooth_signal_savgol,
                     smooth_signal_wavelet,smooth_preprocess)

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
    def peak_pick(data:MSBaseModule,method: str) -> MSBaseModule:
        """
        Perform peak picking on MSI data.
        
        Abstract method that must be implemented by concrete classes to identify
        significant peaks in the mass spectra.
        
        Args:
            Fill this
            
        Returns:
            processed_data (MSI): Processed MSI object with peak-picked data
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """


    @staticmethod
    def peak_pick_spectrum(data:MSBaseModule,method: str) -> MSBaseModule:
        """
        Perform peak picking on MSI data.

        Abstract method that must be implemented by concrete classes to identify
        significant peaks in the mass spectra.

        Args:
            Fill this

        Returns:
            processed_data (MSI): Processed MSI object with peak-picked data

        Raises:
            NotImplementedError: If not implemented by subclass
        """

    @staticmethod
    def tic_normalization(data:MSBaseModule,method: str = 'total_ion_current') -> MSBaseModule:
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
    def peak_alignment(data:MSBaseModule) -> MSBaseModule:
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

    @staticmethod
    def baseline_correction(data:MSBaseModule) -> MSBaseModule:
        """
        Perform baseline correction on MSI data.
        
        Abstract method for removing baseline drift and background signals
        from mass spectra to improve peak detection and quantification.
        
        Args:
            fill this
            
        Returns:
            processed_data (MSI): Processed MSI object with baseline_correction
            
        Raises:
            NotImplementedError: If not implemented by subclass
        """

    @staticmethod
    def noise_reduction(
        data: Union[MSBaseModule, MSImzMLBase],
        method: str = "ma",
        window: int = 2,
        sd: float = 2,
        sd_intensity: float = None,
        p: int = 2,
        coef: np.ndarray = None,
        polyorder: int = 2,
        wavelet: str = 'db4',
        threshold_mode: str = 'soft'
    ) -> Union[MSBaseModule, MSImzMLBase]:
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

        return MSBaseModule(
            mz_list=data.mz_list,
            intensity=smoothed_intensity,
            coordinates=data.coordinates,
        )

    @staticmethod
    def preprocess_pipeline(data:MSBaseModule) -> MSBaseModule:
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
    def estnoise(data:MSBaseModule) -> float:
        """
        Estimate noise level in the MSI data.
        
        Args:
            data (MSBaseModule): Input MSI data
            
        Returns:
            float: Estimated noise level
        """