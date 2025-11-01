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
"""

from typing import Union
from logger import get_logger
import numpy as np
from module.ms_module import MS, MSBaseModule

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
        data: Union[MSBaseModule, MS],
        method: str = "ma",
        window: int = 2,
        sd: float = 2,
        coef: np.ndarray = None,
        polyorder: int = 2,
        wavelet: str = 'db4',
        threshold_mode: str = 'soft'
    ) -> Union[MSBaseModule, MS]:
        """
        Perform noise reduction on MSI data.
        
        Reduces noise in mass spectra while preserving important spectral features
        using various denoising algorithms.
        
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
            Union[MSBaseModule, MS]: Processed MSI object with noise-reduced data
            
        Raises:
            TypeError: If data is not MSBaseModule or MS
            ValueError: If method is not supported
        """

        def smooth_signal_ma(x: np.ndarray, coef=None, window: int = 5) -> np.ndarray:
            """
            Moving average smoothing (or arbitrary convolution kernel smoothing).
            
            Args:
                x : np.ndarray
                    Input signal
                coef : np.ndarray, optional
                    Convolution kernel coefficients (if None, use uniform window)
                window : int
                    Window size (must be a positive integer)
                    
            Returns:
                np.ndarray
                    Smoothed signal
            """
            # If weights not specified, use uniform weights (moving average)
            if coef is None:
                # Ensure window length is odd
                window = window + 1 - window % 2
                coef = np.ones(window)

            # Normalize kernel weights
            coef = coef / np.sum(coef)
            window = len(coef)
            half_window = window // 2

            # Boundary padding: extend using edge values
            xpad = np.pad(x, (half_window, half_window), mode='edge')

            # Convolution filtering
            y = np.convolve(xpad, coef, mode='valid')

            return y

        def smooth_signal_gaussian(x: np.ndarray, sd=None, window=5) -> np.ndarray:
            """
            Gaussian smoothing using NumPy implementation.
            
            Args:
                x : np.ndarray
                    Input signal
                sd : float
                    Standard deviation for Gaussian kernel
                window : int
                    Window size
                    
            Returns:
                np.ndarray
                    Smoothed signal
            """
            from scipy.stats import norm
            # Generate Gaussian kernel
            if sd is None:
                sd = window / 4.0

            half_window = window // 2
            # Generate Gaussian weights
            positions = np.arange(-half_window, half_window + 1)
            coef = norm.pdf(positions, scale=sd)

            return smooth_signal_ma(x, coef=coef)
        
        def smooth_signal_savgol(x: np.ndarray, window: int = 5, polyorder: int = 2) -> np.ndarray:
            """
            Savitzky-Golay filter for signal smoothing.
            
            The Savitzky-Golay filter fits successive sub-sets of adjacent data points 
            with a low-degree polynomial by the method of linear least squares.
            
            Args:
                x : np.ndarray
                    Input signal
                window : int
                    Window size (must be odd and greater than polyorder)
                polyorder : int
                    Polynomial order (must be less than window)
                    
            Returns:
                np.ndarray
                    Smoothed signal
            """
            from scipy.signal import savgol_filter
            # Ensure window is odd
            if window % 2 == 0:
                window += 1
            
            # Ensure polyorder is less than window
            if polyorder >= window:
                polyorder = window - 1
                
            # Minimum window size check
            if window < 3:
                window = 3
                
            # Apply Savitzky-Golay filter
            return savgol_filter(x, window, polyorder)
        
        def smooth_signal_wavelet(x: np.ndarray, wavelet: str = 'db4', threshold_mode: str = 'soft') -> np.ndarray:
            """
            Wavelet denoising for signal smoothing.
            """
            import pywt
            # Ensure writable, contiguous array
            original_length = len(x)
            x = np.array(x, dtype=np.float64, copy=True)
            x = np.ascontiguousarray(x)

            # Perform wavelet decomposition
            coeffs = pywt.wavedec(x, wavelet, mode='symmetric')
            
            # Estimate noise standard deviation using the finest detail coefficients
            sigma = np.median(np.abs(coeffs[-1])) / 0.6745
            
            # Calculate threshold using Donoho-Johnstone threshold
            threshold = sigma * np.sqrt(2 * np.log(len(x)))
            
            # Apply thresholding to all detail coefficients
            coeffs_thresh = list(coeffs)
            coeffs_thresh[1:] = [pywt.threshold(detail, threshold, mode=threshold_mode) 
                                for detail in coeffs[1:]]
            
            # Reconstruct the signal
            reconstructed = pywt.waverec(coeffs_thresh, wavelet, mode='symmetric')

            # Match output length exactly to input
            if len(reconstructed) != original_length:
                if len(reconstructed) > original_length:
                    reconstructed = reconstructed[:original_length]
                else:
                    pad_length = original_length - len(reconstructed)
                    reconstructed = np.pad(reconstructed, (0, pad_length), mode='edge')

            return reconstructed
        
        def _apply_smoothing_single(spectrum: MSBaseModule, method: str, window: int, sd: float, coef, polyorder: int, wavelet: str, threshold_mode: str) -> MSBaseModule:
            """
            Apply smoothing to a single spectrum.
            """
            # Force writable copies to avoid read-only buffer errors
            mz_array = np.array(spectrum.mz_list, dtype=np.float64, copy=True)
            mz_array = np.ascontiguousarray(mz_array)
            intensity_array = np.array(spectrum.intensity, dtype=np.float64, copy=True)
            intensity_array = np.ascontiguousarray(intensity_array)
            
            # Apply smoothing based on method
            if method == "ma":
                smoothed_intensity = smooth_signal_ma(intensity_array, coef=coef, window=window)
            elif method == "gaussian":
                smoothed_intensity = smooth_signal_gaussian(intensity_array, sd=sd, window=window)
            elif method == "savgol":
                smoothed_intensity = smooth_signal_savgol(intensity_array, window=window, polyorder=polyorder)
            elif method == "wavelet":
                smoothed_intensity = smooth_signal_wavelet(intensity_array, wavelet=wavelet, threshold_mode=threshold_mode)
            else:
                raise ValueError(f"Unsupported smoothing method: {method}. Use 'ma', 'gaussian', 'savgol', or 'wavelet'.")
            
            # Ensure smoothed intensity has the same length as mz_array
            if len(smoothed_intensity) != len(mz_array):
                print(f"Warning: {method} produced length mismatch - original: {len(mz_array)}, processed: {len(smoothed_intensity)}")
                min_length = min(len(mz_array), len(smoothed_intensity))
                if len(smoothed_intensity) > len(mz_array):
                    print(f"  -> Truncated processed data to original length: {len(mz_array)}")
                    smoothed_intensity = smoothed_intensity[:len(mz_array)]
                else:
                    print(f"  -> Truncated both m/z and intensity to minimum length: {min_length}")
                    mz_array = mz_array[:min_length]
                    smoothed_intensity = smoothed_intensity[:min_length]
            
            # Create new spectrum with smoothed data (contiguous)
            smoothed_intensity = np.ascontiguousarray(smoothed_intensity)
            smoothed_spectrum = MSBaseModule(
                mz_list=mz_array,
                intensity=smoothed_intensity,
                coordinates=spectrum.coordinates
            )
            return smoothed_spectrum

        # Dispatch based on input type
        if isinstance(data, MSBaseModule):
            return _apply_smoothing_single(data, method, window, sd, coef, polyorder, wavelet, threshold_mode)
        elif isinstance(data, MS):
            # Build a new MS with smoothed spectra
            new_ms = MS()
            for s in data:
                smoothed_s = _apply_smoothing_single(s, method, window, sd, coef, polyorder, wavelet, threshold_mode)
                coords = smoothed_s.get_coordinates()
                x, y, z = coords[0], coords[1], coords[2] if len(coords) > 2 else 0
                new_ms[x, y, z] = smoothed_s
            return new_ms
        else:
            raise TypeError("data must be MSBaseModule or MS")



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
