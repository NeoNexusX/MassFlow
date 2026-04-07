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
from massflow.module import SpectrumImzML, Spectrum
from massflow.tools.logger import get_logger
from massflow.preprocess.helper.peak_align_helper import align_spectrum
from massflow.preprocess.helper.noise_reduction_helper import smoother
from massflow.preprocess.helper.normalizer_helper import normalizer
from massflow.preprocess.helper.baseline_correction_helper import baseline_corrector
from massflow.preprocess.helper.est_noise_helper import estimator
from massflow.preprocess.helper.peak_pick_helper import peak_pick_fun

logger = get_logger("massflow.preprocess")

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
                            method: str = 'origin',
                            relheight: float = 0.012,
                            snr: float = 2.0,
                            return_type: str = 'height',
                            use_numba: bool = True
                            ) -> SpectrumImzML:
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
        peak_intensity,peak_index = peak_pick_fun(intensity, # type: ignore
                                                  index, # type: ignore
                                                  width=width,
                                                  method=method,
                                                  relheight=relheight,
                                                  snr=snr,
                                                  return_type=return_type,
                                                  use_numba=use_numba
                                                  )

        return SpectrumImzML(
            mz_list=peak_index,
            intensity=peak_intensity,
            coordinates=data.coordinate,
        )

    @staticmethod
    def normalization_spectrum(
        data: Union[Spectrum, SpectrumImzML],
        scale_method: str = "none",
        method: str = "tic",
        scale: float = 1.0,
    ) -> SpectrumImzML:
        """Normalize a single spectrum using TIC, RMS, or Median, with optional scaling.

        Parameters:
            data (Spectrum | SpectrumImzML): Spectrum to normalize.
            scale_method (str): 'none' or 'unit' min-max scaling.
            method (str): One of {'tic', 'rms', 'median'}.
            scale (float): Cardinal-like amplitude scaling factor applied after normalization.

        Returns:
            SpectrumImzML: Spectrum with normalized intensity and original coordinates.

        Raises:
            ValueError: If `method` is unsupported or if TIC/RMS/Median is invalid (≤ 0) in helper functions.
        """
        SpectrumPreprocess.base_input_check(data=data)

        intensity = data.intensity
        norm_intensity = normalizer(
            intensity,  # type: ignore
            scale_method=scale_method,
            method=method,
            scale=scale,
        )

        return SpectrumImzML(
            mz_list=data.mz_list,
            intensity=norm_intensity,
            coordinates=data.coordinate,
        )

    @staticmethod
    def peak_align_spectrum(spectrum: Spectrum,
                            reference: np.ndarray,
                            tolerance: float,
                            units: str = 'ppm',
                            ) -> SpectrumImzML:
        """
        Align peaks across spectra in MS data.
        This method can align a single spectrum (Spectrum).

        Parameters:
            spectrum : A single spectrum to align.
            reference (np.ndarray): External reference m/z axis.
                must be provided (for single spectrum).
            units (str): Units for tolerance and resolution ('ppm' or 'mz'). Default is 'ppm'.
            tolerance (float): The tolerance window for peak matching.
                must be provided (for single spectrum).

        Returns:
            SpectrumImzML: The aligned spectrum.

        Raises:
            ValueError: If required parameters (reference, tolerance) are missing for single spectrum alignment.
        """

        if spectrum is None or reference is None or tolerance is None:
            raise ValueError("For single spectrum alignment, 'spectrum', 'reference', and 'tolerance' must be provided.")

        tolerance = tolerance * 1e-6 if units == "ppm" else tolerance

        return align_spectrum(spectrum=spectrum,
                              reference=reference,
                              tolerance=tolerance,
                              units=units)

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
    ) -> SpectrumImzML:
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
            SpectrumImzML: Corrected spectrum (retains original `mz_list` and `coordinates`).

        Raises:
            ValueError: Unsupported `method`.
            TypeError: Invalid input type or non-1D arrays (propagated from helpers).
            ValueError/ImportError: Smoothing dependency/parameter errors when `smooth in {'loess','spline'}` (propagated).

        Notes:
            - `index` is passed to the helper for validation only; current algorithms do not directly use it.
            - Loess smoothing uses `preprocess.peak_alignment._smooth1d` and spline smoothing uses
              `scipy.interpolate.UnivariateSpline`; ensure dependencies are installed and parameters valid.
        """

        SpectrumPreprocess.base_input_check(data=data)

        intensity = data.intensity
        index = data.mz_list

        if intensity is None or index is None:
            raise ValueError("Input spectrum must have both intensity and mz_list data.")

        corrected_intensity, _ = baseline_corrector(
            intensity,
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
        corrected_spectrum = SpectrumImzML(
            mz_list=data.mz_list,
            intensity=corrected_intensity,
            coordinates=data.coordinate,
        )

        return corrected_spectrum

    @staticmethod
    def noise_reduction_spectrum(
        data: Spectrum,
        method: str = "ma",
        window: int = 5,
        sd: Optional[float] = None,
        sd_intensity: Optional[float] = None,
        p: int = 2,
        coef: Optional[np.ndarray] = None,
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0,
        wavelet: str = "db4",
        threshold_mode: str = "soft",
    ) -> SpectrumImzML:
        """Reduce noise for a single Spectrum while preserving spectral features.

        Parameters:
            data (SpectrumBaseModule): Spectrum to denoise.
            method (str): One of {'ma','ma_numba','gaussian','gaussian_numba','savgol','savgol_numba',
            'wavelet','gaussian_ns','gaussian_ns_numba','bi_ns','bi_ns_numba'}.
            window (int): Window size or neighbor count depending on method.
            sd (float, optional): Gaussian scale parameter.
            coef (np.ndarray, optional): Custom kernel for 'ma'.
            polyorder (int): Polynomial order for Savitzky-Golay.
            deriv (int): The order of the derivative to compute. This must be a nonnegative integer. 
            The default is 0, which means to filter the data without differentiating.
            delta (float): The spacing of the samples to which the filter will be applied.
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
            coef=coef,  # type: ignore
            polyorder=polyorder,
            deriv=deriv,
            delta=delta,
            wavelet=wavelet,
            threshold_mode=threshold_mode,
        )

        return SpectrumImzML(
            mz_list=data.mz_list,
            intensity=smoothed_intensity,
            coordinates=data.coordinate,
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
    def calculate_snr_spectrum(spectrum: Spectrum, method="sd") -> float:
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
