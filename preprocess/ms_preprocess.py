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
from module.ms_module import SpectrumBaseModule, SpectrumImzML
from logger import get_logger
from .filter_helper import (
    smoother,
)
from .baseline_correction_helper import asls_baseline, snip_baseline
from .est_noise_helper import estimator
from .peak_pick_helper import peak_pick_fun

logger = get_logger("ms_preprocess")


class MSIPreprocessor:
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
    def peak_pick_spectrum(data: SpectrumBaseModule,
                            width: int = 2,
                            method: str = 'scipy',
                            relheight: float = 0.1,
                            return_type: str = 'height') -> SpectrumBaseModule:
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

        intensity = data.intensity
        index = data.mz_list
        peak_intensity,peak_index = peak_pick_fun(intensity,
                                                    index,
                                                    width=width,
                                                    method=method,
                                                    relheight=relheight,
                                                    return_type=return_type)

        return SpectrumBaseModule(
            mz_list=peak_index,
            intensity=peak_intensity,
            coordinates=data.coordinates,
        )

    @staticmethod
    def normalization_spectrum(
        data: SpectrumBaseModule, method: str = "total_ion_current"
    ) -> SpectrumBaseModule:
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
    def peak_alignment_spectrum(data: SpectrumBaseModule) -> SpectrumBaseModule:
        """
        Align peaks across spectra to compensate for mass calibration drift.

        Notes:
            This is an abstract placeholder in the preprocessing interface.
            Concrete implementations should provide spectrum-wise alignment.

        Parameters:
            data (SpectrumBaseModule): Input spectrum.

        Returns:
            SpectrumBaseModule: Aligned spectrum (implementation-defined).

        Raises:
            NotImplementedError: If not implemented by subclass.
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
        epsilon: float = 1e-3,
    ) -> tuple[Union[np.ndarray, SpectrumBaseModule], np.ndarray]:
        """
        Baseline correction using ASLS or SNIP with optional baseline scaling.

        Parameters:
            data (Union[np.ndarray, SpectrumBaseModule]): Input 1D intensity or spectrum.
            method (str): 'asls' or 'snip'.
            lam (float): ASLS smoothness parameter.
            p (float): ASLS asymmetry parameter.
            niter (int): ASLS iteration count.
            baseline_scale (float): Scale factor in [0,1] applied to estimated baseline.
            m (int): SNIP window half-size.
            decreasing (bool): SNIP decreasing rule.
            epsilon (float): SNIP early-stop tolerance.

        Returns:
            Tuple[Union[np.ndarray, SpectrumBaseModule], np.ndarray]:
                - corrected data (vector or spectrum)
                - estimated baseline (vector)

        Raises:
            ValueError: If `method` is unsupported.
            TypeError: If `data` is neither `np.ndarray` nor `SpectrumBaseModule`.
        """
        #模仿noise_reduction_spectrum仿写，实现分层架构，
        #删除apply_single函数，重写一个baseline_corrector实现apply_single函数的效果，放在baseline_correction_helper里面
        #输入检查模仿smoother的_input_validation函数，集中检查，不用在这里检查
        #这里仅作为spectrum层级的api来处理输入的data对象

        def apply_single(intensity: np.ndarray):
            xi = np.array(intensity, dtype=np.float64, copy=True)
            xi = np.ascontiguousarray(xi)
            # Estimate baseline via module-level functions
            if method == "asls":
                baseline = asls_baseline(xi, lam=lam, p=p, niter=niter)
            elif method == "snip":
                baseline = snip_baseline(
                    xi, m=m, decreasing=decreasing, epsilon=epsilon
                )
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
            corrected_intensity, baseline = apply_single(
                np.array(data.intensity, dtype=np.float64)
            )
            corrected_spectrum = SpectrumBaseModule(
                mz_list=(
                    np.array(data.mz_list, dtype=np.float64, copy=True)
                    if data.mz_list is not None
                    else None
                ),
                intensity=np.ascontiguousarray(corrected_intensity),
                coordinates=data.coordinates,
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
        wavelet: str = "db4",
        threshold_mode: str = "soft",
    ) -> Union[SpectrumBaseModule, SpectrumImzML]:
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
        intensity = data.intensity
        index = data.mz_list

        smoothed_intensity = smoother(
            intensity,
            index=index,
            method=method,
            window=window,
            sd=sd,
            sd_intensity=sd_intensity,
            p=p,
            coef=coef,
            polyorder=polyorder,
            wavelet=wavelet,
            threshold_mode=threshold_mode,
        )

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
        intensity = x.intensity
        index = x.mz_list

        return estimator(
            intensity,
            index,
            nbins=nbins,
            overlap=overlap,
            method=method,
            dynamic=dynamic,
            denoise_method=denoise_method,
        )

    @staticmethod
    def calculate_snr_spectrum(
        spectrum: SpectrumBaseModule, method="sd"
    ) -> float:
        """
        Calculate the signal-to-noise ratio (SNR) for a spectrum.

        Parameters:
            spectrum (SpectrumBaseModule): Input spectrum.
            method (str): Noise estimation method forwarded to `noise_estimation_spectrum`.

        Returns:
            float: SNR computed as quantile-based signal level divided by estimated noise.
        """
        signal_level = np.percentile(spectrum.intensity, 95)

        noise = MSIPreprocessor.noise_estimation_spectrum(spectrum, method=method)

        logger.info(f"SNR: signal_level:{signal_level}, noise:{noise}")
        return signal_level / noise

    @staticmethod
    def preprocess_pipeline(data: SpectrumBaseModule) -> SpectrumBaseModule:
        """
        Composite preprocessing pipeline (placeholder).

        Parameters:
            data (SpectrumBaseModule): Input spectrum to preprocess.

        Returns:
            SpectrumBaseModule: Preprocessed spectrum (implementation-defined).
        """
