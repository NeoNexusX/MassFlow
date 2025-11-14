from typing import Optional
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
from module.ms_module import SpectrumBaseModule, SpectrumImzML, MS
from logger import get_logger
from .peak_alignment import peak_matching, get_reference_mz_axis, get_reference_mz_axis2
from .filter_helper import (
    smoother,
)
from .normalizer_helper import normalizer
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
        data: Union[SpectrumBaseModule, SpectrumImzML],
        scale_method: str = 'none',
        method: str = "tic",
    ) -> Union[SpectrumBaseModule, SpectrumImzML]:
        """
        Normalize a single spectrum using TIC or median, with optional scaling.

        Parameters:
            data (SpectrumBaseModule | SpectrumImzML): Input spectrum whose `intensity` will be normalized.
            scale_method (str): Additional scaling after primary normalization:
                - 'none': no extra scaling
                - 'unit': min-max scaling to [0, 1]
            method (str): Primary normalization method:
                - 'tic': Total Ion Current normalization (sum equals 1)
                - 'median': Median normalization (median equals 1)

        Returns:
            SpectrumBaseModule: A new spectrum that retains the original `mz_list` and `coordinates`,
            with `intensity` replaced by the normalized (and optionally scaled) values.

        Raises:
            ValueError: If `method` is unsupported or if TIC/median is not greater than 0
            (propagated from the lower-level normalizer).
        """
        intensity = data.intensity
        norm_intensity = normalizer(
            intensity,
            scale_method=scale_method,
            method=method
        )

        return SpectrumBaseModule(
            mz_list=data.mz_list,
            intensity=norm_intensity,
            coordinates=data.coordinates,
        )

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
