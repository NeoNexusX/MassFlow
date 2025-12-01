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
from massflow.module.ms_module import SpectrumBaseModule, SpectrumImzML, MS
from massflow.logger import get_logger
from .peak_align_helper import peak_align, get_reference
from .filter_helper import smoother
from .normalizer_helper import normalizer
from .baseline_correction_helper import baseline_corrector
from .est_noise_helper import estimator
from .peak_pick_helper import peak_pick_fun
import os
import time

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
    def base_input_check(data):
        """
        Validate input data.

        Raises:
            ValueError: If input data is invalid or missing required attributes
        """
        if data.intensity is None or data.mz_list is None:
            raise ValueError("Input spectrum must have both intensity and mz_list data.")

    @staticmethod
    def peak_pick_spectrum( data: SpectrumBaseModule,
                            width: int | Sequence[int] = 2,
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
        MSIPreprocessor.base_input_check(data=data)

        intensity = data.intensity
        index = data.mz_list
        peak_intensity,peak_index = peak_pick_fun(  intensity, # type: ignore
                                                    index, # type: ignore
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
        scale_method: str = "none",
        method: str = "tic",
        scale: float = 1.0
    ) -> Union[SpectrumBaseModule, SpectrumImzML]:
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
            scale_method=scale_method,
            method=method,
            scale=scale
        )

        return SpectrumBaseModule(
            mz_list=data.mz_list,
            intensity=norm_intensity,
            coordinates=data.coordinates,
        )

    @staticmethod
    def peak_alignment(ms_data: MS,
                       ref: Optional[np.ndarray] = None,
                       units: str = 'ppm',
                       tolerance: Optional[float] = None,
                       binfun: str = 'median',
                       binratio: int = 2,
                       output_path: Optional[str]=None,
                       ref_method: Optional[str] = None,
                       ref_mzres: Optional[float] = None,
                       ref_mzmaxshift: Optional[float] = None,
                       ref_mzunits: Optional[str] = None) -> MS:
        """
        Peak alignment entry (streaming): dispatch parameters, optionally compute density-based reference,
        delegate to streaming `peak_align`, and return the aligned MS object.

        Returns:
            MS: Aligned collection whose spectra share the reference axis.

        Raises:
            ValueError: When input spectra are empty or density reference computation fails.

        Notes:
            Units mapping follows tolerance semantics: 'ppm'→relative, 'mz'/'absolute'→absolute (Da).
        """
        logger.info(f"peak_alignment_entry: binfun={binfun}, units={units}...")
        if output_path is None:
            ts = time.strftime("%Y%m%d-%H%M%S")
            default_name = f"align_ms_{ts}.h5"
            output_path = default_name
            input_file = getattr(ms_data.meta, 'filepath', None) if ms_data.meta else None
            if input_file:
                input_dir = os.path.dirname(input_file)
                candidate = os.path.join(input_dir, default_name)
                try:
                    tmp = os.path.join(input_dir, ".mf_write_test.tmp")
                    with open(tmp, 'w'):
                        pass
                    os.remove(tmp)
                    output_path = candidate
                except(OSError, IOError, PermissionError):
                    logger.warning(f"Cannot write to input directory {input_dir}, saving result to {output_path}")
        
        logger.info(f"Peak alignment output will be saved to: {output_path}")

        # 1. Handle Density-based Reference (Special Case)
        # If user wants density method, we calculate ref first, then pass to peak_align
        if ref is None and (ref_method == 'density'):
            # We still need to collect mz_list for density calculation
            # Alternatively, peak_align could handle this if we moved get_reference logic inside
            # But keeping it here is fine for separation of concerns
            index_list = [s.mz_list for s in ms_data]
            if len(index_list) == 0:
                raise ValueError('no spectra available')
            mz_all = np.concatenate(index_list)
            
            # ... (Mapping params logic keeps same) ...
            if ref_mzunits is None: ref_mzunits = 'ppm' if units == 'ppm' else 'Da'
            if ref_mzres is None: ref_mzres = float(tolerance) if tolerance is not None else (100.0 if ref_mzunits == 'ppm' else 0.01)
            if ref_mzmaxshift is None: ref_mzmaxshift = float(tolerance) if tolerance is not None else (100.0 if ref_mzunits == 'ppm' else 0.05)
            
            ref = get_reference(mz_all, ref_mzres, ref_mzmaxshift, ref_mzunits)
            ref = ref[np.isfinite(ref)]
            logger.info(f"peak_alignment: density ref computed, peaks={ref.size}")

        # 2. Call peak_align (Now Handles Everything)
        # Pass ms_data directly
        res = peak_align(
            ms_data,
            ref=ref,
            binfun=binfun,
            binratio=binratio,
            tolerance=tolerance,
            units=units,
            output_path=output_path
        )

        # 3. Return the constructed MS object
        ms_aligned = res.ms_aligned
        if ms_aligned.meta is not None and hasattr(ms_aligned.meta, 'continuous'):
            ms_aligned.meta.continuous = True

        logger.info(f"peak_alignment_entry: done, peaks={len(res.ref)}")

        return ms_aligned

    @staticmethod
    def baseline_correction_spectrum(
        data: Union[SpectrumBaseModule, SpectrumImzML],
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
    ) -> tuple[SpectrumBaseModule, np.ndarray]:
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
        corrected_spectrum = SpectrumBaseModule(
            mz_list=data.mz_list,
            intensity=corrected_intensity,
            coordinates=data.coordinates,
        )
        
        return corrected_spectrum, baseline

    @staticmethod
    def noise_reduction_spectrum(
        data: Union[SpectrumBaseModule, SpectrumImzML],
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
        MSIPreprocessor.base_input_check(data=data)

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
        MSIPreprocessor.base_input_check(data=x)

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
        MSIPreprocessor.base_input_check(data=spectrum)

        signal_level = np.percentile(spectrum.intensity, 95)# type: ignore

        noise = MSIPreprocessor.noise_estimation_spectrum(spectrum, method=method)

        logger.info(f"SNR: signal_level:{signal_level}, noise:{noise}")
        return signal_level / noise

    @staticmethod
    def preprocess_pipeline():
        """
        Composite preprocessing pipeline (placeholder).

        Parameters:
            data (SpectrumBaseModule): Input spectrum to preprocess.

        Returns:
            SpectrumBaseModule: Preprocessed spectrum (implementation-defined).
        """

