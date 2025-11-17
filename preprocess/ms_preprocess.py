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
from .peak_align_helper import peak_align, get_reference
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
    def peak_alignment(ms_data: MS,
                       ref: Optional[np.ndarray] = None,
                       units: str = 'ppm',
                       tolerance: Optional[float] = None,
                       binfun: str = 'min',
                       binratio: int = 2,
                       # 2
                       ref_method: Optional[str] = None,
                       ref_mzres: Optional[float] = None,
                       ref_mzmaxshift: Optional[float] = None,
                       ref_mzunits: Optional[str] = None) -> MS:
        """
        Peak alignment entry point: collect per-spectrum peak lists, choose reference axis source,
        call the core `peak_align`, and build a new MS with a unified reference axis.

        Parameters:
            ms_data (MS): unaligned MS collection
            ref (np.ndarray|None): external reference axis; None means internal estimation or external computation via `ref_method`
            units (str): tolerance units, 'ppm' or 'mz'/'absolute'
            tolerance (float|None): tolerance half-window; None means auto-estimation
            binfun (str): resolution aggregation strategy, 'median'|'min'|'max'|'mean'
            binratio (int): tolerance scaling factor, default 2
            ref_method (str|None): reference axis source; 'density' uses `get_reference` to compute the axis
            ref_mzres (float|None): histogram resolution for reference axis; defaults mapped as ppm→100.0, Da→0.01
            ref_mzmaxshift (float|None): peak detection window for the reference axis (gap for `__findpeaks`); defaults mapped as ppm→100.0, Da→0.05
            ref_mzunits (str|None): units used for reference axis computation; defaults derived from `units` (ppm→'ppm'; mz/absolute→'Da')

        Workflow:
            1) Collect m/z and intensity for each spectrum; validate length equality and non-empty
            2) Choose reference axis:
               - If `ref` is provided, use it directly
               - If `ref_method=='density'`, concatenate all spectra m/z and call `get_reference(mz_all, ref_mzres, ref_mzmaxshift, ref_mzunits)`
               - Otherwise, estimate the reference axis internally via the alignment logic
            3) Call `peak_align` to obtain the alignment matrix
            4) Build a new MS: unify all spectra m/z to the reference axis and fill intensities from aligned columns
            5) Update necessary metadata (`processed`)
        """

        logger.info(
            f"peak_alignment: ref={'provided' if ref is not None else 'None'}, units={units}, tolerance={tolerance}, binfun={binfun}, binratio={binratio}, spectra={len(ms_data)}, ref_method={ref_method}")
        # 1) collect per-spectrum index and data lists
        index_list = []
        data_list = []

        # validate equal lengths and non-empty arrays
        for spectrum in ms_data:
            mz = spectrum.mz_list
            it = spectrum.intensity
            if mz is None or it is None:
                logger.error('spectrum has no mz_list/intensity loaded')
                raise ValueError('spectrum has no mz_list/intensity loaded')
            if mz.size != it.size:
                logger.error('lengths of mz_list and intensity must match')
                raise ValueError('lengths of mz_list and intensity must match')
            index_list.append(mz)
            data_list.append(it)

        if ref is None and (ref_method == 'density'):
            if len(index_list) == 0:
                logger.error('peak_alignment: no spectra available for keneng reference computation')
                raise ValueError('no spectra available for reference computation')
            mz_all = np.concatenate(index_list)
            if ref_mzunits is None:
                ref_mzunits = 'ppm' if units == 'ppm' else 'Da'
            if ref_mzres is None:
                ref_mzres = float(tolerance) if tolerance is not None else (100.0 if ref_mzunits == 'ppm' else 0.01)
            if ref_mzmaxshift is None:
                ref_mzmaxshift = float(tolerance) if tolerance is not None else (
                    100.0 if ref_mzunits == 'ppm' else 0.05)
            ref = get_reference(mz_all, ref_mzres, ref_mzmaxshift, ref_mzunits)
            ref = ref[np.isfinite(ref)]
            logger.info(
                f"peak_alignment: keneng ref computed len={ref.size}, ref_units={ref_mzunits}, mzres={ref_mzres}, maxshift={ref_mzmaxshift}, fused_from_tolerance={tolerance is not None}")

        # call peak_align to perform alignment
        res = peak_align(
            index_list,
            data_list,
            ref=ref,
            binfun=binfun,
            binratio=binratio,
            tolerance=tolerance,
            units=units)
        logger.info(
            f"peak_alignment: ref_len={len(res.ref)}, tolerance={res.tolerance}, units={res.units}, spectra={len(index_list)}")

        # build ms_aligned and populate data; copy coordinates and essential metadata tolerance/binfun/binratio
        ms_aligned = MS()
        ms_aligned.meta = ms_data.meta
        for col, spectrum in enumerate(ms_data):
            aligned = res.matrix[:, col]
            new_spec = SpectrumBaseModule(mz_list=res.ref, intensity=aligned, coordinates=spectrum.coordinates,
                                          sorted_by_mz_fun=True)
            ms_aligned.add_spectrum(new_spec)

        logger.info(f"ms_aligned built: B={len(ms_aligned[0])}, S={len(ms_aligned)}")
        ms_aligned.meta.processed = True

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
        window: int = 3,
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
