from typing import Sequence
import numpy as np
from massflow.module import SpectrumImzML, Spectrum
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess

from massflow.tools.logger import get_logger
logger = get_logger("massflow.preprocess")


class BatchPreprocess:
    """
    A class for batch preprocessing of mass spectrometry data.
    """

    @staticmethod
    def peak_align_batch(
        batch_spectra: Sequence[Spectrum],
        reference: np.ndarray,
        tolerance: float,
        units: str = "ppm",
    ) -> Sequence[SpectrumImzML]:
        """
        Align a batch of spectra to a reference m/z axis.

        Parameters:
        - spectra: Sequence of Spectrum objects to be aligned.
        - reference: Reference m/z axis for alignment.
        - tolerance: Tolerance for peak alignment.
        - units: Units for tolerance ('ppm' or 'mz').

        Returns:
        - Aligned spectra as a sequence of SpectrumImzML objects.
        """
        if reference is None or tolerance is None:
            logger.error(
                "Reference m/z axis and tolerance must be provided for alignment."
            )
            raise ValueError("Reference m/z axis and tolerance are required.")

        aligned_spectra = []
        for spectrum in batch_spectra:
            aligned_spectrum = SpectrumPreprocess.peak_align_spectrum(
                spectrum=spectrum, reference=reference, tolerance=tolerance, units=units
            )
            aligned_spectra.append(aligned_spectrum)
        return aligned_spectra

    @staticmethod
    def peak_pick_batch(
        batch_spectra: Sequence[Spectrum],
        width: int | Sequence[int] = 2,
        method: str = "origin",
        relheight: float = 0.012,
        snr: float = 2.0,
        return_type: str = "height",
        use_numba: bool = True,
    ) -> Sequence[SpectrumImzML]:
        """
        Perform peak picking on a batch of spectra.

        Parameters:
        - spectra: Sequence of Spectrum objects to be processed.
        - width: Width parameter for peak picking.
        - method: Method for peak picking ('origin', etc.).
        - relheight: Relative height threshold for peak picking.
        - snr: Signal-to-noise ratio threshold.
        - return_type: Type of return value ('height' or 'area').
        - use_numba: Whether to use Numba acceleration.

        Returns:
        - Processed spectra as a sequence of Spectrum objects.
        """
        picked_spectra = []
        for spectrum in batch_spectra:
            picked_spectrum = SpectrumPreprocess.peak_pick_spectrum(
                data=spectrum,
                width=width,
                method=method,
                relheight=relheight,
                snr=snr,
                return_type=return_type,
                use_numba=use_numba,
            )
            picked_spectra.append(picked_spectrum)
        return picked_spectra

    @staticmethod
    def noise_reduction_batch(
        batch_spectra: Sequence[Spectrum],
        method: str = "ma",
        window: int = 5,
        sd: float | None = None,
        sd_intensity: float | None = None,
        p: int = 2,
        coef: np.ndarray | None = None,
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0,
        wavelet: str = "db4",
        threshold_mode: str = "soft",
    ) -> Sequence[SpectrumImzML]:
        """Perform noise reduction on a batch of spectra.

        Parameters:
        - batch_spectra: Sequence of Spectrum objects to be denoised.
        - method: One of {'ma','ma_numba','gaussian','gaussian_numba','savgol',
            'savgol_numba','wavelet','gaussian_ns','gaussian_ns_numba','bi_ns','bi_ns_numba'}.
        - window: Window size or neighbor count depending on method.
        - sd: Gaussian scale parameter.
        - sd_intensity: Intensity scale for bilateral method.
        - p: Minkowski metric for NS queries.
        - coef: Custom kernel for 'ma'.
        - polyorder: Polynomial order for Savitzky-Golay.
        - deriv: Derivative order for Savitzky-Golay.
        - delta: Sample spacing for Savitzky-Golay.
        - wavelet: Wavelet family for wavelet denoising.
        - threshold_mode: 'soft' or 'hard' thresholding.

        Returns:
        - Sequence of denoised spectra as SpectrumImzML objects.
        """

        # Other methods (including the original savgol) still use the per-spectrum 1D path to keep the API semantics consistent
        denoised_spectra: list[SpectrumImzML] = []
        for spectrum in batch_spectra:
            denoised = SpectrumPreprocess.noise_reduction_spectrum(
                data=spectrum,
                method=method,
                window=window,
                sd=sd,
                sd_intensity=sd_intensity,
                p=p,
                coef=coef,
                polyorder=polyorder,
                deriv=deriv,
                delta=delta,
                wavelet=wavelet,
                threshold_mode=threshold_mode,
            )
            denoised_spectra.append(denoised)

        return denoised_spectra

    @staticmethod
    def normalization_batch(
        batch_spectra: Sequence[Spectrum],
        method: str = "tic",
        scale: float = 1.0,
    ) -> Sequence[SpectrumImzML]:
        """Normalize a batch of spectra.

        Parameters:
        - batch_spectra: Sequence of Spectrum objects to be normalized.
        - method: One of {'tic', 'rms', 'median'}.
        - scale: Cardinal-like amplitude scaling factor applied after normalization.

        Returns:
        - Sequence of normalized spectra as SpectrumImzML objects.
        """
        if not batch_spectra:
            return []

        method_norm = (method or "tic").strip().lower()
        supported_methods = {"tic", "rms", "median"}
        if method_norm not in supported_methods:
            logger.error("normalization_batch only supports tic/rms/median. ")
            raise ValueError("normalization_batch only supports: tic, rms, median. ")

        normalized_spectra: list[SpectrumImzML] = []
        for spectrum in batch_spectra:
            normalized_spectrum = SpectrumPreprocess.normalization_spectrum(
                data=spectrum,
                method=method_norm,
                scale=scale,
            )
            normalized_spectra.append(normalized_spectrum)

        return normalized_spectra

    @staticmethod
    def baseline_correction_batch(
        batch_spectra: Sequence[Spectrum],
        method: str = "asls",
        smooth: str = "none",
        span: float = 0.1,
        s: float | None = 0.0,
        upper: bool = False,
        width: int = 5,
        lam: float = 1e7,
        p: float = 0.01,
        niter: int = 15,
        baseline_scale: float = 1.0,
        m: int | None = None,
        decreasing: bool = True,
    ) -> Sequence[SpectrumImzML]:
        """Perform baseline correction on a batch of spectra.

        Parameters:
        - batch_spectra: Sequence of Spectrum objects to be baseline-corrected.
        - method: One of {'locmin', 'snip', 'asls'}.
        - smooth: LocMin smoothing method; one of {'none', 'loess', 'spline'}.
        - span: Loess span proportion in (0, 1]; used when method='locmin'.
        - s: Spline smoothing target RSS (>= 0); used when method='locmin'.
        - upper: If True, use local maxima as anchors; used when method='locmin'.
        - width: Neighborhood width (>= 3) for extrema detection; used when method='locmin'.
        - lam: ASLS smoothness parameter; used when method='asls'.
        - p: ASLS asymmetry parameter in (0, 1); used when method='asls'.
        - niter: ASLS iteration count (> 0); used when method='asls'.
        - baseline_scale: Scale factor in (0, 1] applied to estimated baseline.
        - m: SNIP window half-size (>= 1); used when method='snip'.
        - decreasing: SNIP decreasing rule; used when method='snip'.

        Returns:
        - Sequence of baseline-corrected spectra.
        """
        if not batch_spectra:
            return []

        corrected_spectra: list[SpectrumImzML] = []
        for spectrum in batch_spectra:
            corrected = SpectrumPreprocess.baseline_correction_spectrum(
                data=spectrum,
                method=method,
                smooth=smooth,
                span=span,
                s=s,
                upper=upper,
                width=width,
                lam=lam,
                p=p,
                niter=niter,
                baseline_scale=baseline_scale,
                m=m,
                decreasing=decreasing,
            )
            corrected_spectra.append(corrected)

        return corrected_spectra
