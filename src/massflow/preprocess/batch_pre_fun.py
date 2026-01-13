from typing import Sequence
import numpy as np
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.spectrum import Spectrum
from massflow.tools.logger import get_logger
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.preprocess.helper.filter_helper import smoother

logger = get_logger("batch_pre_fun")

class BatchPreprocess:
    """
    A class for batch preprocessing of mass spectrometry data.
    """
    @staticmethod
    def peak_align_batch(batch_spectra: Sequence[Spectrum],
                         ref: np.ndarray,
                         tolerance: float,
                         units: str = "ppm") -> Sequence[SpectrumImzML]:
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
        if ref is None or tolerance is None:
            logger.error("Reference m/z axis and tolerance must be provided for alignment.")
            raise ValueError("Reference m/z axis and tolerance are required.")

        aligned_spectra = []
        for spectrum in batch_spectra:
            aligned_spectrum = SpectrumPreprocess.peak_align_spectrum(
                spectrum=spectrum,
                ref=ref,
                tolerance=tolerance,
                units=units
            )
            aligned_spectra.append(aligned_spectrum)
        return aligned_spectra

    @staticmethod
    def peak_pick_batch(batch_spectra: Sequence[Spectrum],
                        width: int | Sequence[int] = 2,
                        method: str = 'scipy',
                        relheight: float = 0.1,
                        return_type: str = 'height') -> Sequence[SpectrumImzML]:
        """
        Perform peak picking on a batch of spectra.

        Parameters:
        - spectra: Sequence of Spectrum objects to be processed.
        - width: Width parameter for peak picking.
        - method: Method for peak picking ('scipy', etc.).
        - relheight: Relative height threshold for peak picking.
        - return_type: Type of return value ('height' or 'area').

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
                return_type=return_type
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
        - method: One of {'ma','gaussian','savgol','wavelet','ma_ns','gaussian_ns','bi_ns'}.
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
        if not batch_spectra:
            return []

    
        # Specialized 2D Numba acceleration path for savgol_numba:
        # At the batch level, first stack spectra into a (n_spectra, n_mz) matrix,
        # then call the smoother once to trigger the parallel implementation
        # in savgol_batch_jit from noise_reduction_numba.
        if method == "savgol_numba":
            intensities = np.array([s.intensity for s in batch_spectra], dtype=np.float64)

            smoothed = smoother(
                intensities,
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

            denoised_spectra: list[SpectrumImzML] = []
            for i, spectrum in enumerate(batch_spectra):
                denoised_spectra.append(
                    SpectrumImzML(
                        coordinates=spectrum.coordinate,
                        mz_list=spectrum.mz_list,
                        intensity=smoothed[i],
                    )
                )

            return denoised_spectra

        # 其他 method（包括原始 savgol）仍然走单谱 1D 路径，保持 API 语义一致
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
