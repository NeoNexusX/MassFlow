from typing import Sequence
import numpy as np
from massflow.module.spectrum_imzml import SpectrumImzML
from massflow.module.spectrum import Spectrum
from massflow.tools.logger import get_logger
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess

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
