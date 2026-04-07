from typing import Sequence
import numpy as np
from massflow.module import SpectrumImzML, Spectrum
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.preprocess.helper.noise_reduction_helper import smoother
from massflow.preprocess.helper.normalizer_helper import normalizer

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
        method: str = "ma_numba",
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
        - method: One of {'ma_numba', 'gaussian_numba', 'savgol_numba'}.
        - window: Smoothing window size.
        - sd: Gaussian scale parameter (for gaussian_numba).
        - sd_intensity: Reserved for API compatibility.
        - p: Reserved for API compatibility.
        - coef: Reserved for API compatibility.
        - polyorder: Polynomial order for Savitzky-Golay.
        - deriv: Derivative order for Savitzky-Golay.
        - delta: Sample spacing for Savitzky-Golay.
        - wavelet: Reserved for API compatibility.
        - threshold_mode: Reserved for API compatibility.

        Returns:
        - Sequence of denoised spectra as SpectrumImzML objects.
        """
        if not batch_spectra:
            return []

        method_norm = (method or "").strip().lower()
        supported_methods = {"savgol_numba", "gaussian_numba", "ma_numba"}
        if method_norm not in supported_methods:
            raise ValueError("noise_reduction_batch only supports: savgol_numba, gaussian_numba, ma_numba")

        lengths = np.array([s.intensity.size for s in batch_spectra], dtype=np.int64)  # type: ignore
        total_len = int(np.sum(lengths))
        intensity_flat = np.empty(total_len, dtype=np.float32)

        offset = 0
        for spectrum, valid_len in zip(batch_spectra, lengths):
            end = offset + int(valid_len)
            intensity_flat[offset:end] = spectrum.intensity.astype(np.float32, copy=False)  # type: ignore
            offset = end

        smoothed_flat = smoother(
            intensity_flat,
            method=method_norm,
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
            lengths=lengths,
        )

        denoised_spectra: list[SpectrumImzML] = []
        offset = 0
        for spectrum, valid_len in zip(batch_spectra, lengths):
            end = offset + int(valid_len)
            denoised_intensity = smoothed_flat[offset:end].copy()
            denoised_spectra.append(
                SpectrumImzML(
                    coordinates=spectrum.coordinate,
                    mz_list=spectrum.mz_list,
                    intensity=denoised_intensity,
                )
            )
            offset = end

        return denoised_spectra

    @staticmethod
    def normalization_batch(
        batch_spectra: Sequence[Spectrum],
        scale_method: str = "none",
        method: str = "tic",
        scale: float = 1.0,
    ) -> Sequence[SpectrumImzML]:
        """Normalize a batch of spectra.

        Parameters:
        - batch_spectra: Sequence of Spectrum objects to be normalized.
        - scale_method: 'none' or 'unit' min-max scaling.
        - method: One of {'tic', 'rms', 'median'}.
        - scale: Cardinal-like amplitude scaling factor applied after normalization.

        Returns:
        - Sequence of normalized spectra as SpectrumImzML objects.
        """
        if not batch_spectra:
            return []

        # -----------------------------------------------------------------
        # Path A: High Performance Batching (Numba)
        # Try to use the direct Numba backend if available and method is supported.
        # -----------------------------------------------------------------
        numba_supported_methods = {"tic_numba", "rms_numba", "median_numba"}

        # Check availability of Numba implementation locally
        if method in numba_supported_methods:
            # 1. Construct 2D Matrix with padding
            n_spectra = len(batch_spectra)

            # Calculate lengths for each spectrum
            lengths = np.array([s.intensity.size for s in batch_spectra], dtype=np.int32)  # type: ignore
            max_len = np.max(lengths) if n_spectra > 0 else 0

            if max_len > 0:
                intensities_padded = np.zeros((n_spectra, max_len), dtype=np.float32)

                # Fill the matrix
                for i, s in enumerate(batch_spectra):
                    valid_len = lengths[i]
                    if valid_len > 0:
                        intensities_padded[i, :valid_len] = s.intensity.astype(np.float32, copy=False)  # type: ignore

                # 2. Call the Numba kernel directly (normalizer_numba)
                # It accepts the 2D matrix and 'lengths' for efficient batch processing
                normalized_matrix = normalizer(
                    intensities_padded,
                    method=method,
                    scale=scale,
                    scale_method=scale_method,
                    lengths=lengths,
                )

                # 3. Reconstruct SpectrumImzML objects
                normalized_spectra_numba: list[SpectrumImzML] = []
                for i, spectrum in enumerate(batch_spectra):
                    valid_len = lengths[i]

                    # Slice back to original length
                    processed_intensity = normalized_matrix[i, :valid_len]

                    normalized_spectra_numba.append(
                        SpectrumImzML(
                            coordinates=spectrum.coordinate,
                            mz_list=spectrum.mz_list,
                            intensity=processed_intensity,
                        )
                    )
                return normalized_spectra_numba

        # -----------------------------------------------------------------
        # Path B: Fallback (Sequential Loop)
        # Delegate to SpectrumPreprocess.normalization_spectrum for each spectrum.
        # -----------------------------------------------------------------
        normalized_spectra: list[SpectrumImzML] = []
        for spectrum in batch_spectra:
            # Re-use the single spectrum processing logic
            normalized_spectrum = SpectrumPreprocess.normalization_spectrum(
                data=spectrum, method=method, scale=scale, scale_method=scale_method
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
