from dataclasses import dataclass
from typing import Optional
import numpy as np
from massflow.preprocess.helper.baseline_correction_helper import baseline_corrector
from massflow.preprocess.helper.est_noise_helper import estimator
from massflow.preprocess.helper.noise_reduction_helper import smoother
from massflow.preprocess.helper.normalizer_helper import normalizer
from massflow.preprocess.helper.peak_align_parallel import align_spectra_parallel
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess")


@dataclass(slots=True)
class FlatBatchResult:
    """Result container for flat-array preprocessing steps."""
    mz_data: np.ndarray | None
    intensity: np.ndarray
    lengths: np.ndarray


class FlatPreprocess:
    """Flat-array preprocessing helpers."""

    @staticmethod
    def baseline_reduction_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
        method: str = "locmin_numba",
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
    ) -> FlatBatchResult:
        """Perform flat-mode baseline reduction and return corrected flat intensity.

        Supported methods are limited to:
        - ``locmin_numba``
        - ``snip_numba``
        """
        method_norm = (method or "").strip().lower()
        supported_methods = {"locmin_numba", "snip_numba"}
        if method_norm not in supported_methods:
            raise ValueError("baseline_reduction_flat only supports: locmin_numba, snip_numba")

        baselined_intensity, _ = baseline_corrector(
            intensity=intensity,
            method=method_norm,
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
            lengths=lengths,
        )

        return FlatBatchResult(mz_data=mz_data, intensity=baselined_intensity, lengths=lengths)

    @staticmethod
    def noise_reduction_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
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
    ) -> FlatBatchResult:
        method_norm = (method or "").strip().lower()
        supported_methods = {"savgol_numba", "gaussian_numba", "ma_numba"}
        if method_norm not in supported_methods:
            raise ValueError("noise_reduction_flat only supports: savgol_numba, gaussian_numba, ma_numba")

        smoothed_intensity = smoother(
            intensity=intensity,
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
        return FlatBatchResult(mz_data=mz_data, intensity=np.asarray(smoothed_intensity), lengths=lengths)

    @staticmethod
    def peak_pick_flat() -> FlatBatchResult:
        ...

    @staticmethod
    def est_nosise_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
        nbins: int = 1,
        overlap: float = 0.2,
        method: str = "sd",
        denoise_method: str = "gaussian_numba",
        dynamic: bool = False,
    ) -> FlatBatchResult:
        """Estimate flat-mode noise level and return estimated flat intensity."""
        index = mz_data if mz_data.size == intensity.size else None
        noise_estimation = estimator(
            intensity=intensity,
            indexes=index,
            lengths=lengths,
            nbins=nbins,
            overlap=overlap,
            dynamic=dynamic,
            method=method,
            denoise_method=denoise_method,
        )

        return FlatBatchResult(mz_data=mz_data, intensity=noise_estimation, lengths=lengths)

    @staticmethod
    def peak_align_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
        reference: Optional[np.ndarray] = None,
        tolerance: Optional[float] = None,
        units: str = "ppm",
    ) -> FlatBatchResult:
        if reference is None or tolerance is None:
            logger.error("Reference m/z axis and tolerance must be provided for alignment.")
            raise ValueError("Reference m/z axis and tolerance are required.")

        aligned_2d = align_spectra_parallel(
            mz_data=mz_data,
            intensity=intensity,
            lengths=lengths,
            reference=reference,
            tolerance=tolerance,
            units=units,
        )

        n_spec = int(lengths.size)
        ref_len = int(reference.size)
        out_lengths = np.full(n_spec, ref_len, dtype=np.int32)
        out_intensity = aligned_2d.reshape(-1)

        return FlatBatchResult(mz_data=reference, intensity=out_intensity, lengths=out_lengths)

    @staticmethod
    def normalization_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
        method: str = "tic_numba",
        scale: float | None = None,
        mz_flat: np.ndarray | None = None,
        ref: float | None = None,
        ref_tolerance: float = 0.1,
    ) -> FlatBatchResult:
        """Perform flat-mode normalization and return processed flat intensity.

        Supported methods are limited to:
        - ``tic_numba``
        - ``rms_numba``
        - ``ref_numba``

        """

        method_norm = (method or "").strip().lower()
        supported_methods = {"tic_numba", "rms_numba", "ref_numba"}
        if method_norm not in supported_methods:
            raise ValueError("normalization_flat only supports: tic_numba, rms_numba, ref_numba")

        normalized_intensity = normalizer(
            intensity=intensity,
            method=method_norm,
            scale=scale,
            mz_flat=mz_flat,
            ref=ref,
            ref_tolerance=ref_tolerance,
            lengths=lengths,
        )

        return FlatBatchResult(mz_data=mz_data, intensity=normalized_intensity, lengths=lengths)
