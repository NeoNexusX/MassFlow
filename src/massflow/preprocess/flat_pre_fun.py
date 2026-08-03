"""Flat-array preprocessing functions for mass spectrometry data."""

from typing import Literal, Optional
from dataclasses import dataclass

import numpy as np

from massflow.preprocess.helper.baseline_correction_helper import baseline_corrector
from massflow.preprocess.helper.est_noise_helper import estimator
from massflow.preprocess.helper.noise_reduction_helper import smoother
from massflow.preprocess.helper.normalizer_helper import normalizer
from massflow.preprocess.helper.peak_pick_helper import peak_picker
from massflow.preprocess.helper.peak_align_helper import peak_aligner
from massflow.preprocess.helper.isotope_filter_helper import filter_peaks_by_isotopes
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
        upper: bool = False,
        width: int = 5,
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
            upper=upper,
            width=width,
            niter=niter,
            baseline_scale=baseline_scale,
            m=m,
            decreasing=decreasing,
            lengths=lengths,
            return_baseline=False,
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
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0,
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
            polyorder=polyorder,
            deriv=deriv,
            delta=delta,
            lengths=lengths,
        )
        return FlatBatchResult(mz_data=mz_data, intensity=np.asarray(smoothed_intensity), lengths=lengths)

    @staticmethod
    def peak_pick_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
        width: int = 5,
        method: str = "mad",
        snr: float = 2.0,
        return_type: str = "height",
        prominence: float | None = None,
        relheight: float | None = None,
        nbins: int = 1,
        overlap: float = 0.5,
    ) -> FlatBatchResult:

        if method not in {"sd", "mad", "quantile", "diff"}:
            raise ValueError("method must be one of 'sd', 'mad', 'quantile', 'diff'")

        if nbins < 1:
            raise ValueError("nbins must be >= 1")
        if overlap < 0.0 or overlap > 0.99:
            raise ValueError("overlap must be in [0.0, 0.99]")
        if return_type not in {"height", "area"}:
            raise ValueError("return_type must be 'height' or 'area'")

        picked_mz, picked_intensity, picked_lengths = peak_picker(
            mz_data=mz_data,
            intensity=intensity,
            lengths=lengths,
            width=width,
            method=method,
            snr=snr,
            return_type=return_type,
            prominence=prominence,
            relheight=relheight,
            nbins=nbins,
            overlap=overlap,
        )

        return FlatBatchResult(mz_data=picked_mz, intensity=picked_intensity, lengths=picked_lengths)

    @staticmethod
    def isotope_filter_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
        mz_tolerance: float = 10.0,
        units: Literal["da", "ppm"] = "ppm",
        min_charge: int = 1,
        max_charge: int = 3,
        min_isotope_count: int = 3,
        max_isotope_count: int = 10,
        use_decreasing_model: bool = True,
        start_intensity_check: int = 2,
    ) -> FlatBatchResult:
        """Remove isotope peaks independently from each centroided spectrum."""
        filtered_mz, filtered_intensity, filtered_lengths = filter_peaks_by_isotopes(
            mz_data,
            intensity,
            lengths,
            mz_tolerance=mz_tolerance,
            units=units,
            min_charge=min_charge,
            max_charge=max_charge,
            min_isotope_count=min_isotope_count,
            max_isotope_count=max_isotope_count,
            use_decreasing_model=use_decreasing_model,
            start_intensity_check=start_intensity_check,
        )
        return FlatBatchResult(
            mz_data=filtered_mz,
            intensity=filtered_intensity,
            lengths=filtered_lengths,
        )

    @staticmethod
    def est_noise_flat(
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
        noise_estimation = estimator(
            intensity=intensity,
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

        aligned_2d = peak_aligner(
            mz_data=mz_data,
            intensity=intensity,
            lengths=lengths,
            reference=reference,
            tolerance=tolerance,
            units=units,
        )

        n_spec = int(lengths.size)
        ref_len = int(reference.size)
        out_lengths = np.full(n_spec, ref_len, dtype=np.uint32)
        out_intensity = aligned_2d.reshape(-1)

        return FlatBatchResult(mz_data=np.asarray(reference, dtype=mz_data.dtype), intensity=out_intensity, lengths=out_lengths)

    @staticmethod
    def normalization_flat(
        mz_data: np.ndarray,
        intensity: np.ndarray,
        lengths: np.ndarray,
        method: str = "tic_numba",
        scale: float | None = None,
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

        if ref is None and method_norm == "ref_numba":
            ref = float(mz_data[mz_data.size // 2]) if mz_data is not None and mz_data.size > 0 else None

        normalized_intensity = normalizer(
            intensity=intensity,
            method=method_norm,
            scale=scale,
            mz_flat=mz_data if method_norm == "ref_numba" else None,
            ref=ref,
            ref_tolerance=ref_tolerance,
            lengths=lengths,
        )

        return FlatBatchResult(mz_data=mz_data, intensity=normalized_intensity, lengths=lengths)
