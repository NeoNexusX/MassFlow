from dataclasses import dataclass
from typing import Optional
import numpy as np

from massflow.preprocess.helper.noise_reduction_helper import smoother
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
    def baseline_correction_flat() -> FlatBatchResult:
        ...

    @staticmethod
    def noise_reduction_flat(
        mz_data: np.ndarray | None,
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
    def normalization_flat() -> FlatBatchResult:
        ...

    @staticmethod
    def peak_pick_flat() -> FlatBatchResult:
        ...

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

        aligned_2d = align_spectra_parallel( # 后续改成输出直接就是 1d 的形式，减少内存占用 TODO
            mz_data=mz_data,
            intensity=intensity,
            lengths=lengths,
            reference=reference,
            tolerance=tolerance,
            units=units,
        )

        ref_len = int(reference.size)
        out_intensity = aligned_2d.reshape(-1)

        return FlatBatchResult(mz_data=reference, intensity=out_intensity, lengths=ref_len)
