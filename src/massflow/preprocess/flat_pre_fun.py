from typing import Optional
import numpy as np

from massflow.preprocess.helper.noise_reduction_helper import smoother
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess")


class FlatPreprocess:
    """Flat-array preprocessing helpers."""

    @staticmethod
    def noise_reduction_flat(
        intensity: np.ndarray,
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
        lengths: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Perform flat-mode noise reduction and return processed flat intensity.

        Supported methods are limited to:
        - ``savgol_numba``
        - ``gaussian_numba``
        - ``ma_numba``
        """
        method_norm = (method or "").strip().lower()
        supported_methods = {"savgol_numba", "gaussian_numba", "ma_numba"}
        if method_norm not in supported_methods:
            raise ValueError("noise_reduction_flat only supports: savgol_numba, gaussian_numba, ma_numba")

        return smoother(
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
