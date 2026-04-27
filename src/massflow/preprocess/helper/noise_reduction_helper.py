from typing import Optional

import numpy as np

from massflow.preprocess.numba.noise_reduction_numba import (
    smooth_signal_gaussian_numba,
    smooth_signal_ma_numba,
    smooth_signal_savgol_numba,
)
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.noise_reduction")


def _input_validation(
    intensity: np.ndarray,
    *,
    window: Optional[int] = None,
    sd: Optional[float] = None,
) -> None:
    """Validate shared flat smoothing inputs."""
    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1 or intensity.size == 0:
        logger.error("intensity must be a non-empty 1D numpy array")
        raise TypeError("intensity must be a non-empty 1D numpy array")

    if window is not None and (not isinstance(window, int) or window <= 0):
        logger.error("window must be a positive integer")
        raise ValueError("window must be a positive integer")

    if sd is not None and sd <= 0:
        logger.error("sd must be positive")
        raise ValueError("sd must be positive")


def smoother(
    intensity: Optional[np.ndarray] = None,
    *,
    flat: Optional[np.ndarray] = None,
    method: str = "ma_numba",
    window: int = 5,
    sd: Optional[float] = None,
    polyorder: int = 2,
    deriv: int = 0,
    delta: float = 1.0,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Flat-mode smoothing dispatcher.

    Only the parallel flat numba backends used by ``FlatPreprocess`` are
    supported here: ``ma_numba``, ``gaussian_numba`` and ``savgol_numba``.
    Parameters kept for the public wrapper signature are accepted but ignored
    when they are not meaningful for the selected backend.
    """

    method_norm = (method or "ma_numba").strip().lower()
    input_signal = flat if flat is not None else intensity
    if input_signal is None:
        raise TypeError("Either flat or intensity must be provided")

    _input_validation(input_signal, window=window, sd=sd)

    if method_norm == "ma_numba":
        return smooth_signal_ma_numba(
            intensity=input_signal,
            window=window,
            lengths=lengths,
        )

    if method_norm == "gaussian_numba":
        return smooth_signal_gaussian_numba(
            intensity=input_signal,
            window=window,
            sd=sd,
            lengths=lengths,
        )

    if method_norm == "savgol_numba":
        return smooth_signal_savgol_numba(
            intensity=input_signal,
            window=window,
            polyorder=polyorder,
            deriv=deriv,
            delta=delta,
            lengths=lengths,
        )

    raise ValueError("smoother only supports: ma_numba, gaussian_numba, savgol_numba")
