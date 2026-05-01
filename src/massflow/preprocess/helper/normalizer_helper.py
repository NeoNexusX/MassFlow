from typing import Optional

import numpy as np

from massflow.preprocess.numba.normalization_numba import normalizer_numba
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.normalization")


def _input_validation(
    intensity: np.ndarray,
    scale: Optional[float] = None,
    ref_tolerance: float = 0.1,
    method_norm: Optional[str] = None,
    mz_flat: Optional[np.ndarray] = None,
    ref: Optional[float] = None,
) -> None:
    """Validate shared flat normalization inputs."""
    if intensity is None or not isinstance(intensity, np.ndarray):
        logger.error("intensity must be a numpy array")
        raise TypeError("intensity must be a numpy array")

    if intensity.ndim != 1 or intensity.size == 0:
        logger.error("intensity must be a non-empty 1D array")
        raise ValueError("intensity must be a non-empty 1D array")

    if scale is not None and (not np.isfinite(scale) or float(scale) < 0.0):
        logger.error("scale must be a finite non-negative number")
        raise ValueError("scale must be a finite non-negative number")

    if not np.isfinite(ref_tolerance) or float(ref_tolerance) < 0.0:
        logger.error("ref_tolerance must be a finite non-negative number")
        raise ValueError("ref_tolerance must be a finite non-negative number")

    if method_norm == "ref_numba":
        if mz_flat is None:
            logger.error("mz_flat is required when method='ref_numba'")
            raise ValueError("mz_flat is required when method='ref_numba'")
        if ref is None:
            logger.error("ref is required when method='ref_numba'")
            raise ValueError("ref is required when method='ref_numba'")


def normalizer(
    intensity: np.ndarray,
    method: str = "tic_numba",
    scale: Optional[float] = None,
    mz_flat: Optional[np.ndarray] = None,
    ref: Optional[float] = None,
    ref_tolerance: float = 0.1,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Flat-mode normalization dispatcher.

    The current preprocessing API routes only the parallel flat numba backends:
    ``tic_numba``, ``rms_numba`` and ``ref_numba``.
    """
    method_norm = (method or "tic_numba").strip().lower()
    if method_norm not in {"tic_numba", "rms_numba", "ref_numba"}:
        raise ValueError("normalizer only supports: tic_numba, rms_numba, ref_numba")

    _input_validation(
        intensity,
        scale=scale,
        ref_tolerance=ref_tolerance,
        method_norm=method_norm,
        mz_flat=mz_flat,
        ref=ref,
    )

    return normalizer_numba(
        intensity=intensity,
        method=method_norm.replace("_numba", ""),
        scale=scale,
        mz_flat=mz_flat,
        ref=ref,
        ref_tolerance=ref_tolerance,
        lengths=lengths,
    )
