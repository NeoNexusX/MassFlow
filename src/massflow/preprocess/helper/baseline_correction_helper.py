from typing import Optional

import numpy as np

from massflow.preprocess.numba.baseline_correction_numba import (
    baseline_locmin_numba,
    baseline_snip_numba,
)
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.baseline")


def _input_validation(intensity: np.ndarray, baseline_scale: float = 1.0) -> None:
    """Validate shared flat baseline-correction inputs."""
    if not isinstance(intensity, np.ndarray):
        logger.error("intensity must be a numpy array")
        raise TypeError("intensity must be a numpy array")

    if intensity.ndim != 1 or intensity.size == 0:
        logger.error("intensity must be a non-empty 1D array")
        raise ValueError("intensity must be a non-empty 1D array")

    if not np.isfinite(baseline_scale) or not 0.0 < float(baseline_scale) <= 1.0:
        logger.error("baseline_scale must be a finite number in (0,1]")
        raise ValueError("baseline_scale must be a finite number in (0,1]")


def baseline_corrector(
    intensity: np.ndarray,
    method: str = "locmin_numba",
    smooth: str = "none",
    span: float = 0.1,
    upper: bool = False,
    width: int = 5,
    niter: int = 15,
    baseline_scale: float = 1.0,
    m: Optional[int] = None,
    decreasing: bool = True,
    lengths: Optional[np.ndarray] = None,
    return_baseline: bool = True,
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Flat-mode baseline estimation and correction.

    The current preprocessing pipeline supports only the parallel flat backends:
    ``locmin_numba`` and ``snip_numba``. Legacy single-spectrum Python methods
    are intentionally not routed from this helper.
    """
    _input_validation(intensity, baseline_scale=baseline_scale)
    target_dtype = intensity.dtype
    xi = np.asarray(intensity, dtype=np.float64)

    method_norm = (method or "locmin_numba").strip().lower()
    if method_norm == "locmin_numba":
        baseline = baseline_locmin_numba(
            intensity=xi,
            lengths=lengths,
            width=width,
            upper=upper,
            span=span,
            niter=niter,
            smooth=smooth,
        )
    elif method_norm == "snip_numba":
        baseline = baseline_snip_numba(
            intensity=xi,
            lengths=lengths,
            m=m,
            decreasing=decreasing,
        )
    else:
        raise ValueError("baseline_corrector only supports: locmin_numba, snip_numba")

    if baseline_scale != 1.0:
        baseline *= baseline_scale

    corrected = xi - baseline
    np.maximum(corrected, 0.0, out=corrected)

    return (
        corrected.astype(target_dtype, copy=False),
        baseline.astype(target_dtype, copy=False) if return_baseline else None,
    )
