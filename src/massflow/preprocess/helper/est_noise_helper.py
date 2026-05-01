"""
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root

Noise estimation utilities.
"""
from typing import Optional
import numpy as np
from massflow.tools.logger import get_logger
from massflow.preprocess.numba.est_noise_numba import estimate_flat_numba

logger = get_logger("preprocesss")

def _input_validation(intensity,nbins,overlap):

    if intensity is None or intensity.ndim != 1:
        logger.error("x.intensity must be a 1D numpy array")
        raise TypeError("x.intensity must be a 1D numpy array")
    if nbins < 1:
        logger.error("nbins must be >= 1")
        raise ValueError("nbins must be >= 1")
    if overlap > 1 or overlap <0:
        logger.error("overlap must be in [0,1]")
        raise ValueError("overlap must be in [0,1]")


def _method_to_code(method: str) -> int:
    if method == 'sd':
        return 0
    if method == 'mad':
        return 1
    if method == 'quantile':
        return 2
    if method == 'diff':
        return 3
    raise ValueError(f"Unknown noise estimation method: {method}")


def estimator(
    intensity: np.ndarray,
    lengths: Optional[np.ndarray] = None,
    nbins: int = 1,
    overlap: float = 0.5,
    dynamic: bool = False,
    method: str = 'sd',
    denoise_method: str = 'gaussian_numba',
):
    """
    Estimate noise level in the MSI data.
    """
    _input_validation(intensity, nbins, overlap)
    lengths_arr = lengths.astype(np.int64) if lengths is not None else np.array([intensity.size], dtype=np.int64)

    method_code = _method_to_code(method)

    noise_flat = estimate_flat_numba(
        intensity=intensity,
        lengths=lengths_arr,
        nbins=nbins,
        dynamic=dynamic,
        overlap=overlap,
        method_code=method_code,
        denoise_method=denoise_method,
        floor_value=0.001,
    )

    return noise_flat
