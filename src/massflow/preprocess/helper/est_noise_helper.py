"""
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root

Noise estimation utilities.
"""
from typing import Optional
import numpy as np
from scipy import stats
from massflow.tools.logger import get_logger
from massflow.preprocess.numba.est_noise_numba import estimate_flat_numba

logger = get_logger("preprocesss")

def _input_validation(intensity,index,nbins,overlap):

    if intensity is None or intensity.ndim != 1:
        logger.error("x.intensity must be a 1D numpy array")
        raise TypeError("x.intensity must be a 1D numpy array")
    if index is not None and (index.ndim !=1):
        logger.error("index must be a 1D numpy array of the same length as intensity")
        raise TypeError("index must be a 1D numpy array of the same length as intensity")
    if nbins < 1:
        logger.error("nbins must be >= 1")
        raise ValueError("nbins must be >= 1")
    if overlap > 1 or overlap <0:
        logger.error("overlap must be in [0,1]")
        raise ValueError("overlap must be in [0,1]")


def estimation_fun(data,method: str = 'sd'):
    """
    Return a callable for noise estimation from a vector.
    """

    if method == 'sd':
        return np.nanstd(data)
    elif method == 'mad':
        return stats.median_abs_deviation(data)
    elif method == 'quantile':
        return np.nanquantile(data,0.95)
    elif method == 'diff':
        return np.nanmean(np.abs(data - np.nanmean(data)))
    else:
        raise ValueError(f"Unknown noise estimation method: {method}")


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


def estimator(intensity: np.ndarray,
              indexes: Optional[np.ndarray],
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
    _input_validation(intensity, indexes, nbins, overlap)
    lengths_arr = lengths.astype(np.int64) if lengths is not None else np.array([intensity.size], dtype=np.int64)

    method_code = _method_to_code(method)

    noise_flat = estimate_flat_numba(
        intensity=intensity,
        indexes=indexes,
        lengths=lengths_arr,
        nbins=nbins,
        dynamic=dynamic,
        overlap=overlap,
        method_code=method_code,
        denoise_method=denoise_method,
        floor_value=0.001,
    )

    return noise_flat
