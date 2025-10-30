from math import log
import numpy as np
from typing import Optional
from logger import get_logger
from module.ms_module import MSBaseModule

logger = get_logger(__name__)

def smooth_signal_ma(x: MSBaseModule, coef: Optional[np.ndarray] = None, window: int = 5) -> np.ndarray:
    """
    Apply moving-average (or arbitrary kernel) smoothing to an MSI spectrum.

    Parameters:
        x (MSBaseModule): Input MSI spectrum object containing `intensity`.
        coef (Optional[np.ndarray]): Convolution kernel coefficients. If None, a
            uniform window is used. When provided, the full length of `coef` is
            used as the window size.
        window (int): Window size for uniform moving average when `coef` is None.
            Must be a positive integer. If even, it is adjusted to the next odd
            integer to preserve center alignment.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as the input.

    Raises:
        ValueError: If `window` is not a positive integer when `coef` is None.
        TypeError: If `x` is not an `MSBaseModule` instance or `intensity` is not 1D.
    """
    # Basic validation
    intensity = x.intensity
    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        raise TypeError("x.intensity must be a 1D numpy array")

    # If weights not specified, use uniform weights (moving average)
    if coef is None:
        if not isinstance(window, int) or window <= 0:
            raise ValueError("window must be a positive integer when coef is None")
        # Ensure window length is odd
        window = window + 1 - window % 2
        coef = np.ones(window, dtype=float)

    # Normalize kernel weights
    coef = coef.astype(float)
    coef = coef / np.sum(coef)
    window = len(coef)
    half_window = window // 2

    # Boundary padding: extend using edge values
    xpad = np.pad(intensity, (half_window, half_window), mode="edge")

    # Convolution filtering
    y = np.convolve(xpad, coef, mode="valid")

    return y

def smooth_signal_gaussian(x: MSBaseModule, sd: Optional[float] = None, window: int = 5) -> np.ndarray:
    """
    Apply Gaussian smoothing to an MSI spectrum using a discrete Gaussian kernel.

    Parameters:
        x (MSBaseModule): Input MSI spectrum object containing `intensity`.
        sd (Optional[float]): Standard deviation of the Gaussian. If None, it is
            set to `window / 4.0` to give a reasonable spread.
        window (int): Kernel window size (number of samples). Should be a
            positive integer. If even, it is adjusted to the next odd integer to
            preserve center alignment.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as the input.

    Raises:
        ValueError: If `window` is not a positive integer.
        TypeError: If `x` is not an `MSBaseModule` instance or `intensity` is not 1D.
    """
    intensity = x.intensity
    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        raise TypeError("x.intensity must be a 1D numpy array")
    if not isinstance(window, int) or window <= 0:
        raise ValueError("window must be a positive integer")

    # Generate Gaussian kernel
    if sd is None:
        sd = window / 4.0

    # Ensure window length is odd for symmetric kernel
    window = window + 1 - window % 2
    half_window = window // 2

    # Create Gaussian weights centered at zero
    positions = np.arange(-half_window, half_window + 1, dtype=float)
    try:
        # Prefer scipy.stats.norm if available for numerical stability
        from scipy.stats import norm  # type: ignore
        coef = norm.pdf(positions, scale=sd)
    except Exception:
        # Fallback to numpy-only implementation if SciPy is unavailable
        coef = np.exp(-0.5 * (positions / sd) ** 2)

    return smooth_signal_ma(x, coef=coef)

def smooth_ns_signal_ma(
    x: MSBaseModule, 
    k: int= 5, 
    p: int= 1,
    metric: str= "mse",):
    """
    Apply moving-average (or arbitrary kernel) smoothing to an MSI spectrum with noise suppression.
    ns (NS): Non-uniform sampling of the spectrum.
    """

    intensity = x.intensity
    mz_list = x.mz_list

    if k < 1 or p < 1:
        logger.error("p and k must be a positive integer")
        raise ValueError("k must be a positive integer")
    
    if mz_list is None or len(mz_list) != len(intensity):
        logger.warning("mz_list must be provided and have the same length as intensity , using np list index as mz_list")
        mz_list = np.arange(len(intensity))
    
    weights = np.ones(k)
    weights = weights / np.sum(weights)

    from scipy.spatial import cKDTree
    tree = cKDTree(mz_list.reshape(-1, 1))

    if len(intensity) < k:
        logger.warning("spectrum length must be greater than k")
        k = len(intensity)
    
    dists, idxs = tree.query(mz_list.reshape(-1, 1), k=k, p=p)

    # Ensure idxs is (N, k)
    if np.ndim(idxs) == 1:
        idxs = idxs.reshape(-1, 1)

    neigh_vals = intensity[idxs]  # shape: (N, k)

    y = np.sum(neigh_vals * weights, axis=1)

    return








