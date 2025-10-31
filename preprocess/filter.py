from typing import Optional
from logger import get_logger
from module.ms_module import MSBaseModule
import numpy as np

logger = get_logger(__name__)


def smooth_signal_ma(
    x: MSBaseModule, coef: Optional[np.ndarray] = None, window: int = 5):
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

def smooth_signal_gaussian(
    x: MSBaseModule, sd: Optional[float] = None, window: int = 5):
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
    except ImportError:
        # Fallback to numpy-only implementation if SciPy is unavailable
        coef = np.exp(-0.5 * (positions / sd) ** 2)

    return smooth_signal_ma(x, coef=coef)

def smooth_ns_signal_pre(
    x: MSBaseModule,
    k: int = 5,
    p: int = 1,
):
    """
    Apply moving-average smoothing with neighborhood search and noise suppression.

    Parameters:
        x (MSBaseModule): Input MSI spectrum object. Must provide 1D `intensity`
            and optional `mz_list` aligned to the same length.
        k (int): Number of nearest neighbors used for smoothing. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.
        weights (Optional[np.ndarray]): Weight array with shape (k,), where k is the number of neighbors.
            If None, a uniform weight is used.

    Returns:
        np.ndarray: Smoothed intensity array with shape (N,), where N is the
            length of `intensity`.

    Raises:
        ValueError: If `k` or `p` is not a positive integer.
        TypeError: If `x.intensity` is not a 1D numpy array.
    """

    intensity = x.intensity
    mz_list = x.mz_list

    # Basic validation
    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        logger.error("x.intensity must be a 1D numpy array ")
        raise TypeError("x.intensity must be a 1D numpy array ")

    if k < 1 or p < 1:
        logger.error("k and p must be positive integers")
        raise ValueError("k and p must be positive integers")

    if mz_list is None or len(mz_list) != len(intensity):
        logger.warning(
            "mz_list must be provided and have the same length as intensity , using np list index as mz_list"
        )
        mz_list = np.arange(len(intensity))

    from scipy.spatial import cKDTree

    tree = cKDTree(mz_list.reshape(-1, 1))

    if len(intensity) < k:
        logger.warning("spectrum length must be greater than k")
        k = len(intensity)

    dists, idxs = tree.query(mz_list.reshape(-1, 1), k=k, p=p)

    # Ensure idxs is (N, k) for consistent neighbor aggregation
    if np.ndim(idxs) == 1:
        idxs = idxs.reshape(-1, 1)

    # Impute NaNs in intensity using mean of its k neighbors
    nan_idx = np.where(np.isnan(intensity))[0]
    if nan_idx.size > 0:
        # shape: (M, k)
        neigh = intensity[idxs[nan_idx]]
        fill_vals = np.nanmean(neigh, axis=1)
        # Fallback for rows whose neighbors are all NaN
        global_median = np.nanmedian(intensity)
        if np.isnan(global_median):
            global_median = 0.0
        fill_vals = np.where(np.isnan(fill_vals), global_median, fill_vals)
        intensity[nan_idx] = fill_vals

    neigh_intensity = intensity[idxs]  # shape: (N, k)

    return  neigh_intensity, dists, idxs


def smooth_ns_signal_calaulate(
    neigh_intensity: np.ndarray,
    weights: np.ndarray,
    axis: int = 1):
    """
    Calculate the smoothed intensity using the given neighborhood intensity and weights.

    Parameters:
        neigh_intensity (np.ndarray): Neighborhood intensity array with shape (N, k),
            where N is the number of data points and k is the number of neighbors.
        weights (np.ndarray): Weight array with shape (k,), where k is the number of neighbors.
        axis (int, optional): Axis along which to perform the sum. Default is 1.

    Returns:
        np.ndarray: Smoothed intensity array with shape (N,), where N is the number of data points.
    """
    return np.sum(neigh_intensity * weights, axis=axis)


def smooth_ns_signal_ma(
    x: MSBaseModule,
    k: int = 5,
    p: int = 1,
):
    """
    Apply moving-average smoothing with neighborhood search and noise suppression.

    Parameters:
        x (MSBaseModule): Input MSI spectrum object. Must provide 1D `intensity`
            and optional `mz_list` aligned to the same length.
        k (int): Number of nearest neighbors used for smoothing. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.

    Returns:
        np.ndarray: Smoothed intensity array with shape (N,), where N is the
            length of `intensity`.

    Raises:
        ValueError: If `k` or `p` is not a positive integer.
        TypeError: If `x.intensity` is not a 1D numpy array.
    """
    neigh_intensity, _, _ = smooth_ns_signal_pre(x, k=k, p=p)
    weights = np.ones(k, dtype=float)
    weights = weights / np.sum(weights)
    smoothed_intensity = smooth_ns_signal_calaulate(neigh_intensity, weights, axis=1)

    return smoothed_intensity


def smooth_ns_signal_gauss(
    x: MSBaseModule,
    k: int = 5,
    p: int = 1,
    sd: float = None,
):
    """
    Apply Gaussian smoothing with neighborhood search and noise suppression.
    """

    neigh_intensity, weights = smooth_ns_signal_pre(x, k=k, p=p)
    
