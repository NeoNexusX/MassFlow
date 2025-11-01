from typing import Optional
import numpy as np
from scipy import stats
from logger import get_logger
from module.ms_module import MSBaseModule

logger = get_logger(__name__)


def smooth_signal_ma(
    x: MSBaseModule,
    coef: Optional[np.ndarray] = None,
    window: int = 5,
):
    """
    Apply moving-average (or arbitrary kernel) smoothing to a 1D spectrum.

    Parameters:
        x (MSBaseModule): Spectrum with 1D `intensity` array and aligned `mz_list`.
        coef (Optional[np.ndarray]): 1D convolution kernel. If None, uses a
            uniform kernel with length `window`.
        window (int): Window size used when `coef` is None. Must be > 0. If even,
            it will be adjusted to the next odd number for center alignment.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as the input.

    Raises:
    ValueError: If `window` <= 0 when `coef` is None.
    TypeError: If `intensity` is not a 1D numpy array.
    """
    # Basic validation
    intensity = x.intensity
    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        logger.error("x.intensity must be a 1D numpy array ")
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
    x: MSBaseModule,
    sd: Optional[float] = None,
    window: int = 5,
):
    """
    Apply Gaussian smoothing to a 1D spectrum using a discrete Gaussian kernel.

    Parameters:
    x (MSBaseModule): Spectrum with 1D `intensity` array and aligned `mz_list`.
        sd (Optional[float]): Standard deviation of the Gaussian. If None, it is
            set to `window / 4.0` to give a reasonable spread.
        window (int): Kernel window size (number of samples). Should be a
            positive integer. If even, it is adjusted to the next odd integer to
            preserve center alignment.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as the input.

    Raises:
        ValueError: If `window` is not a positive integer.
    TypeError: If `intensity` is not 1D.
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
    Prepare neighborhood search (kNN on m/z) for NS-based smoothing.

    Parameters:
        x (MSBaseModule): Input MSI spectrum object. Must provide 1D `intensity`
            and optional `mz_list` aligned to the same length.
        k (int): Number of nearest neighbors used for smoothing. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.
            (Weights are applied in subsequent functions.)

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
            - neigh_intensity: shape (N, k) neighbor intensities
            - dists: shape (N, k) neighbor distances (Minkowski p)
            - idxs: shape (N, k) neighbor indices

    Raises:
        ValueError: If `k` or `p` is not a positive integer.
    TypeError: If `x.intensity` is not a 1D numpy array.
    """

    intensity = x.intensity
    mz_list = x.mz_list

    # Basic validation
    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        logger.error("x.intensity must be a 1D numpy array")
        raise TypeError("x.intensity must be a 1D numpy array")

    if k < 1 or p < 1:
        logger.error("k and p must be positive integers")
        raise ValueError("k and p must be positive integers")

    if mz_list is None or len(mz_list) != len(intensity):
        logger.warning(
            "mz_list must be provided and match intensity length; using np.arange as fallback mz_list"
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

def smooth_ns_signal_calculate(
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
    smoothed_intensity = smooth_ns_signal_calculate(neigh_intensity, weights, axis=1)

    return smoothed_intensity

def smooth_ns_signal_gaussian(
    x: MSBaseModule,
    k: int = 5,
    p: int = 1,
    sd: float = None,
):
    """
    Apply Gaussian smoothing with neighborhood search and noise suppression.
    """

    neigh_intensity, dists, _ = smooth_ns_signal_pre(x, k=k, p=p)

    if len(dists.shape)<2:
        dists = dists.reshape(-1,1)
    dists_max = np.max(dists, axis=1)

    sd = np.median(dists_max) / 2.0 if sd is None else sd

    # dist_ = np.exp(-dists**2 / (2 * sd**2))
    exponent = -0.5 * (dists / sd) ** 2
    exponent = np.clip(exponent, -88.0, 0.0)  # float64 下约为 exp(-700) ~ 5e-305
    weights = np.exp(exponent)

    # calculate row-wise normalized weights avoid divide by zero
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0.0, 1.0, row_sums)
    weights = weights / row_sums

    smoothed_intensity = smooth_ns_signal_calculate(neigh_intensity, weights, axis=1)

    return smoothed_intensity

def smooth_ns_signal_bi(
    x: MSBaseModule,
    k: int = 5,
    p: int = 2,
    sd_dist: float = None,
    sd_intensity: float = None
):
    """
    Apply Bilateral Gaussian smoothing with neighborhood search and noise suppression.
    """
    intensity = x.intensity
    neigh_intensity, dists, _ = smooth_ns_signal_pre(x, k=k, p=p)

    if len(dists.shape)<2:
        dists = dists.reshape(-1,1)
    dists_max = np.max(dists, axis=1)

    sd_dist = np.median(dists_max) / 2.0 if sd_dist is None else sd_dist

    # dist_ = np.exp(-dists**2 / (2 * sd**2))
    exponent = -0.5 * (dists / sd_dist) ** 2
    lower = np.log(np.finfo(exponent.dtype if hasattr(exponent, "dtype") else np.float64).tiny)
    exponent = np.clip(exponent, lower, 0.0)
    s_weights = np.exp(exponent)

    # calculate row-wise normalized weights avoid divide by zero
    srow_sums = s_weights.sum(axis=1, keepdims=True)
    srow_sums = np.where(srow_sums == 0.0, 1.0, srow_sums)
    s_weights = s_weights / srow_sums

    # Intensity weights (based on neighbor intensity differences)
    sd_intensity = stats.median_abs_deviation(intensity, nan_policy="omit", scale="normal") if sd_intensity is None else sd_intensity

    diff = neigh_intensity - intensity.reshape(-1, 1)  # shape (N, k)
    iexponent = -0.5 * (diff / sd_intensity) ** 2
    ilower = np.log(np.finfo(iexponent.dtype if hasattr(iexponent, "dtype") else np.float64).tiny)
    iexponent = np.clip(iexponent, ilower, 0.0)
    a_weights = np.exp(iexponent)

    # Multiply weights and normalize row-wise to avoid division by zero
    combined = s_weights * a_weights
    row_sums = combined.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0.0, 1.0, row_sums)
    weights = combined / row_sums

    smoothed_intensity = smooth_ns_signal_calculate(neigh_intensity, weights, axis=1)

    return smoothed_intensity


def smooth_preprocess(data:MSBaseModule):
    """ A general preprocess pipeline for MS data smoothing
    """
    pass