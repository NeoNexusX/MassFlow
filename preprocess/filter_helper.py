"""
Author: MassFlow Development Team Bionet/NeoNexus lyk
License: See LICENSE file in project root
"""
from typing import Optional
import numpy as np
from scipy import stats
from logger import get_logger
from module.ms_module import SpectrumBaseModule

logger = get_logger(__name__)

def _input_validation(
    intensity:np.ndarray,
    index: Optional[np.ndarray] = None):
    """
    Validate input parameters for smoothing functions.
    
    Parameters:
        intensity (np.ndarray): 1D intensity array to be preprocessed.
        index (Optional[np.ndarray]): 1D index array (e.g., m/z values). If None, 
            will be generated as np.arange(len(intensity)).
    
    Raises:
        TypeError: If intensity is not a numpy array.
        ValueError: If intensity is not 1D or empty.
    """
    # Validate intensity array
    if intensity is None:
        logger.error("intensity must be a numpy array")
        raise TypeError("intensity must be a numpy array")

    elif intensity.ndim != 1 or intensity.size == 0:
        logger.error("intensity must be a non-empty 1D array")
        raise ValueError("intensity must be a non-empty 1D array")

    if index is not None and (index.ndim != 1 or index.size != intensity.size):
        logger.error("index must be a 1D array with the same length as intensity")
        raise ValueError("index must be a 1D array with the same length as intensity")

def smooth_signal_ma(
    intensity:np.ndarray,
    coef: Optional[np.ndarray] = None,
    window: int = 5,
):
    """
    Apply moving-average (or arbitrary kernel) smoothing to a 1D spectrum.

    Parameters:
        intensity (np.ndarray): 1D intensity array to be preprocessed.
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
    intensity:np.ndarray,
    sd: Optional[float] = None,
    window: int = 5,
):
    """
    Gaussian smoothing for a 1D intensity array using a discrete kernel.

    Parameters:
        intensity (np.ndarray): 1D intensity array to be smoothed.
        sd (Optional[float]): Standard deviation of the Gaussian. If None, set to
            `window / 4.0` for a reasonable spread.
        window (int): Kernel window size (number of samples). If even, it is
            adjusted to the next odd integer to preserve center alignment.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as the input.

    Raises:
        ValueError: If `window` <= 0.
        TypeError: If `intensity` is not a 1D array.
    """

    # Ensure window length is odd for symmetric kernel
    window = window + 1 - window % 2
    half_window = window // 2

    # Generate Gaussian kernel
    if sd is None:
        sd = window / 4.0

    # Create Gaussian weights centered at zero
    positions = np.arange(-half_window, half_window + 1, dtype=float)
    try:
        # Prefer scipy.stats.norm if available for numerical stability
        from scipy.stats import norm  # type: ignore

        coef = norm.pdf(positions, scale=sd)
    except ImportError:
        # Fallback to numpy-only implementation if SciPy is unavailable
        coef = np.exp(-0.5 * (positions / sd) ** 2)

    return smooth_signal_ma(intensity, coef=coef)

def smooth_signal_savgol(
        intensity:np.ndarray,
        window: int = 5,
        polyorder: int = 2
):
    """
    Savitzky-Golay filter for signal smoothing.
    
    The Savitzky-Golay filter fits successive sub-sets of adjacent data points 
    with a low-degree polynomial by the method of linear least squares.
    
    Args:
        intensity : np.ndarray
            Input signal
        window : int
            Window size (must be odd and greater than polyorder)
        polyorder : int
            Polynomial order (must be less than window)
            
    Returns:
        np.ndarray
            Smoothed signal
    """
    from scipy.signal import savgol_filter

    # Minimum window size check
    if window < 3:
        window = 3

    # Ensure window is odd
    if window % 2 == 0:
        window += 1

    # Ensure polyorder is less than window
    if polyorder >= window:
        polyorder = window - 1

    # Apply Savitzky-Golay filter
    return savgol_filter(intensity, window, polyorder)

def smooth_signal_wavelet(
        intensity:np.ndarray,
        wavelet: str = 'db4',
        threshold_mode: str = 'soft'
):
    """
    Wavelet-based denoising for 1D intensity smoothing.

    Parameters:
        intensity (np.ndarray): 1D intensity array to be denoised.
        wavelet (str): Wavelet family (e.g., 'db4', 'haar').
        threshold_mode (str): Thresholding mode, 'soft' or 'hard'.

    Returns:
        np.ndarray: Denoised intensity array with the same length as input.

    Raises:
        ImportError: If `pywt` is not available.
    """
    import pywt
    original_length = len(intensity)

    # Perform wavelet decomposition
    coeffs = pywt.wavedec(intensity, wavelet, mode='symmetric')

    # Estimate noise standard deviation using the finest detail coefficients
    sigma = np.median(np.abs(coeffs[-1])) / 0.6745

    # Calculate threshold using Donoho-Johnstone threshold
    threshold = sigma * np.sqrt(2 * np.log(len(intensity)))

    # Apply thresholding to all detail coefficients
    coeffs_thresh = list(coeffs)
    coeffs_thresh[1:] = [pywt.threshold(detail, threshold, mode=threshold_mode)
                        for detail in coeffs[1:]]

    # Reconstruct the signal
    reconstructed = pywt.waverec(coeffs_thresh, wavelet, mode='symmetric')

    # Match output length exactly to input
    if len(reconstructed) != original_length:
        if len(reconstructed) > original_length:
            reconstructed = reconstructed[:original_length]
        else:
            pad_length = original_length - len(reconstructed)
            reconstructed = np.pad(reconstructed, (0, pad_length), mode='edge')

    return reconstructed

def smooth_ns_signal_pre(
    intensity:np.ndarray,
    index:np.ndarray,
    k: int = 5,
    p: int = 1,
):
    """
    Prepare neighborhood search (kNN over `index`) for NS-based smoothing.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity` (e.g., m/z).
        k (int): Number of nearest neighbors used for smoothing. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
            - neigh_intensity: shape (N, k) neighbor intensities
            - dists: shape (N, k) neighbor distances (Minkowski p)
            - idxs: shape (N, k) neighbor indices

    Raises:
        ValueError: If `k` or `p` is not positive.
        ValueError: If `index` is None or length mismatches `intensity` (fallback applied).
    """

    if k < 1 or p < 1:
        logger.error("k and p must be positive integers")
        raise ValueError("k and p must be positive integers")

    if index is None or len(index) != len(intensity):
        logger.warning(
            "mz_list must be provided and match intensity length; using np.arange as fallback mz_list"
        )
        index = np.arange(len(intensity))

    from scipy.spatial import cKDTree

    tree = cKDTree(index.reshape(-1, 1))

    if len(intensity) < k:
        logger.warning("spectrum length must be greater than k")
        k = len(intensity)

    dists, idxs = tree.query(index.reshape(-1, 1), k=k, p=p)

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
    axis: int = 1
):
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
    intensity:np.ndarray,
    index:np.ndarray,
    k: int = 5,
    p: int = 1,
):
    """
    Moving-average smoothing with neighborhood search and row-wise normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity`.
        k (int): Number of neighbors. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.

    Returns:
        np.ndarray: Smoothed intensity array of shape (N,).

    Raises:
        ValueError: If `k` or `p` is not positive.
    """
    neigh_intensity, _, _ = smooth_ns_signal_pre(intensity, index, k=k, p=p)
    weights = np.ones(k, dtype=float)
    weights = weights / np.sum(weights)
    smoothed_intensity = smooth_ns_signal_calculate(neigh_intensity, weights, axis=1)

    return smoothed_intensity

def smooth_ns_signal_gaussian(
    intensity:np.ndarray,
    index:np.ndarray,
    k: int = 5,
    p: int = 1,
    sd: float = None,
):
    """
    Gaussian-weighted neighborhood smoothing with distance-based weights.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity`.
        k (int): Number of neighbors. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.
        sd (float, optional): Gaussian scale over neighbor distances; auto-estimated if None.

    Returns:
        np.ndarray: Smoothed intensity array of shape (N,).

    Raises:
        ValueError: If `k` or `p` is not positive.
    """

    neigh_intensity, dists, _ = smooth_ns_signal_pre(intensity, index, k=k, p=p)

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
    intensity:np.ndarray,
    index:np.ndarray,
    k: int = 5,
    p: int = 2,
    sd_dist: float = None,
    sd_intensity: float = None
):
    """
    Bilateral Gaussian smoothing combining spatial and intensity similarity.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity`.
        k (int): Number of neighbors. Must be >= 1.
        p (int): Minkowski metric for KD-tree queries (distance).
        sd_dist (float, optional): Spatial Gaussian scale over neighbor distances.
        sd_intensity (float, optional): Intensity Gaussian scale; MAD-based if None.

    Returns:
        np.ndarray: Smoothed intensity array of shape (N,).

    Raises:
        ValueError: If `k` or `p` is not positive.
    """
    neigh_intensity, dists, _ = smooth_ns_signal_pre(intensity, index, k=k, p=p)

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

def smooth_preprocess(data:SpectrumBaseModule):
    """
    Basic preprocessing pipeline for MS data smoothing.

    Parameters:
        data (SpectrumBaseModule): Spectrum object whose `intensity` will be sanitized.

    Returns:
        SpectrumBaseModule: The same object with non-negative intensity values.
    """
    intensity = data.intensity.copy()
    data.intensity = None # clear intensity to avoid confusion
    intensity[intensity < 0] = 0

    data.intensity = intensity
    return data

def smoother(intensity:np.ndarray,
            index:Optional[np.ndarray]=None,
            method: str = "ma",
            window: int = 5,
            sd: float = None,
            sd_intensity: float = None,
            p: int = 2,
            coef: np.ndarray = None,
            polyorder: int = 2,
            wavelet: str = 'db4',
            threshold_mode: str = 'soft'):

    """
    Unified smoothing entry for multiple methods.

    Parameters:
        intensity (np.ndarray): 1D intensity array.
        index (Optional[np.ndarray]): 1D coordinate array aligned with `intensity` for NS methods.
        method (str): One of {'ma','gaussian','savgol','wavelet','ma_ns','gaussian_ns','bi_ns'}.
        window (int): Window size or neighbor count depending on method.
        sd (float, optional): Gaussian scale parameter for relevant methods.
        sd_intensity (float, optional): Intensity scale for bilateral method.
        p (int): Minkowski metric for NS queries.
        coef (np.ndarray, optional): Custom kernel for moving-average.
        polyorder (int): Polynomial order for Savitzky-Golay.
        wavelet (str): Wavelet family for wavelet denoising.
        threshold_mode (str): 'soft' or 'hard' for wavelet thresholding.

    Returns:
        np.ndarray: Smoothed intensity array.

    Raises:
        ValueError: If `method` is unsupported or invalid parameter combinations are provided.
    """
    # # Basic validation
    _input_validation(intensity)

    # Choose smoothing method
    if method == "ma":
        return  smooth_signal_ma(intensity, coef=coef, window=window)
    elif method == "gaussian":
        return  smooth_signal_gaussian(intensity, sd=sd, window=window)
    elif method == "ma_ns":
        return  smooth_ns_signal_ma(intensity,index=index, p=p, k=window)
    elif method == "savgol":
        return smooth_signal_savgol(intensity, window=window, polyorder=polyorder)
    elif method == "wavelet":
        return smooth_signal_wavelet(intensity, wavelet=wavelet, threshold_mode=threshold_mode)
    elif method == "gaussian_ns":
        return smooth_ns_signal_gaussian(intensity,index=index, sd=sd,p=p,k=window)
    elif method == "bi_ns":
        return smooth_ns_signal_bi(intensity,index=index, sd_dist=sd, sd_intensity=sd_intensity, p=p,k=window)
    else:
        supported = "ma, gaussian, savgol, wavelet, ma_ns, gaussian_ns, bi_ns"
        logger.error(f"Unsupported smoothing method: {method}. Use one of: {supported}.")
        raise ValueError(f"Unknown smoothing method: {method}")
