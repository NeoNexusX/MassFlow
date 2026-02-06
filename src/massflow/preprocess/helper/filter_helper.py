"""
Author: MassFlow Development Team Bionet/NeoNexus lyk
License: See LICENSE file in project root
"""
from typing import Optional
import numpy as np
import pywt
from scipy.spatial import cKDTree #type: ignore
from scipy.signal import savgol_filter
from scipy.stats import norm
from scipy import stats
from scipy import signal, linalg
from massflow.tools.logger import get_logger
from massflow.module.spectrum import Spectrum
from massflow.preprocess.numba.noise_reduction_numba import (
    smooth_signal_savgol_numba,
    smooth_ns_signal_ma_numba,
    smooth_ns_signal_gaussian_numba,
    smooth_ns_signal_bi_numba,
    smooth_signal_ma_numba,
    smooth_signal_gaussian_numba,
    smooth_signal_ma_loop
)

logger = get_logger("preprocesss")

def _input_validation(
    intensity: np.ndarray,
    index: Optional[np.ndarray] = None,
    *,
    window: Optional[int] = None,
    sd: Optional[float] = None,
    k: Optional[int] = None,
    p: Optional[int] = None,
) -> None:
    """
    Unified validation and normalization for smoothing-related parameters.

    Performs validation only; does not normalize or mutate values.

    Raises:
        TypeError: `intensity` is not a NumPy array.
        ValueError: `intensity` is not 1D or empty; invalid `window/k/p/sd/...`;
            or `threshold_mode` is not 'soft'/'hard'.
    """
    # Validate intensity (must be non-empty 1D NumPy array)
    if not isinstance(intensity, np.ndarray):
        logger.error("intensity must be a numpy array")
        raise TypeError("intensity must be a numpy array")
    if intensity.ndim != 1 or intensity.size == 0:
        logger.error("intensity must be a non-empty 1D array")
        raise ValueError("intensity must be a non-empty 1D array")

    # Index check (if provided, must be 1D and match intensity length)
    if index is not None:
        if not isinstance(index, np.ndarray) or index.ndim != 1 or index.size != intensity.size:
            logger.error("index must be a 1D array with the same length as intensity")
            raise ValueError("index must be a 1D array with the same length as intensity")

    # Validate window (if provided): positive integer
    if window is not None:
        if not isinstance(window, int) or window <= 0:
            logger.error("window must be a positive integer")
            raise ValueError("window must be a positive integer")

    # Validate sd (if provided): positive float
    if sd is not None:
        if sd <= 0:
            logger.error("sd must be positive")
            raise ValueError("sd must be positive")

    # Validate k/p (if provided): positive integers; k must not exceed length
    if k is not None:
        if not isinstance(k, int) or k < 1 or k > intensity.size:
            logger.error("k must be a positive integer not exceeding intensity length")
            raise ValueError("k must be a positive integer not exceeding intensity length")
    if p is not None:
        if not isinstance(p, int) or p < 1:
            logger.error("p must be a positive integer")
            raise ValueError("p must be a positive integer")

    # No normalization; return None on success
    return None

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
        window (int): Window size used when `coef` is None. Must be a positive
            integer.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as the input.

    Raises:
        ValueError: `window` must be a positive integer when `coef` is None.
        TypeError: `intensity` must be 1D NumPy array; `coef` if provided must be 1D.
    """
    # Validate and prepare kernel
    if coef is None:
        if not isinstance(window, int) or window <= 0:
            logger.error("window must be a positive integer when coef is None")
            raise ValueError("window must be a positive integer when coef is None")
        # Force odd window length for centered alignment
        window = window + 1 - (window % 2)
        coef = np.ones(window, dtype=float)
    else:
        if not isinstance(coef, np.ndarray) or coef.ndim != 1 or coef.size == 0:
            logger.error("coef must be a non-empty 1D numpy array")
            raise TypeError("coef must be a non-empty 1D numpy array")

    # Normalize kernel weights
    coef = coef.astype(float)
    coef = coef / np.sum(coef)
    window = len(coef) #type: ignore
    half_window = window // 2

    # Boundary padding: extend using edge values
    xpad = np.pad(intensity, (half_window, half_window), mode="edge")

    # Convolution filtering
    y = np.convolve(xpad, coef, mode="valid") #type: ignore

    # Store as float32 to reduce memory footprint
    return y.astype(np.float32)

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
        window (int): Kernel window size (number of samples). Must be a
            positive integer.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as the input.

    Raises:
        ValueError: `window` must be a positive integer; `sd` must be positive if provided.
        TypeError: `intensity` must be a 1D array.
    """
    # Validate window and sd
    if not isinstance(window, int) or window <= 0:
        logger.error("window must be a positive integer")
        raise ValueError("window must be a positive integer")
    if sd is not None and sd <= 0:
        logger.error("sd must be positive")
        raise ValueError("sd must be positive")

    # Force odd window length for centered alignment
    window = window + 1 - (window % 2)

    half_window = window // 2

    # Generate Gaussian kernel
    if sd is None:
        sd = window / 4.0

    # Create Gaussian weights centered at zero
    positions = np.arange(-half_window, half_window + 1, dtype=float)

    coef = norm.pdf(positions, scale=sd)

    return smooth_signal_ma(intensity, coef=coef)

def smooth_signal_savgol(
        intensity:np.ndarray,
        window: int = 5,
        polyorder: int = 3,
        deriv: int = 0,
        delta: float = 1.0
):
    """
    Savitzky-Golay filter for signal smoothing.
    
    The Savitzky-Golay filter fits successive sub-sets of adjacent data points
    with a low-degree polynomial by linear least squares.
    
    Args:
        intensity : np.ndarray
            Input signal
        window : int
            Window size (must be a positive integer and greater than polyorder; SciPy requires odd window length).
        polyorder : int
            Polynomial order (must be less than window)
        deriv : int
            The order of the derivative to compute. This must be a nonnegative integer. The default is 0, which means to filter
            the data without differentiating.
        delta : float
            The spacing of the samples to which the filter will be applied.
            This is only used if deriv > 0. Default is 1.0.    
    Returns:
        np.ndarray
            Smoothed signal
    """

    if window < 3:
        window = 3

    # Ensure window is odd
    if window % 2 == 0:
        window += 1

    # Ensure polyorder is less than window
    if polyorder >= window:
        logger.error("polyorder must be less than window")
        raise ValueError("polyorder must be less than window")

    # Apply Savitzky-Golay filter
    return savgol_filter(intensity, window, polyorder, deriv=deriv, delta=delta)

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

    # Validate threshold mode
    if threshold_mode not in {"soft", "hard"}:
        logger.error("threshold_mode must be 'soft' or 'hard' for wavelet denoising")
        raise ValueError("threshold_mode must be 'soft' or 'hard' for wavelet denoising")

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
    intensity: np.ndarray,
    index: np.ndarray,
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
        ValueError: If `index` is None or length mismatches `intensity`.
    """

    # Enforce positivity
    if k < 1 or p < 1:
        logger.error("k and p must be positive integers")
        raise ValueError("k and p must be positive integers")

    if index is None or len(index) != len(intensity):
        logger.warning(
            "mz_list must be provided and match intensity length; using np.arange as fallback mz_list"
        )
        index = np.arange(len(intensity))
    else:
        logger.debug(f"using provided mz_list for neighborhood search with k={k} and p={p}\r\n{index[:5]}")

    if not isinstance(k, int):
        logger.error("k must be an integer")
        raise ValueError("k must be an integer")
    if not isinstance(p, int):
        logger.error("p must be an integer")
        raise ValueError("p must be an integer")

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
    smoothed_intensity = neigh_intensity.mean(axis=1, dtype=np.float32)

    return smoothed_intensity

def smooth_ns_signal_gaussian(
    intensity:np.ndarray,
    index:np.ndarray,
    k: int = 5,
    p: int = 1,
    sd: Optional[float] = None,
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
    if sd <= 0:
        logger.error("sd must be positive")
        raise ValueError("sd must be positive")

    # dist_ = np.exp(-dists**2 / (2 * sd**2))
    exponent = -0.5 * (dists / sd) ** 2
    exponent = np.clip(exponent, -88.0, 0.0)  # Numerical stability: clamp exponent to avoid underflow
    weights = np.exp(exponent)
    # calculate row-wise normalized weights avoid divide by zero
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0.0, 1.0, row_sums)
    weights = weights / row_sums

    smoothed_intensity = smooth_ns_signal_calculate(neigh_intensity, weights, axis=1)
    smoothed_intensity = smoothed_intensity.astype(np.float32)

    return smoothed_intensity

def smooth_ns_signal_bi(
    intensity:np.ndarray,
    index:np.ndarray,
    k: int = 5,
    p: int = 2,
    sd_dist: Optional[float] = None,
    sd_intensity: Optional[float] = None
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
    if sd_dist <= 0:
        logger.error("sd_dist must be positive")
        raise ValueError("sd_dist must be positive")

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
    sd_intensity = stats.median_abs_deviation(intensity, nan_policy="omit", scale="normal") if sd_intensity is None else sd_intensity #type: ignore
    if sd_intensity <= 0:   #type: ignore
        logger.error("sd_intensity must be positive")
        raise ValueError("sd_intensity must be positive")

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
    smoothed_intensity = smoothed_intensity.astype(np.float32)

    return smoothed_intensity

def smooth_preprocess(data:Spectrum):
    """
    Basic preprocessing pipeline for MS data smoothing.

    Parameters:
        data (SpectrumBaseModule): Spectrum object whose `intensity` will be sanitized.

    Returns:
        SpectrumBaseModule: The same object with non-negative intensity values.
    """
    intensity = data.intensity.copy()   #type: ignore
    data.intensity = None # clear intensity to avoid confusion
    intensity[intensity < 0] = 0

    data.intensity = intensity
    return data

#不对，这里没有库没，自己写的算到什么时候去了？？
def smooth_signal_loess(intensity: np.ndarray, window: int = 11):
    """
    Loess smoothing (quadratic) with tri-cubic weighting only.

    Parameters:
        intensity (np.ndarray): 1D intensity array to be smoothed.
        window (int): Local regression window size (odd, >=3). Auto-adjusted if too small/large.

    Returns:
        np.ndarray: Smoothed intensity array with the same length as input.
    """
    _input_validation(intensity=intensity, window=window)

    if window < 3:
        window = 3
    if window % 2 == 0:
        window += 1

    n = len(intensity)
    if window > n:
        adj = n if (n % 2 == 1) else max(3, n - 1)
        window = adj

    y = np.asarray(intensity, dtype=np.float64)
    leny = y.size

    halfw = np.floor((window / 2.)).astype(int)
    window = int(2. * halfw + 1.)
    x1 = np.arange(1. - halfw, (halfw - 1.) + 1)

    weight = (1. - np.divide(np.abs(x1), halfw) ** 3.) ** 1.5
    v = (np.vstack((np.hstack(weight), np.hstack(weight * x1), np.hstack(weight * x1 * x1)))).transpose()
    q, _ = linalg.qr(v, mode='economic')  #type: ignore

    alpha = np.dot(q[halfw - 1,], q.transpose())  #type: ignore
    yhat = signal.lfilter(alpha * weight, 1, y)
    yhat[int(halfw + 1) - 1:-halfw] = yhat[int(window - 1) - 1:-1]  #type: ignore

    x1 = np.arange(1., (window - 1.) + 1)
    v = (np.vstack((np.hstack(np.ones([1, window - 1])), np.hstack(x1), np.hstack(x1 * x1)))).transpose() #type: ignore

    for j in np.arange(1, (halfw) + 1):
        weight = (1. - np.divide(np.abs((np.arange(1, window) - j)), window - j) ** 3.) ** 1.5
        w = (np.kron(np.ones((3, 1)), weight)).transpose()
        q, _ = linalg.qr(v * w, mode='economic')  #type: ignore

        alpha = np.dot(q[j - 1,], q.transpose())  #type: ignore
        alpha = alpha * weight
        yhat[int(j) - 1] = np.dot(alpha, y[:int(window) - 1])  #type: ignore
        yhat[int(-j)] = np.dot(alpha, y[np.arange(leny - 1, leny - window, -1, dtype=int)])  #type: ignore

    return yhat

def smoother(intensity:np.ndarray,
            index:Optional[np.ndarray]=None,
            method: str = "ma",
            window: int = 5,
            sd: Optional[float] = None,
            sd_intensity: Optional[float] = None,
            p: int = 2,
            coef: Optional[np.ndarray] = None,
            polyorder: int = 2,
            deriv: int = 0,
            delta: float = 1.0,
            wavelet: str = 'db4',
            threshold_mode: str = 'soft',
            lengths: Optional[np.ndarray] = None,
            numba_max_threads: Optional[int] = None):
    """
    Unified smoothing entry for multiple methods.

    Parameters:
        intensity (np.ndarray): 1D intensity array.
        index (Optional[np.ndarray]): 1D coordinate array aligned with `intensity` for NS methods.
        method (str): One of {'ma','ma_numba','ma_loop','gaussian','gaussian_numba','savgol','savgol_numba','wavelet',
        'ma_ns','ma_ns_numba','gaussian_ns','gaussian_ns_numba','bi_ns','bi_ns_numba'}.
        window (int): Window size or neighbor count depending on method.
        sd (float, optional): Gaussian scale parameter for relevant methods.
        sd_intensity (float, optional): Intensity scale for bilateral method.
        p (int): Minkowski metric for NS queries.
        coef (np.ndarray, optional): Custom kernel for moving-average.
        polyorder (int): Polynomial order for Savitzky-Golay.
        deriv (int): The order of the derivative to compute. This must be a nonnegative integer. The default is 0, which means to filter
            the data without differentiating.
        delta (float): The spacing of the samples to which the filter will be applied.
            This is only used if deriv > 0. Default is 1.0.
        wavelet (str): Wavelet family for wavelet denoising.
        threshold_mode (str): 'soft' or 'hard' for wavelet thresholding.
        lengths (Optional[np.ndarray]): Array of valid lengths for each spectrum in a 2D batch.
        numba_max_threads (Optional[int]): Maximum number of threads for Numba parallel execution.

    Returns:
        np.ndarray: Smoothed intensity array.

    Raises:
        ValueError: If `method` is unsupported or parameter combinations are invalid.
        TypeError: If `intensity` is invalid or `coef` is not 1D.
    """
    # Normalize method and validate supported set
    method_norm = (method or "ma").strip().lower()

    supported = {
        "ma", "gaussian", "savgol", "wavelet",
        "ma_ns", "gaussian_ns", "bi_ns",
        "savgol_numba", "ma_ns_numba", "gaussian_ns_numba", "bi_ns_numba",
        "ma_numba", "gaussian_numba",
        "ma_loop",
    }
    if method_norm not in supported:
        logger.error(
            "Unsupported smoothing method: %s. Use one of: ma, gaussian, savgol, wavelet, ma_ns, gaussian_ns, bi_ns, "
            "savgol_numba, ma_ns_numba, gaussian_ns_numba, bi_ns_numba, ma_numba, gaussian_numba, ma_loop",
            method,
        )
        raise ValueError(f"Unsupported smoothing method: {method}")

    # Validate shared parameters only (no normalization)
    if method_norm == "ma":
        _input_validation(intensity=intensity, window=window if coef is None else None)
        return smooth_signal_ma(intensity, window=window, coef=coef)

    elif method_norm == "gaussian":
        _input_validation(intensity=intensity, window=window, sd=sd)
        return smooth_signal_gaussian(intensity, window=window, sd=sd)

    elif method_norm == "savgol":
        _input_validation(intensity=intensity, window=window)
        return smooth_signal_savgol(intensity, window=window, polyorder=polyorder, deriv=deriv, delta=delta)

    elif method_norm == "savgol_numba":
        # _input_validation(intensity=intensity, window=window)
        # Note: savgol_numba supports 2D batch inputs (n_spectra, n_points); the 1D-only input validation cannot be used here.
        return smooth_signal_savgol_numba(intensity, window=window, polyorder=polyorder, deriv=deriv, delta=delta, lengths=lengths, numba_max_threads=numba_max_threads)

    elif method_norm == "ma_numba":
        # Supports 2D batch inputs
        return smooth_signal_ma_numba(intensity, window=window, lengths=lengths, numba_max_threads=numba_max_threads)

    elif method_norm == "ma_loop":
        return smooth_signal_ma_loop(intensity, window=window, lengths=lengths, numba_max_threads=numba_max_threads)

    elif method_norm == "gaussian_numba":
        # Supports 2D batch inputs
        return smooth_signal_gaussian_numba(intensity, window=window, sd=sd, lengths=lengths, numba_max_threads=numba_max_threads)

    elif method_norm == "wavelet":
        _input_validation(intensity=intensity)
        return smooth_signal_wavelet(intensity, wavelet=wavelet, threshold_mode=threshold_mode)

    elif method_norm == "ma_ns":
        _input_validation(intensity=intensity, index=index, k=window, p=p)
        return smooth_ns_signal_ma(intensity, index=index, k=window, p=p)    #type: ignore

    elif method_norm == "gaussian_ns":
        _input_validation(intensity=intensity, index=index, k=window, p=p, sd=sd)
        return smooth_ns_signal_gaussian(intensity, index=index, k=window, p=p, sd=sd)  #type: ignore

    elif method_norm == "bi_ns":
        _input_validation(intensity=intensity, index=index, k=window, p=p)
        return smooth_ns_signal_bi(intensity, index=index, k=window, p=p, sd_dist=sd, sd_intensity=sd_intensity)  #type: ignore

    elif method_norm == "ma_ns_numba":
        _input_validation(intensity=intensity, index=index, k=window, p=p)
        return smooth_ns_signal_ma_numba(intensity, index=index, k=window, p=p, numba_max_threads=numba_max_threads)  #type: ignore

    elif method_norm == "gaussian_ns_numba":
        _input_validation(intensity=intensity, index=index, k=window, p=p, sd=sd)
        return smooth_ns_signal_gaussian_numba(intensity, index=index, k=window, p=p, sd=sd, numba_max_threads=numba_max_threads)  #type: ignore

    else:
        _input_validation(intensity=intensity, index=index, k=window, p=p)
        return smooth_ns_signal_bi_numba(intensity, index=index, k=window, p=p, sd_dist=sd, sd_intensity=sd_intensity, numba_max_threads=numba_max_threads)  #type: ignore
