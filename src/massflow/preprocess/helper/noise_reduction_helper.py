"""
Author: MassFlow Development Team Bionet/NeoNexus lyk
License: See LICENSE file in project root
"""

from inspect import signature
from typing import Any, Callable, Optional
import numpy as np
import pywt
from scipy.signal import savgol_filter
from scipy.stats import norm
from scipy import stats
from scipy import signal as scipy_signal, linalg
from massflow.preprocess.numba.noise_reduction_numba import ns_signal_pre
from massflow.tools.logger import get_logger
from massflow.module import Spectrum
from massflow.preprocess.numba.noise_reduction_numba import (
    smooth_signal_savgol_numba,
    smooth_ns_signal_gaussian_numba,
    smooth_ns_signal_bi_numba,
    smooth_signal_gaussian_numba,
    smooth_signal_ma_numba,
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
    """
    # Validate intensity (must be non-empty 1D NumPy array)
    if (not isinstance(intensity, np.ndarray)
        or intensity.ndim != 1
        or intensity.size == 0
    ):
        logger.error("intensity must be a numpy array")
        raise TypeError("intensity must be a non-empty 1D numpy array")

    # Index check (if provided, must be 1D and match intensity length)
    if index is not None:
        if (not isinstance(index, np.ndarray)
            or index.ndim != 1
            or index.size != intensity.size
        ):
            logger.error("index must be a 1D array with the same length as intensity")
            raise ValueError( "index must be a 1D array with the same length as intensity")

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


def smooth_signal_ma(
    intensity: np.ndarray,
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
    window = len(coef)  # type: ignore
    half_window = window // 2

    # Boundary padding: extend using edge values
    xpad = np.pad(intensity, (half_window, half_window), mode="edge")

    # Convolution filtering
    y = np.convolve(xpad, coef, mode="valid")  # type: ignore

    # Store as float32 to reduce memory footprint
    return y.astype(np.float32)


def smooth_signal_gaussian(
    intensity: np.ndarray,
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
    intensity: np.ndarray,
    window: int = 5,
    polyorder: int = 3,
    deriv: int = 0,
    delta: float = 1.0,
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
    intensity: np.ndarray, wavelet: str = "db4", threshold_mode: str = "soft"
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
        raise ValueError(
            "threshold_mode must be 'soft' or 'hard' for wavelet denoising"
        )

    original_length = len(intensity)

    # Perform wavelet decomposition
    coeffs = pywt.wavedec(intensity, wavelet, mode="symmetric")

    # Estimate noise standard deviation using the finest detail coefficients
    sigma = np.median(np.abs(coeffs[-1])) / 0.6745

    # Calculate threshold using Donoho-Johnstone threshold
    threshold = sigma * np.sqrt(2 * np.log(len(intensity)))

    # Apply thresholding to all detail coefficients
    coeffs_thresh = list(coeffs)
    coeffs_thresh[1:] = [
        pywt.threshold(detail, threshold, mode=threshold_mode) for detail in coeffs[1:]
    ]

    # Reconstruct the signal
    reconstructed = pywt.waverec(coeffs_thresh, wavelet, mode="symmetric")

    # Match output length exactly to input
    if len(reconstructed) != original_length:
        if len(reconstructed) > original_length:
            reconstructed = reconstructed[:original_length]
        else:
            pad_length = original_length - len(reconstructed)
            reconstructed = np.pad(reconstructed, (0, pad_length), mode="edge")

    return reconstructed


def smooth_ns_signal_gaussian(
    intensity: np.ndarray,
    index: np.ndarray,
    window: int = 5,
    p: int = 1,
    sd: Optional[float] = None,
):
    """
    Gaussian-weighted neighborhood smoothing with distance-based weights.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity`.
        window (int): Number of neighbors. Must be >= 1.
        p (int): Minkowski metric parameter for KD-tree query. Must be >= 1.
        sd (float, optional): Gaussian scale over neighbor distances; auto-estimated if None.

    Returns:
        np.ndarray: Smoothed intensity array of shape (N,).

    Raises:
        ValueError: If `k` or `p` is not positive.
    """

    neigh_intensity, dists, _ = ns_signal_pre(intensity, index, k=window, p=p)

    if len(dists.shape) < 2:
        dists = dists.reshape(-1, 1)
    dists_max = np.max(dists, axis=1)

    sd = np.median(dists_max) / 2.0 if sd is None else sd
    if sd is None or sd <= 0:
        logger.error("sd must be positive")
        raise ValueError("sd must be positive")

    # dist_ = np.exp(-dists**2 / (2 * sd**2))
    exponent = -0.5 * (dists / sd) ** 2
    exponent = np.clip(
        exponent, -88.0, 0.0
    )  # Numerical stability: clamp exponent to avoid underflow
    weights = np.exp(exponent)
    # calculate row-wise normalized weights avoid divide by zero
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0.0, 1.0, row_sums)
    weights = weights / row_sums

    smoothed_intensity = np.sum(neigh_intensity * weights, axis=1)
    smoothed_intensity = smoothed_intensity.astype(np.float32)

    return smoothed_intensity


def smooth_ns_signal_bi(
    intensity: np.ndarray,
    index: np.ndarray,
    window: int = 5,
    p: int = 2,
    sd: Optional[float] = None,
    sd_intensity: Optional[float] = None,
):
    """
    Bilateral Gaussian smoothing combining spatial and intensity similarity.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity`.
        window (int): Number of neighbors. Must be >= 1.
        p (int): Minkowski metric for KD-tree queries (distance).
        sd_dist (float, optional): Spatial Gaussian scale over neighbor distances.
        sd_intensity (float, optional): Intensity Gaussian scale; MAD-based if None.

    Returns:
        np.ndarray: Smoothed intensity array of shape (N,).

    Raises:
        ValueError: If `k` or `p` is not positive.
    """
    neigh_intensity, dists, _ = ns_signal_pre(intensity, index, k=window, p=p)

    if len(dists.shape) < 2:
        dists = dists.reshape(-1, 1)
    dists_max = np.max(dists, axis=1)

    sd = np.median(dists_max) / 2.0 if sd is None else sd
    if sd is None or sd <= 0:
        logger.error("sd_dist must be positive")
        raise ValueError("sd_dist must be positive")

    # dist_ = np.exp(-dists**2 / (2 * sd**2))
    exponent = -0.5 * (dists / sd) ** 2
    lower = np.log(
        np.finfo(exponent.dtype if hasattr(exponent, "dtype") else np.float64).tiny
    )
    exponent = np.clip(exponent, lower, 0.0)
    s_weights = np.exp(exponent)

    # calculate row-wise normalized weights avoid divide by zero
    srow_sums = s_weights.sum(axis=1, keepdims=True)
    srow_sums = np.where(srow_sums == 0.0, 1.0, srow_sums)
    s_weights = s_weights / srow_sums

    # Intensity weights (based on neighbor intensity differences)
    sd_intensity = stats.median_abs_deviation(intensity, nan_policy="omit", scale="normal") if sd_intensity is None else sd_intensity  # type: ignore
    if sd_intensity <= 0:  # type: ignore
        logger.error("sd_intensity must be positive")
        raise ValueError("sd_intensity must be positive")

    diff = neigh_intensity - intensity.reshape(-1, 1)  # shape (N, k)
    iexponent = -0.5 * (diff / sd_intensity) ** 2
    ilower = np.log(
        np.finfo(iexponent.dtype if hasattr(iexponent, "dtype") else np.float64).tiny
    )
    iexponent = np.clip(iexponent, ilower, 0.0)
    a_weights = np.exp(iexponent)

    # Multiply weights and normalize row-wise to avoid division by zero
    combined = s_weights * a_weights
    row_sums = combined.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0.0, 1.0, row_sums)
    weights = combined / row_sums

    smoothed_intensity = np.sum(neigh_intensity * weights, axis=1)
    smoothed_intensity = smoothed_intensity.astype(np.float32)

    return smoothed_intensity


def smooth_preprocess(data: Spectrum):
    """
    Basic preprocessing pipeline for MS data smoothing.

    Parameters:
        data (SpectrumBaseModule): Spectrum object whose `intensity` will be sanitized.

    Returns:
        SpectrumBaseModule: The same object with non-negative intensity values.
    """
    intensity = data.intensity.copy()  # type: ignore
    data.intensity = None  # clear intensity to avoid confusion
    intensity[intensity < 0] = 0

    data.intensity = intensity
    return data


def smooth_signal_loess(intensity: np.ndarray, window: int = 5):
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

    halfw = np.floor((window / 2.0)).astype(int)
    window = int(2.0 * halfw + 1.0)
    x1 = np.arange(1.0 - halfw, (halfw - 1.0) + 1)

    weight = (1.0 - np.divide(np.abs(x1), halfw) ** 3.0) ** 1.5
    v = (
        np.vstack(
            (np.hstack(weight), np.hstack(weight * x1), np.hstack(weight * x1 * x1))
        )
    ).transpose()
    q, _ = linalg.qr(v, mode="economic")  # type: ignore

    alpha = np.dot(q[halfw - 1,], q.transpose())  # type: ignore
    yhat = scipy_signal.lfilter(alpha * weight, 1, y)
    yhat[int(halfw + 1) - 1 : -halfw] = yhat[int(window - 1) - 1 : -1]  # type: ignore

    x1 = np.arange(1.0, (window - 1.0) + 1)
    v = (np.vstack((np.hstack(np.ones([1, window - 1])), np.hstack(x1), np.hstack(x1 * x1)))).transpose()  # type: ignore

    for j in np.arange(1, (halfw) + 1):
        weight = (
            1.0 - np.divide(np.abs((np.arange(1, window) - j)), window - j) ** 3.0
        ) ** 1.5
        w = (np.kron(np.ones((3, 1)), weight)).transpose()
        q, _ = linalg.qr(v * w, mode="economic")  # type: ignore

        alpha = np.dot(q[j - 1,], q.transpose())  # type: ignore
        alpha = alpha * weight
        yhat[int(j) - 1] = np.dot(alpha, y[: int(window) - 1])  # type: ignore
        yhat[int(-j)] = np.dot(alpha, y[np.arange(leny - 1, leny - window, -1, dtype=int)])  # type: ignore

    return yhat


def smoother(
    intensity: Optional[np.ndarray] = None,
    index: Optional[np.ndarray] = None,
    *,
    flat: Optional[np.ndarray] = None,
    method: str = "ma",
    window: int = 5,
    sd: Optional[float] = None,
    sd_intensity: Optional[float] = None,
    p: int = 2,
    coef: Optional[np.ndarray] = None,
    polyorder: int = 2,
    deriv: int = 0,
    delta: float = 1.0,
    wavelet: str = "db4",
    threshold_mode: str = "soft",
    lengths: Optional[np.ndarray] = None,
):
    """
    Unified smoothing entry for multiple methods.

    Parameters:
        intensity (Optional[np.ndarray]): Legacy 1D intensity input.
        flat (Optional[np.ndarray]): Flat 1D intensity input. When provided, it is preferred over `intensity`.
        index (Optional[np.ndarray]): 1D coordinate array aligned with `intensity` for NS methods.
        method (str): One of {'ma','ma_numba','gaussian','gaussian_numba','savgol','savgol_numba','wavelet',
        'gaussian_ns','gaussian_ns_numba','bi_ns','bi_ns_numba'}.
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
        lengths (Optional[np.ndarray]): Array of valid lengths for each spectrum in flat mode.

    Returns:
        np.ndarray: Smoothed intensity array.

    Raises:
        ValueError: If `method` is unsupported or parameter combinations are invalid.
        TypeError: If input intensity is invalid or `coef` is not 1D.
    """
    def _dispatch_with_supported_kwargs(func: Callable[..., Any], **kwargs: Any) -> np.ndarray:
        supported = signature(func).parameters
        filtered_kwargs = {name: value for name, value in kwargs.items() if name in supported}
        result = func(**filtered_kwargs)
        return np.asarray(result)

    # Normalize method and validate supported set
    method_norm = (method or "ma").strip().lower()

    input_signal = flat if flat is not None else intensity
    if input_signal is None:
        raise TypeError("Either flat or intensity must be provided")

    _input_validation(input_signal, index, window=window, sd=sd, k=window if "ns" in method_norm else None, p=p)

    method_map: dict[str, Callable[..., Any]] = {
        "ma": smooth_signal_ma,
        "gaussian": smooth_signal_gaussian,
        "savgol": smooth_signal_savgol,
        "savgol_numba": smooth_signal_savgol_numba,
        "ma_numba": smooth_signal_ma_numba,
        "gaussian_numba": smooth_signal_gaussian_numba,
        "wavelet": smooth_signal_wavelet,
        "gaussian_ns": smooth_ns_signal_gaussian,
        "bi_ns": smooth_ns_signal_bi,
        "gaussian_ns_numba": smooth_ns_signal_gaussian_numba,
        "bi_ns_numba": smooth_ns_signal_bi_numba,
    }

    target_func = method_map.get(method_norm)
    if target_func is None:
        raise ValueError(f"Unsupported smoothing method: {method}")

    return _dispatch_with_supported_kwargs(
        target_func,
        intensity=input_signal,
        flat=input_signal,
        index=index,
        window=window,
        sd=sd,
        sd_intensity=sd_intensity,
        p=p,
        coef=coef,
        polyorder=polyorder,
        deriv=deriv,
        delta=delta,
        wavelet=wavelet,
        threshold_mode=threshold_mode,
        lengths=lengths,
    )
