from typing import Any, Callable, Optional
import numpy as np
from massflow.tools.logger import get_logger
from massflow.preprocess.numba.normalization_numba import normalizer_numba
from massflow.tools.funs import dispatch_with_supported_kwargs

logger = get_logger("preprocesss")


def _input_validation(
    intensity: np.ndarray,
    index: Optional[np.ndarray] = None,
    scale: Optional[float] = None,
    ref_tolerance: float = 0.1,
    method_norm: Optional[str] = None,
    mz_flat: Optional[np.ndarray] = None,
    ref: Optional[float] = None,
):
    """
    Validate input parameters for normalization functions.
    
    Parameters:
        intensity (np.ndarray): 1D intensity array to be preprocessed.
        index (Optional[np.ndarray]): 1D index array (e.g., m/z values).
        scale (Optional[float]): Scale factor to validate.
        ref_tolerance (float): Ref matching tolerance to validate.
        method_norm (Optional[str]): Normalized method name.
        mz_flat (Optional[np.ndarray]): Flat m/z array used by ref_numba.
        ref (Optional[float]): Reference m/z value used by ref_numba.

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

def tic_normalize(
        intensity: np.ndarray,
    scale: Optional[float] = None
):
    """
    Total Ion Current (TIC) normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale (Optional[float]): Amplitude scaling factor applied after normalization.
            If None, defaults to current spectrum length.

    Returns:
        np.ndarray: TIC-normalized intensity array. Sum equals `scale`.

    Raises:
        ValueError: If TIC (sum of intensity) is not greater than 0.
    """
    tic = float(np.sum(intensity))
    # Apply TIC normalization
    if tic > 0.0:
        norm_intensity = intensity / tic
    else:
        logger.error("TIC value is not greater than 0, cannot normalize data")
        raise ValueError("TIC value is not greater than 0, cannot normalize data")
    # Apply amplitude scaling (cardinal-style)
    scale_val = float(intensity.size) if scale is None else float(scale)
    norm_intensity = norm_intensity * scale_val
    return norm_intensity

def rms_normalize(
        intensity: np.ndarray,
    scale: Optional[float] = None
):
    """
    Root Mean Square (RMS) normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale (Optional[float]): Amplitude scaling factor applied after normalization.
            If None, defaults to current spectrum length.

    Returns:
        np.ndarray: RMS-normalized intensity array. RMS equals `scale`.

    Raises:
        ValueError: If RMS of intensity is not greater than 0.

    Notes:
        Matches the provided R implementation: b = sqrt(mean(x^2, na.rm=TRUE));
        if b > 0: y = scale * x / b; otherwise an error is raised.
        Here `na.rm=TRUE` is implemented via np.nanmean.
    """
    # Compute RMS with NaN ignored
    b = float(np.sqrt(np.nanmean(np.square(intensity))))
    if not np.isfinite(b) or b <= 0.0:
        logger.error("RMS value is not greater than 0, cannot normalize data")
        raise ValueError("RMS value is not greater than 0, cannot normalize data")
    norm_intensity = intensity / b
    # Apply amplitude scaling (cardinal-style)
    scale_val = float(intensity.size) if scale is None else float(scale)
    norm_intensity = norm_intensity * scale_val
    return norm_intensity

def normalizer(
    intensity: np.ndarray,
    method: str = "tic",
    scale: Optional[float] = None,
    mz_flat: Optional[np.ndarray] = None,
    ref: Optional[float] = None,
    ref_tolerance: float = 0.1,
    lengths: Optional[np.ndarray] = None,
):
    """
    Unified normalization dispatcher.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        method (str): Normalization method.
            - Python backend: 'tic', 'rms'
            - Numba backend: 'tic_numba', 'rms_numba', 'ref_numba'
        scale (Optional[float]): Amplitude scaling factor (applied after normalization).
            If None, default is per-spectrum length for tic/rms, and 1 for ref.
        mz_flat (Optional[np.ndarray]): Flat mz array used in ref_numba mode.
        ref (Optional[float]): Reference mz value used in ref_numba mode.
        ref_tolerance (float): Matching tolerance used in ref_numba mode.
        numba_max_threads (Optional[int]): Thread cap when using the numba backend.
        lengths: Optional[np.ndarray]: Optional array of lengths for normalization.
    Returns:
        np.ndarray: Normalized intensity array.

    Raises:
        ValueError: If `method` is not supported.
    """
    # Normalize and validate method
    method_norm = (method or "tic").strip().lower()

    method_map: dict[str, Callable[..., Any]] = {
        "tic": tic_normalize,
        "rms": rms_normalize,
        "tic_numba": normalizer_numba,
        "rms_numba": normalizer_numba,
        "ref_numba": normalizer_numba,
    }

    target_func = method_map.get(method_norm)
    if target_func is None:
        logger.error("Unsupported normalization method: %s. Use one of: tic, rms, ref_numba", method)
        raise ValueError(f"Unsupported normalization method: {method}")

    _input_validation(
        intensity,
        scale=scale,
        ref_tolerance=ref_tolerance,
        method_norm=method_norm,
        mz_flat=mz_flat,
        ref=ref,
    )

    base_method = method_norm.replace("_numba", "")
    return dispatch_with_supported_kwargs(
        target_func,
        intensity=intensity,
        method=base_method,
        scale=scale,
        mz_flat=mz_flat,
        ref=ref,
        ref_tolerance=ref_tolerance,
        lengths=lengths,
    )
