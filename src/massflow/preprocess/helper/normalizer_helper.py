import numpy as np
from massflow.tools.logger import get_logger
from typing import Optional
from massflow.preprocess.numba.normalization_numba import normalizer_numba

logger = get_logger("preprocesss")


def _input_validation(
    intensity:np.ndarray,
    index: Optional[np.ndarray] = None):
    """
    Validate input parameters for normalization functions.
    
    Parameters:
        intensity (np.ndarray): 1D intensity array to be preprocessed.
        index (Optional[np.ndarray]): 1D index array (e.g., m/z values).

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

def tic_normalize(
        intensity: np.ndarray,
        scale_method: str = 'none',
        scale: float = 1.0
):
    """
    Total Ion Current (TIC) normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale_method (str): Optional scaling after normalization:
            - 'none': no additional scaling
            - 'unit': min-max scale to [0, 1]
        scale (float): Cardinal-like amplitude scaling factor applied after normalization.

    Returns:
        np.ndarray: TIC-normalized intensity array. Sum equals `scale` (if scale_method='none').

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
    norm_intensity = norm_intensity * float(scale)
    # Apply scaling method
    norm_intensity = apply_scaling(norm_intensity, scale_method)

    return norm_intensity

def rms_normalize(
        intensity: np.ndarray,
        scale_method: str = 'none',
        scale: float = 1.0
):
    """
    Root Mean Square (RMS) normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale_method (str): Optional scaling after normalization:
            - 'none': no additional scaling
            - 'unit': min-max scale to [0, 1]
        scale (float): Cardinal-like amplitude scaling factor applied after normalization.

    Returns:
        np.ndarray: RMS-normalized intensity array. RMS equals `scale` (if scale_method='none').

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
    norm_intensity = norm_intensity * float(scale)
    # Apply optional scaling
    norm_intensity = apply_scaling(norm_intensity, scale_method)
    return norm_intensity

def median_normalize(
        intensity: np.ndarray,
        scale_method: str = 'none',
        scale: float = 1.0
):
    """
    Median normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale_method (str): Optional scaling after normalization:
            - 'none': no additional scaling
            - 'unit': min-max scale to [0, 1]
        scale (float): Cardinal-like amplitude scaling factor applied after normalization.

    Returns:
        np.ndarray: Median-normalized intensity array. Median equals `scale` (if scale_method='none').

    Raises:
        ValueError: If median of intensity is not greater than 0.
    """
    med = float(np.median(intensity))
    # Apply median normalization
    if med > 0.0:
        norm_intensity = intensity / med
    else:
        logger.error("Median value is not greater than 0, cannot normalize data")
        raise ValueError("Median value is not greater than 0, cannot normalize data")
    # Apply amplitude scaling (cardinal-style)
    norm_intensity = norm_intensity * float(scale)
    # Apply scaling method
    norm_intensity = apply_scaling(norm_intensity, scale_method)

    return norm_intensity

def apply_scaling(
        intensity: np.ndarray,
        scale_method: str
):
    """
    Apply scaling transformation to intensity data.

    Parameters:
        intensity (np.ndarray): Input intensity array after primary normalization.
        scale_method (str): Scaling method to apply:
                          - 'none': No additional scaling
                          - 'unit': Scale to 0-1 range using min-max normalization

    Returns:
        np.ndarray: Scaled intensity array.

    Raises:
        ValueError: If scale_method is not supported.
    """
    method_norm = (scale_method or 'none').strip().lower()
    if method_norm == 'none':
        return intensity
    elif method_norm == 'unit':
        # Min-max scaling to [0, 1] range
        if intensity.size == 0:
            return intensity
        intensity_min = np.min(intensity)
        intensity_max = np.max(intensity)
        if intensity_max - intensity_min > 0:
            return (intensity - intensity_min) / (intensity_max - intensity_min)
        else:
        # If all values are the same, return original values
            return intensity

    else:
        logger.error(f"Unsupported scale_method: {method_norm}. Supported methods are: 'none', 'unit'")
        raise ValueError(f"Unsupported scale_method: {method_norm}. Supported methods are: 'none', 'unit'")

def normalizer(
    intensity: np.ndarray,
    scale_method: str = 'none',
    method: str = "tic",
    scale: float = 1.0,
    numba_max_threads: Optional[int] = None,
    lengths: Optional[np.ndarray] = None,
):
    """
    Unified normalization dispatcher.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale_method (str): Optional scaling after normalization:
            - 'none': no additional scaling
            - 'unit': min-max scale to [0, 1]
        method (str): Normalization method.
            - Python backend: 'tic', 'rms', 'median'
            - Numba backend: 'tic_numba', 'rms_numba', 'median_numba'
        scale (float): Cardinal-like amplitude scaling factor (applied after normalization).
        numba_max_threads (Optional[int]): Thread cap when using the numba backend.
        lengths: Optional[np.ndarray]: Optional array of lengths for normalization.
    Returns:
        np.ndarray: Normalized (and optionally scaled) intensity array.

    Raises:
        ValueError: If `method` is not supported.
    """
    # Basic validation
    if not np.isfinite(scale) or float(scale) < 0.0:
        logger.error("scale must be a finite non-negative number")
        raise ValueError("scale must be a finite non-negative number")

    # Normalize and validate method
    method_norm = (method or "tic").strip().lower()
    scale_method_norm = (scale_method or 'none').strip().lower()
    
    # Check if Numba backend is requested via suffix
    is_numba = method_norm.endswith("_numba")
    base_method = method_norm.replace("_numba", "") if is_numba else method_norm
    
    supported = {"tic", "rms", "median"}
    if base_method not in supported:
        logger.error("Unsupported normalization method: %s. Use one of: tic, rms, median (suffix _numba for acceleration)", method)
        raise ValueError(f"Unsupported normalization method: {method}")

    if is_numba:
        return normalizer_numba(
            intensity,
            method=base_method,
            scale=scale,
            scale_method=scale_method_norm,
            numba_max_threads=numba_max_threads,
            lengths=lengths,
        )
    else:
        _input_validation(intensity)
        # Python Backend Dispatch
        if base_method == "tic":
            return tic_normalize(intensity, scale_method=scale_method_norm, scale=scale)
        elif base_method == "rms":
            return rms_normalize(intensity, scale_method=scale_method_norm, scale=scale)
        else:  # base_method == "median"
            return median_normalize(intensity, scale_method=scale_method_norm, scale=scale)
