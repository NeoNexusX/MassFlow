import numpy as np
from typing import Optional
from logger import get_logger

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

def tic_normalize(
        intensity: np.ndarray,
        scale_method: str = 'none'
):
    """
    Total Ion Current (TIC) normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale_method (str): Optional scaling after normalization:
            - 'none': no additional scaling
            - 'unit': min-max scale to [0, 1]

    Returns:
        np.ndarray: TIC-normalized intensity array. Sum equals 1 when input sum > 0.

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
    # Apply scaling method
    norm_intensity = apply_scaling(norm_intensity, scale_method)
    
    return norm_intensity

def median_normalize(
        intensity: np.ndarray,
        scale_method: str = 'none'
):
    """
    Median normalization.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale_method (str): Optional scaling after normalization:
            - 'none': no additional scaling
            - 'unit': min-max scale to [0, 1]

    Returns:
        np.ndarray: Median-normalized intensity array. Median equals 1 when input median > 0.

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
    if scale_method == 'none':
        return intensity
    elif scale_method == 'unit':
        # Min-max scaling to [0, 1] range
        if intensity.size == 0:
            return intensity
        intensity_min = np.min(intensity)
        intensity_max = np.max(intensity)
        if intensity_max - intensity_min > 0:
            return (intensity - intensity_min) / (intensity_max - intensity_min)
        else:
            # If all values are the same, keep original
            return np.zeros_like(intensity)
        
    else:
        logger.error(f"Unsupported scale_method: {scale_method}. "
                        f"Supported methods are: 'none', 'unit'")
        raise ValueError(f"Unsupported scale_method: {scale_method}. "
                        f"Supported methods are: 'none', 'unit'")

def normalizer(intensity: np.ndarray,
               scale_method: str = 'none',
               method: str = "tic"):
    """
    Unified normalization dispatcher.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array.
        scale_method (str): Optional scaling after normalization:
            - 'none': no additional scaling
            - 'unit': min-max scale to [0, 1]
        method (str): Primary normalization method:
            - 'tic': Total Ion Current normalization
            - 'median': Median normalization

    Returns:
        np.ndarray: Normalized (and optionally scaled) intensity array.

    Raises:
        ValueError: If `method` is not supported.
    """
    # # Basic validation
    _input_validation(intensity)
    # Choose normalization method
    if method == "tic":
        return tic_normalize(intensity, scale_method=scale_method)
    elif method == "median":
        return  median_normalize(intensity, scale_method=scale_method)
    else:
        supported = "tic, median"
        logger.error(f"Unsupported normalization method: {method}. Use one of: {supported}.")
        raise ValueError(f"Unknown normalization method: {method}")