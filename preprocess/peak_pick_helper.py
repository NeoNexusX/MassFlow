"""
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from scipy.signal import find_peaks
import numpy as np
from logger import get_logger
from .est_noise_helper import estimator
from typing import List
from scipy.integrate import simpson
logger = get_logger(__name__)


def _input_validation(intensity,index):
    """Input validation for peak pick function."""
    if intensity.ndim != 1:
        logger.error("Input data must be a 1D array")
        raise ValueError("Input data must be a 1D array")
    if index.ndim != 1:
        logger.error("Index data must be a 1D array")
        raise ValueError("Index data must be a 1D array")


def _mask_peak_props(props: dict, mask: np.ndarray) -> dict:
    """
    Apply a boolean mask to all 1D numpy arrays in the `props` dict returned by `find_peaks`.

    Args:
        props (Dict[str, np.ndarray]): Dictionary of peak properties aligned to `peaks`.
        mask (np.ndarray): 1D boolean mask of length equal to `len(peaks)` before filtering.

    Returns:
        Dict[str, np.ndarray]: A new dict where arrays aligned to `peaks` are filtered by `mask`.

    Raises:
        ValueError: If `mask` is not a 1D boolean array.
    """
    if mask.dtype != bool or mask.ndim != 1:
        raise ValueError("mask must be a 1D boolean array")
    out = {}
    for k, v in props.items():
        if isinstance(v, np.ndarray) and v.ndim == 1 and v.shape[0] == mask.shape[0]:
            out[k] = v[mask]
        else:
            out[k] = v
    return out


def peak_pick_fun(intensity: np.ndarray,
                  index: np.ndarray,
                  width: int = 2,
                  prominence: float | None = None,
                  snr: float  = 3,
                  noise: str = "wavelet",
                  relheight: float = 0.01,
                  method: str = 'scipy',
                  return_type: str = 'height'):

    _input_validation(intensity,index)

    if method == 'scipy':
        peaks, props = find_pick_scipy(intensity,
                                        index,
                                        width=width,
                                        prominence=prominence,
                                        snr=snr,
                                        noise=noise,
                                        relheight=relheight)
    else:
        logger.error("method must be 'scipy'")
        raise ValueError("method must be 'scipy'")

    if return_type == 'height':
        return intensity[peaks],index[peaks]

    elif return_type == 'area':
        return compute_peak_areas(intensity,index,peaks,props),index[peaks]
    else:
        logger.error("type must be 'height' or 'area'")
        raise ValueError("type must be 'height' or 'area'")


def find_pick_scipy(
    intensity: np.ndarray,
    index: np.ndarray,
    width: int | List[int] = 6,
    distance: int | None = 1,
    prominence: float | None = None,
    snr: float  = 2,
    noise: str = "wavlete",
    relheight: float = 0.3):
    """
    Peak picking using `scipy.signal.find_peaks` with SNR-based filtering.

    Parameters:
        intensity (np.ndarray): 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity`.
        width (int | List[int]): Required width(s) of peaks.
        distance (int | None): Required minimal horizontal distance between peaks.
        prominence (float | None): Required prominence of peaks; enables base indices in props.
        snr (float): Signal-to-noise threshold applied after peak detection.
        noise (str): Noise estimation method identifier used by `estimator`.
        relheight (float): Relative height ratio for `height` constraint.

    Returns:
        Tuple[np.ndarray, dict]: Indices of selected peaks and properties dict from `find_peaks`.

    Raises:
        ValueError: If inputs are invalid or method constraints fail downstream.
    """

    max_height = np.nanmax(intensity)
    relheight = relheight * max_height

    #width,relheight,prominence
    peaks, props = find_peaks(intensity,
                              width=width,
                              distance=distance,
                              height=[relheight,max_height],
                              prominence=prominence)

    #noise est
    noise_estimation = estimator(intensity,index,denoise_method=noise)

    snr_selection = (intensity[peaks] / noise_estimation) > snr
    peaks = peaks[snr_selection]
    props = _mask_peak_props(props, snr_selection)

    #snr selection
    return peaks, props

def _interp(arr: np.ndarray, idx_f: float, n: int | None = None) -> float:
    """
    Linear interpolation over an array for a fractional index.

    Parameters:
        arr (np.ndarray): Array to sample from.
        idx_f (float): Fractional index in [0, len(arr)-1].
        n (int | None): Optional length hint; if provided, must match arr length.

    Returns:
        float: Interpolated value at idx_f.

    Raises:
        ValueError: If n is provided and does not match arr length.
    """
    m = int(arr.shape[0])
    if n is not None and n != m:
        raise ValueError("Provided n does not match arr length")
    if m == 0:
        return float("nan")

    idx_f = float(np.clip(idx_f, 0, m - 1))
    i0 = int(np.floor(idx_f))
    i1 = min(i0 + 1, m - 1)
    alpha = idx_f - i0
    return float(arr[i0] + alpha * (arr[i1] - arr[i0]))

def compute_peak_areas(intensity: np.ndarray,
                       index: np.ndarray,
                       peaks: np.ndarray,
                       props: dict,
                       boundary: str = "ips",
                       ) -> np.ndarray:
    """
    Compute peak areas via piecewise linear interpolation and Simpson integration.

    Parameters:
        intensity (np.ndarray): 1D intensity array.
        index (np.ndarray): 1D coordinate array aligned with `intensity`.
        peaks (np.ndarray): Peak indices returned by `find_peaks`.
        props (dict): Peak properties from `find_peaks`; must contain
            `left_ips/right_ips` when `boundary='ips'` or
            `left_bases/right_bases` when `boundary='bases'`.
        boundary (str): One of {'ips','bases'} selecting integration boundaries.

    Returns:
        np.ndarray: Areas per peak, same length as `peaks`.

    Raises:
        ValueError: If required boundary keys are missing from `props` or boundary is invalid.
    """

    n = len(intensity)

    if boundary == "bases":
        if "left_bases" not in props or "right_bases" not in props:
            raise ValueError("Missing `left_bases/right_bases` in props; enable `prominence` during peak picking.")
        left_f = np.asarray(props["left_bases"], dtype=float)
        right_f = np.asarray(props["right_bases"], dtype=float)
    elif boundary == "ips":
        if "left_ips" not in props or "right_ips" not in props:
            raise ValueError("Missing `left_ips/right_ips` in props; compute widths in `find_peaks`.")
        left_f = np.asarray(props["left_ips"], dtype=float)
        right_f = np.asarray(props["right_ips"], dtype=float)
    else:
        raise ValueError("boundary must be one of {'ips','bases'}")

    areas = np.zeros(len(peaks), dtype=float)
    for k, (lf, rf) in enumerate(zip(left_f, right_f)):
        if rf <= lf:
            areas[k] = 0.0
            continue
        xl = _interp(index, lf,n)
        xr = _interp(index, rf,n)
        yl = _interp(intensity, lf,n)
        yr = _interp(intensity, rf,n)
        logger.debug(f"li={lf},ri={rf},peak={peaks[k]}")
        li = int(np.ceil(lf))
        ri = int(np.floor(rf))
        xs = [xl]
        ys = [yl]

        #add points between lf and rf
        if ri >= li:
            xs.extend(index[li:ri + 1].tolist())
            logger.debug(f"intensity[li:ri + 1]={intensity[li:ri + 1]}")
            ys.extend(intensity[li:ri + 1].tolist())

        xs.append(xr)
        ys.append(yr)
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)

        areas[k] = float(simpson(ys, xs))
        logger.debug(f"xs={xs}\r\nys={ys},area={areas[k]}")
    return areas
