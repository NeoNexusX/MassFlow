"""
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root
"""
from scipy.signal import find_peaks
import numpy as np
from logger import get_logger
from module.ms_module import SpectrumBaseModule

logger = get_logger(__name__)

def peak_pick_fun(method: str = 'scipy'):
    if method == 'scipy':
        return peak_pick_scipy



def peak_pick_scipy(
    x: np.ndarray,
    width: int = 5,
    prominence: float | None = None,
    snr: float  = 2,
    noise: str = "quant",
    relheight: float = 0.005,
    wlen: int | None = None):
    """Peak pick using scipy.signal.find_peaks."""

    if x.ndim != 1:
        logger.error("Input data must be a 1D array")
        raise ValueError("Input data must be a 1D array")

    if prominence is None:
        prominence = 0.005 * np.nanmax(x)

    peaks, props = find_peaks(x,
                              distance=width,
                              prominence=prominence,
                              relheight=relheight,
                              )
    return peaks, props


