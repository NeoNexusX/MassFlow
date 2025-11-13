"""
Author: MassFlow Development Team Bionet/NeoNexus
License: See LICENSE file in project root

Noise estimation utilities.
"""
from typing import Optional
import numpy as np
from scipy import stats
from scipy.interpolate import InterpolatedUnivariateSpline
from preprocess.filter_helper import smoother
from logger import get_logger

logger = get_logger(__name__)

def _input_validation(intensity,index,nbins,overlap):

    if intensity is None or intensity.ndim != 1:
        logger.error("x.intensity must be a 1D numpy array")
        raise TypeError("x.intensity must be a 1D numpy array")
    if index is not None and (index.ndim !=1 or len(index)!=len(intensity)):
        logger.error("index must be a 1D numpy array of the same length as intensity")
        raise TypeError("index must be a 1D numpy array of the same length as intensity")
    if nbins < 1:
        logger.error("nbins must be >= 1")
        raise ValueError("nbins must be >= 1")
    if overlap > 1 or overlap <0:
        logger.error("overlap must be in [0,1]")
        raise ValueError("overlap must be in [0,1]")


def _linear_sse(data: np.ndarray, i: int, j: int) -> float:
    """Compute SSE of a first-degree polynomial fit in the interval [i, j].

    Args:
        data (np.ndarray): 1D signal array to fit.
        i (int): Lower index (inclusive).
        j (int): Upper index (inclusive).

    Returns:
        float: Sum of squared residuals for the linear fit. NaNs are masked.

    Raises:
        ValueError: If ``i > j``.
    """
    if i > j:
        raise ValueError("Invalid interval: i > j")
    n = len(data)
    i = max(0, i)
    j = min(n - 1, j)
    y = data[i:j + 1]
    if y.size == 0:
        return 0.0
    mask = ~np.isnan(y)
    if not np.any(mask):
        return 0.0
    y = y[mask]
    x = np.arange(i, j + 1)[mask].astype(float)
    # Closed-form linear regression (least squares)
    x_mean = np.nanmean(x)
    y_mean = np.nanmean(y)
    x_centered = x - x_mean
    y_centered = y - y_mean
    denom = np.dot(x_centered, x_centered)
    if denom == 0.0:
        # Degenerate case (single point): best fit is mean
        residual = y - y_mean
        return float(np.dot(residual, residual))
    slope = np.dot(x_centered, y_centered) / denom
    intercept = y_mean - slope * x_mean
    y_hat = slope * x + intercept
    residual = y - y_hat
    return float(np.dot(residual, residual))


def _bins_sse(data: np.ndarray, lowers: np.ndarray, uppers: np.ndarray) -> np.ndarray:
    """Compute SSE for each bin defined by ``lowers``/``uppers``.

    Args:
        data (np.ndarray): 1D signal array to fit.
        lowers (np.ndarray): Lower indices for bins (inclusive).
        uppers (np.ndarray): Upper indices for bins (inclusive).

    Returns:
        np.ndarray: SSE values per bin.
    """
    return np.array([
        _linear_sse(data, int(li), int(ui))
        for li, ui in zip(lowers, uppers)
    ], dtype=float)


def _merge_best_pair(lowers: np.ndarray, uppers: np.ndarray, sse: np.ndarray):
    """Merge the adjacent pair whose combined SSE is minimal.

    Args:
        lowers (np.ndarray): Current bin lower bounds.
        uppers (np.ndarray): Current bin upper bounds.
        sse (np.ndarray): Current per-bin SSEs.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Updated ``lowers``, ``uppers`` after one merge.
    """
    if lowers.size <= 1:
        return lowers, uppers
    pair_scores = sse[:-1] + sse[1:]
    p = int(np.argmin(pair_scores))
    new_lowers = np.delete(lowers, p + 1)
    new_uppers = np.delete(uppers, p)
    # Merge bins p and p+1
    new_lowers[p] = lowers[p]
    new_uppers[p] = uppers[p + 1]
    return new_lowers, new_uppers


def _best_split_index(data: np.ndarray, i: int, j: int) -> int:
    """Find a split point in [i, j] that minimizes SSE(left) + SSE(right).

    Args:
        data (np.ndarray): 1D signal array to fit.
        i (int): Lower index (inclusive).
        j (int): Upper index (inclusive).

    Returns:
        int: Split index ``t`` such that the left bin is ``[i, t-1]`` and right is ``[t, j]``.
             Returns ``-1`` if no valid split is found.
    """
    if j - i + 1 < 2:
        return -1
    best_t = -1
    best_score = np.inf
    for t in range(i + 1, j + 1):
        s_left = _linear_sse(data, i, t - 1)
        s_right = _linear_sse(data, t, j)
        score = s_left + s_right
        if score < best_score:
            best_score = score
            best_t = t
    return best_t


def _split_worst_bin(data: np.ndarray, lowers: np.ndarray, uppers: np.ndarray, sse: np.ndarray):
    """Split the bin with the highest SSE using the best split point.

    Args:
        data (np.ndarray): 1D signal array to fit.
        lowers (np.ndarray): Current bin lower bounds.
        uppers (np.ndarray): Current bin upper bounds.
        sse (np.ndarray): Current per-bin SSEs.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Updated lowers, uppers after one split.

    Raises:
        ValueError: If no valid split point can be found for any bin.
    """
    order = np.argsort(-sse)  # descending by SSE
    for k in order:
        i, j = int(lowers[k]), int(uppers[k])
        t = _best_split_index(data, i, j)
        if t != -1:
            new_lowers = np.insert(lowers, k + 1, t)
            new_uppers = np.insert(uppers, k, t - 1)
            # Adjust original bin upper/lower to match the split
            new_lowers[k] = i
            new_uppers[k] = t - 1
            new_lowers[k + 1] = t
            new_uppers[k + 1] = j
            return new_lowers, new_uppers
    raise ValueError("Failed to find a valid split for the worst bin.")


def _findbins(data: np.ndarray,
              nbins: int = 1,
              dynamic: bool = False,
              niter: int = 10,
              overlap: float = 0.5,
              limits_only: bool = False):
    if data.ndim != 1:
        logger.error("Input data must be a 1D array")
        raise ValueError(f"Input data must be a 1D array, got shape {data.shape}")
    if nbins < 1:
        logger.warning(f"nbins={nbins} is not >= 1, clipped to 1")
        nbins = 1
    if overlap < 0.0 or overlap >= 0.99:
        logger.warning(f"overlap={overlap} is not in [0.0, 0.99), clipped to 0.5")
        overlap = 0.5

    # Unified output extras for dynamic mode; remain None in static mode
    sse = None
    trace = None

    n = len(data)
    nbins = min(nbins, n)  # clip to data length

    if overlap > 0 and not dynamic:
        # Bin size such that total covered length equals n after accounting for overlap 
        width = n / (nbins * (1 - overlap) + overlap)
        # Successive bin start step after overlap
        step = width * (1 - overlap)
        starts = np.arange(nbins, dtype=float) * step
        # 0-based inclusive boundaries; use floor for lower and ceil for upper for robust coverage
        lower = np.clip(np.floor(starts).astype(int), 0, n - 1)
        upper = np.clip(np.ceil(starts + width).astype(int) - 1, 0, n - 1)
        # Warn if any bins collapsed due to extreme parameters
        if np.any(upper < lower):
            logger.warning("Collapsed bins detected under high overlap/nbins; consider adjusting parameters.")

    elif overlap == 0 :
        edges = np.floor(np.linspace(0, n, num=nbins + 1, dtype=float)).astype(int)
        # 0-based inclusive boundaries; use floor for lower and ceil for upper for robust coverage
        lower = np.clip(edges[:-1], 0, n - 1)
        upper = np.clip(edges[1:] - 1, 0, n - 1)

    elif dynamic:
        # Dynamic binning ignores overlap and starts from equal-width bins
        if nbins < 3:
            logger.warning(f"nbins={nbins} is not >= 3, clipped to 3")
            nbins = 3

        # Equal-width initial bins
        edges = np.floor(np.linspace(0, n, num=nbins + 1, dtype=float)).astype(int)
        lower = np.clip(edges[:-1], 0, n - 1)
        upper = np.clip(edges[1:] - 1, 0, n - 1)

        # Initialize SSE and trace
        sse = _bins_sse(data, lower, upper)
        trace = [float(np.sum(sse))]

        # Optimization loop: merge-best then split-worst
        for _ in range(int(max(1, niter))):
            merged_lower, merged_upper = _merge_best_pair(lower, upper, sse)
            merged_sse = _bins_sse(data, merged_lower, merged_upper)
            try:
                split_lower, split_upper = _split_worst_bin(data, merged_lower, merged_upper, merged_sse)
            except ValueError:
                # If split fails, stop early
                break
            new_sse = _bins_sse(data, split_lower, split_upper)
            new_score = float(np.sum(new_sse))

            if new_score < trace[-1]:
                # Accept update
                lower, upper = split_lower, split_upper
                sse = new_sse
                trace.append(new_score)
            else:
                trace.append(trace[-1])
                break

        # Do not return here; fall through to unified return block below

    # Unified return block: build meta and optionally include dynamic extras
    meta = {
        "lower": lower.astype(int),
        "upper": upper.astype(int),
        "size": (upper - lower + 1).astype(int),
    }
    if sse is not None:
        meta["sse"] = sse.astype(float)
    if trace is not None:
        meta["trace"] = np.array(trace, dtype=float)
    if limits_only:
        return meta
    else:
        bins = [data[int(li):int(ui) + 1] for li, ui in zip(lower, upper)]
        return bins, meta


def estimation_fun(data,method: str = 'sd'):
    """
    Return a callable for noise estimation from a vector.
    """

    if method == 'sd':
        return np.nanstd(data)
    elif method == 'mad':
        return stats.median_abs_deviation(data)
    elif method == 'quantile':
        return np.nanquantile(data,0.95)
    elif method == 'diff':
        return np.nanmean(np.abs(data - np.nanmean(data)))
    else:
        raise ValueError(f"Unknown noise estimation method: {method}")


def estimator(intensity: np.ndarray,
              indexes: np.ndarray,
              nbins: int = 1,
              overlap: float = 0.5,
              dynamic: bool = False,
              method: str = 'sd',
              denoise_method: str = 'bi_ns'
    ):
    """
    Estimate noise level in the MSI data.
    """
    #basic input validation
    _input_validation(intensity,indexes,nbins,overlap)

    # Smooth signal (neighborhood-search Gaussian) and compute absolute residuals
    smoothed = smoother(intensity,indexes,method=denoise_method)

    residuals = np.abs(smoothed - intensity)

    #find bins
    if nbins > 1:
        bins, meta = _findbins(residuals, nbins=nbins, overlap=overlap,dynamic=dynamic)
        lower = meta["lower"].astype(int)
        upper = meta["upper"].astype(int)

        #core location
        midpoints = (lower + upper) // 2

        #noise_estimation for all bins
        noise_estimations = []
        for bin_data in bins:
            noise_estimations.append(estimation_fun(bin_data,method=method))

        #spline for noise line
        rank_spline = int(max(1, min(3, len(midpoints) - 1))) #此处要验证
        spline_fn = InterpolatedUnivariateSpline(midpoints, noise_estimations, k=rank_spline)
        noise_estimation = spline_fn(np.arange(len(residuals), dtype=float))
    else:
        #normal noise_estimation part
        noise_estimation = estimation_fun(residuals,method=method)
        noise_estimation = max(noise_estimation,0.1)

    return noise_estimation


def calculate_snr(intensity: np.ndarray,
                  indexes: Optional[np.ndarray] = None,
                  method="mad",
                  denoise_method: str = "bi_ns",):
    """
    pass
    """
    signal_level = np.percentile(intensity, 95)

    noise = estimator(intensity,
                      indexes,
                      method=method,
                      denoise_method=denoise_method)

    logger.info(f"SNR: signal_level:{signal_level}, noise:{noise}")
    return signal_level / noise
