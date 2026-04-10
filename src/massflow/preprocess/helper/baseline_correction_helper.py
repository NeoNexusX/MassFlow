from typing import Optional
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.interpolate import UnivariateSpline
from massflow.preprocess.numba.baseline_correction_numba import (
    baseline_locmin_numba,
    baseline_snip_numba,
    local_maxima_numba,
)
from massflow.preprocess.numba.noise_reduction_numba import smooth_lowess_numba
from massflow.tools import dispatch_with_supported_kwargs
from massflow.tools.logger import get_logger

logger = get_logger("preprocesss")


def _input_validation(intensity: np.ndarray, baseline_scale: float = 1.0):
    """
    Validate input parameters for baseline_correction functions.

    Parameters:
        intensity (np.ndarray): 1D intensity array to be preprocessed.
        baseline_scale (float): Scale factor applied to estimated baseline.

    Raises:
        TypeError: If intensity is not a numpy array.
        ValueError: If intensity is not 1D/empty, or baseline_scale is invalid.
    """
    # Validate intensity array
    if not isinstance(intensity, np.ndarray):
        logger.error("intensity must be a numpy array")
        raise TypeError("intensity must be a numpy array")

    elif intensity.ndim != 1 or intensity.size == 0:
        logger.error("intensity must be a non-empty 1D array")
        raise ValueError("intensity must be a non-empty 1D array")

    if not np.isfinite(baseline_scale) or not 0.0 < float(baseline_scale) <= 1.0:
        logger.error("baseline_scale must be a finite number in (0,1]")
        raise ValueError("baseline_scale must be a finite number in (0,1]")


def asls_baseline(
    intensity: np.ndarray,
    lam: float = 1e5,
    p: float = 0.001,
    niter: int = 15
) -> np.ndarray:
    """
    Asymmetric Least Squares (ASLS) baseline correction using sparse solver.

    Parameters:
        intensity (np.ndarray): 1D input signal.
        lam (float): Smoothness parameter; larger -> smoother baseline.
        p (float): Asymmetry parameter (0<p<1).
        niter (int): Iteration count.

    Returns:
        np.ndarray: Estimated baseline with the same dtype as input `intensity`
            (internal computations use float64).
    """
    # Parameter validation
    if not np.isfinite(lam) or lam <= 0:
        logger.error("lam must be a positive finite number for ASLS")
        raise ValueError("lam must be a positive finite number for ASLS")
    if not np.isfinite(p) or not 0.0 < p < 1.0:
        logger.error("p must be in (0,1) for ASLS")
        raise ValueError("p must be in (0,1) for ASLS")
    if not isinstance(niter, (int, np.integer)) or niter < 1:
        logger.error("niter must be a positive integer for ASLS")
        raise ValueError("niter must be a positive integer for ASLS")

    target_dtype = intensity.dtype
    n = intensity.size
    if n == 0:
        return np.array([], dtype=target_dtype)

    ones_n_minus_2 = np.ones(n - 2)
    d = diags(
        diagonals=[ones_n_minus_2, -2.0 * ones_n_minus_2, ones_n_minus_2],
        offsets=[0, 1, 2],  # type: ignore
        shape=(n - 2, n),
        dtype=np.float64,
    )
    # d^T d is a 5-diagonal sparse matrix
    dt_d = d.T @ d

    w = np.ones(n, dtype=np.float64)
    baseline = intensity.copy()

    for _ in range(int(max(1, niter))):
        w_mat = diags(w, 0, shape=(n, n), dtype=np.float64)
        z_mat = (w_mat + lam * dt_d).tocsc()
        rhs = w * intensity
        baseline = spsolve(z_mat, rhs)
        # asymmetric weights: points above baseline get smaller weight
        w = p * (intensity > baseline) + (1.0 - p) * (intensity <= baseline)

    return np.asarray(baseline, dtype=target_dtype)


def snip_baseline(
    intensity: np.ndarray,
    m: Optional[int] = None,
    decreasing: bool = True
) -> np.ndarray:
    """
    Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP) baseline estimation.

    This implementation follows the SNIP algorithm using vectorized NumPy operations:" to accurately reflect the actual implementation language
    - Uses a working buffer `y` and a temporary buffer `z`.
    - For each window size `p`, updates `z[i] = min(y[i], (y[i-p] + y[i+p]) / 2)`,
      then writes `z[i]` back to `y[i]` for indices in [p, n - p).
    - Supports decreasing (`p = m .. 1`) or increasing (`p = 1 .. m`) iteration order.

    Parameters:
        intensity (np.ndarray): Input 1D spectrum.
        m (Optional[int]): Window half-size. If None, auto-selects based on spectrum length.
        decreasing (bool): If True, iterate `p` in decreasing order.

    Returns:
        np.ndarray: Estimated baseline.

    Raises:
        ValueError: If `intensity` is not a 1D array.
    """
    # Validate input shape and length
    if intensity.ndim != 1:
        raise ValueError("intensity must be a 1D array")
    target_dtype = intensity.dtype
    n = int(intensity.size)
    if n == 0:
        return np.array([], dtype=target_dtype)
    if n < 3:
        # Too short to apply SNIP; return a copy as baseline
        return np.asarray(intensity, dtype=target_dtype)

    # Validate m when provided
    if m is not None and not isinstance(m, (int, np.integer)):
        logger.error("m must be an integer for SNIP baseline")
        raise TypeError("m must be an integer for SNIP baseline")
    if m is not None and m < 1:
        logger.error("m must be a positive integer for SNIP baseline")
        raise ValueError("m must be a positive integer for SNIP baseline")

    # Auto-select window: 10% of length, lower-bounded by 10 and capped at 100
    # Triple protection to keep indices valid and prevent out-of-bounds access
    m_in = int(m) if m is not None else min(100, max(10, n // 10))
    m_eval = max(1, min(m_in, (n - 1) // 2))

    # Working buffer baseline and temporary buffer z
    baseline = intensity.copy()
    z = np.empty_like(baseline)

    if decreasing:
        # Iterate p from m down to 1
        for p in range(m_eval, 0, -1):
            # Vectorized update for indices [p, n - p)
            cur = baseline[p : n - p]
            left = baseline[0 : n - 2 * p]
            right = baseline[2 * p : n]
            clip_val = 0.5 * (left + right)
            z[p : n - p] = np.minimum(cur, clip_val)
            baseline[p : n - p] = z[p : n - p]
    else:
        # Iterate p from 1 up to m
        for p in range(1, m_eval + 1):
            cur = baseline[p : n - p]
            left = baseline[0 : n - 2 * p]
            right = baseline[2 * p : n]
            clip_val = 0.5 * (left + right)
            z[p : n - p] = np.minimum(cur, clip_val)
            baseline[p : n - p] = z[p : n - p]

    return np.asarray(baseline, dtype=target_dtype)


def locmin_baseline(
    intensity: np.ndarray,
    smooth: str = "none",
    span: float = 0.1,
    s: Optional[float] = 0.0,
    upper: bool = False,
    width: int = 5,
) -> np.ndarray:
    """
    Baseline estimation by interpolation from local extrema.

    - Detect local minima (or maxima if `upper=True`) using a windowed rule (width).
    - Include both endpoints to ensure full-range coverage.
    - Linearly interpolate across anchor positions (R's approx behavior).
    - Optionally smooth the interpolated baseline using loess or spline.

    Parameters:
        intensity (np.ndarray): Input 1D spectrum.
        smooth (str): 'none' | 'loess' | 'spline'. Default 'none'.
        span (float): Loess span proportion (0 < span <= 1). Default 0.1.
        s (float): Spline smoothing target residual sum of squares; 0.0 means interpolation.
        upper (bool): If True, anchor at local maxima; otherwise local minima.
        width (int): Neighborhood width for extrema detection (default 5).

    Returns:
        np.ndarray: Estimated baseline with the same dtype as input `intensity`.
    """
    target_dtype = intensity.dtype
    n = intensity.size
    if n == 0:
        return np.array([], dtype=target_dtype)
    if n < 3:
        val = np.nanmax(intensity) if upper else np.nanmin(intensity)
        if not np.isfinite(val):
            val = 0.0
        return np.full((n,), float(val), dtype=target_dtype)

    # Validate smoothing parameters
    smooth_kind = (smooth or "none").strip().lower()
    if smooth_kind not in {"none", "loess", "spline"}:
        logger.error("smooth must be one of 'none', 'loess', or 'spline'")
        raise ValueError("smooth must be one of 'none', 'loess', or 'spline'")
    if smooth_kind == "loess":
        if not np.isfinite(span) or not (0.0 < span <= 1.0): # pylint: disable=superfluous-parens
            logger.error("span must be a finite number in (0,1] for loess smoothing")
            raise ValueError("span must be a finite number in (0,1] for loess smoothing")
    if smooth_kind == "spline": #有限且非负数
        if s is not None and (not np.isfinite(s) or s < 0.0):
            logger.error("s must be a finite number >= 0 for spline smoothing")
            raise ValueError("s must be a finite number >= 0 for spline smoothing")
    # Fill NaNs for robust processing
    med = np.nanmedian(intensity)
    if not np.isfinite(med):
        med = 0.0
    y_filled = np.where(np.isfinite(intensity), intensity, med)

    # Local extrema mask using numba local maxima: minima are maxima on -signal.
    width_eff = 3 if width < 3 else int(width)
    extrema_input = y_filled if upper else -y_filled
    mask, _ = local_maxima_numba(extrema_input, width=width_eff)
    locs_inner = np.where(mask)[0]

    # Always include endpoints
    locs = np.concatenate(([0], locs_inner, [n - 1]))
    locs = np.unique(locs)

    if locs.size >= 2:
        x = np.arange(n, dtype=np.float64)
        baseline = np.interp(x, locs, y_filled[locs])

        smooth_kind = (smooth or "none").strip().lower()
        if smooth_kind == "loess":
            span = float(span)
            _, baseline_smooth = smooth_lowess_numba(y=baseline, f=span)
            baseline = baseline_smooth
        elif smooth_kind == "spline":
            s_val = float(s) if s is not None else 0.0
            spl = UnivariateSpline(x, baseline, s=max(0.0, s_val))
            baseline = spl(x)

        return np.asarray(baseline, dtype=target_dtype)

    # Fallback when insufficient anchors: constant baseline
    val = np.nanmax(intensity) if upper else np.nanmin(intensity)
    if not np.isfinite(val):
        val = 0.0
    return np.full((n,), float(val), dtype=target_dtype)


def baseline_corrector(
    intensity: np.ndarray,
    method: str = "asls",
    smooth: str = "none",
    span: float = 0.1,
    s: Optional[float] = 0.0,
    upper: bool = False,
    width: int = 5,
    lam: float = 1e7,
    p: float = 0.01,
    niter: int = 15,
    baseline_scale: float = 1.0,
    m: Optional[int] = None,
    decreasing: bool = True,
    lengths: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Baseline estimation and correction for a single 1D intensity array.

    Parameters:
        intensity (np.ndarray): Input 1D intensity array maybe flat.
        method (str): Baseline method key. Supported keys are 'locmin', 'locmin_numba','snip', 'snip_numba', and 'asls'.
            Unknown keys fall back to 'asls'.
        smooth (str): LocMin smoothing method ('none', 'gaussian', etc.), used when method='locmin'.
        span (float): LocMin smoothing span (fraction of data), used when method='locmin'.
        s (float, optional): LocMin smoothing strength, used when method='locmin'.
        upper (bool): LocMin baseline type; if True, finds upper baseline, used when method='locmin'.
        width (int): LocMin local minima window width, used when method='locmin'.
        lam (float): ASLS smoothness parameter, used when method='asls'.
        p (float): ASLS asymmetry parameter, used when method='asls'.
        niter (int): ASLS iteration count, used when method='asls'.
        baseline_scale (float): Scale factor in (0,1] applied to estimated baseline.
        m (int): SNIP window half-size, used when method in {'snip','snip_numba'}.
        decreasing (bool): SNIP decreasing rule, used when method in {'snip','snip_numba'}.
        lengths (np.ndarray, optional): Segment lengths for flat input. Used when
            method='locmin_numba' to process concatenated spectra correctly.

    Notes:
        - Baseline function arguments are dispatched via `_dispatch_with_supported_kwargs`.
        - Only parameters present in the selected baseline function signature are forwarded.

    Returns:
        tuple[np.ndarray, np.ndarray]: (corrected, scaled_baseline), returned with
            the same dtype as input `intensity`; internal computations use float64.

    Raises:
        TypeError: If `intensity` is not a numpy array.
        ValueError: If `intensity` is invalid or `baseline_scale` is out of range.
    """
    _input_validation(intensity, baseline_scale=baseline_scale)
    target_dtype = intensity.dtype
    xi = np.asarray(intensity, dtype=np.float64)

    # Normalize method name (no explicit whitelist validation)
    method_norm = (method or "locmin").strip().lower()

    # Estimate baseline via method dispatch with kwargs filtered by target signature.
    method_map = {
        "locmin": locmin_baseline,
        "locmin_numba": baseline_locmin_numba,
        "snip": snip_baseline,
        "snip_numba": baseline_snip_numba,
        "asls": asls_baseline,
    }

    baseline_func = method_map.get(method_norm, asls_baseline)
    baseline = dispatch_with_supported_kwargs(
        baseline_func,
        intensity=xi,
        smooth=smooth,
        span=span,
        s=s,
        upper=upper,
        width=width,
        lam=lam,
        p=p,
        niter=niter,
        m=m,
        decreasing=decreasing,
        lengths=lengths,
    )

    # Scale baseline and subtract
    scaled_baseline = baseline_scale * baseline
    corrected = xi - scaled_baseline
    corrected = np.maximum(corrected, 0.0)

    corrected_out = corrected.astype(target_dtype, copy=False)
    scaled_baseline_out = scaled_baseline.astype(target_dtype, copy=False)

    return corrected_out, scaled_baseline_out
