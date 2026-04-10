import numpy as np
from numba import jit, njit, prange
from massflow.preprocess.numba.noise_reduction_numba import _lowess_core
from massflow.tools import lengths_to_offsets, prepare_flat_inputs

@njit(cache=True, fastmath=True)
def _median_of_finite(values: np.ndarray) -> float:
    """Return median of finite entries; fallback to 0.0 when no finite value exists."""
    n = values.size
    finite_count = 0
    for i in range(n):
        if np.isfinite(values[i]):
            finite_count += 1

    if finite_count == 0:
        return 0.0

    finite_vals = np.empty(finite_count, dtype=np.float64)
    k = 0
    for i in range(n):
        v = values[i]
        if np.isfinite(v):
            finite_vals[k] = v
            k += 1

    med = np.median(finite_vals)
    if not np.isfinite(med):
        return 0.0
    return med


@njit(cache=True, fastmath=True)
def _constant_extrema_baseline(y: np.ndarray, upper: bool) -> np.ndarray:
    """Build constant baseline using finite extrema with robust fallback."""
    n = y.size
    out = np.empty(n, dtype=np.float64)
    if n == 0:
        return out

    has_finite = False
    val = 0.0
    for i in range(n):
        yi = y[i]
        if np.isfinite(yi):
            if not has_finite:
                val = yi
                has_finite = True
            elif upper and yi > val:
                val = yi
            elif (not upper) and yi < val:
                val = yi

    if not has_finite or (not np.isfinite(val)):
        val = 0.0

    for i in range(n):
        out[i] = val
    return out


@njit(cache=True, fastmath=True)
def baseline_locmin_core(
    intensity: np.ndarray,
    width: int,
    upper: bool,
    apply_lowess: int,
    span: float,
    iter_count: int,
) -> np.ndarray:
    """Single-spectrum locmin core. Computes one pixel spectrum baseline in float64."""
    n = intensity.size
    if n == 0:
        return np.empty(0, dtype=np.float64)
    if n < 3:
        return _constant_extrema_baseline(intensity, upper)

    med = _median_of_finite(intensity)
    y_filled = np.empty(n, dtype=np.float64)
    for i in range(n):
        yi = intensity[i]
        y_filled[i] = yi if np.isfinite(yi) else med

    width_eff = 3 if width < 3 else int(width)
    extrema_input = np.empty(n, dtype=np.float64)
    if upper:
        for i in range(n):
            extrema_input[i] = y_filled[i]
    else:
        for i in range(n):
            extrema_input[i] = -y_filled[i]

    mask_u8 = np.zeros(n, dtype=np.uint8)
    _local_maxima_core_numba(extrema_input, width_eff, mask_u8)

    baseline = np.empty(n, dtype=np.float64)
    left = 0
    for right in range(1, n):
        if mask_u8[right] != 1 and right != (n - 1):
            continue

        y_left = y_filled[left]
        y_right = y_filled[right]
        slope = (y_right - y_left) / float(right - left)
        for i in range(left, right + 1):
            baseline[i] = y_left + slope * float(i - left)
        left = right

    if apply_lowess == 1:
        x = np.arange(n, dtype=np.float64)
        delta = 0.01 * (x[n - 1] - x[0]) if n > 1 else 0.0
        baseline = _lowess_core(x, baseline, span, iter_count, delta)

    return baseline


@njit(parallel=True, cache=True, fastmath=True)
def _baseline_locmin_jit(
    flat: np.ndarray,
    lengths: np.ndarray,
    width: int,
    upper: bool,
    apply_lowess: int,
    span: float,
    iter_count: int,
) -> np.ndarray:
    """Flat batch entry. Dispatch each pixel spectrum to baseline_locmin_core in parallel."""
    out = np.empty(flat.size, dtype=np.float64)
    offsets = lengths_to_offsets(lengths)

    for idx in prange(lengths.size):  # pylint: disable=not-an-iterable
        start = offsets[idx]
        end = offsets[idx + 1]
        if end <= start:
            continue
        out[start:end] = baseline_locmin_core(
            flat[start:end],
            width,
            upper,
            apply_lowess,
            span,
            iter_count,
        )

    return out

@njit(cache=True, fastmath=True)
def _snip_1d_core(intensity: np.ndarray, m_eval: int, decreasing: bool) -> np.ndarray:
    """SNIP core routine that iteratively clips baseline estimates by window p."""
    n = intensity.size
    baseline = intensity.copy()
    z = np.empty_like(baseline)

    if decreasing:
        p_start = m_eval
        p_end = 0
        p_step = -1
    else:
        p_start = 1
        p_end = m_eval + 1
        p_step = 1

    for p in range(p_start, p_end, p_step):
        for i in range(p, n - p):
            left = baseline[i - p]
            right = baseline[i + p]
            clip_val = 0.5 * (left + right)
            cur = baseline[i]
            if cur <= clip_val:
                z[i] = cur
            else:
                z[i] = clip_val

        for i in range(p):
            z[i] = baseline[i]
        for i in range(n - p, n):
            z[i] = baseline[i]

        for i in range(n):
            baseline[i] = z[i]

    return baseline


@njit(parallel=True, cache=True, fastmath=True)
def baseline_snip_jit(
    flat: np.ndarray,
    lengths: np.ndarray,
    m_eval: int,
    decreasing: bool,
) -> np.ndarray:
    """Parallel SNIP dispatcher for flat spectra."""
    out = np.empty(flat.size, dtype=np.float64)
    offsets = lengths_to_offsets(lengths)

    for idx in prange(lengths.size):  # pylint: disable=not-an-iterable
        start = offsets[idx]
        end = offsets[idx + 1]
        seg_len = end - start

        if seg_len <= 0:
            continue
        if seg_len < 3:
            out[start:end] = flat[start:end]
            continue

        m_seg = m_eval
        m_cap = (seg_len - 1) // 2
        if m_seg > m_cap:
            m_seg = m_cap
        if m_seg < 1:
            m_seg = 1

        out[start:end] = _snip_1d_core(flat[start:end], m_seg, decreasing)

    return out


@jit(nopython=True, cache=True)
def _local_maxima_core_numba(
    intensity: np.ndarray, width: int, buffer: np.ndarray
) -> int:
    """Numba core routine mirroring the C local_maxima loop behavior."""
    n = int(intensity.size)
    nmax = 0
    r = abs(int(width) // 2)

    for i in range(n):
        # Keep endpoints excluded, matching the Python helper behavior.
        if i < r or i >= n - r:
            continue

        a = i - r
        b = i + r
        is_max = False

        for j in range(a, b + 1):
            if intensity[i] > intensity[j]:
                is_max = True
            if j < i and intensity[j] >= intensity[i]:
                is_max = False
                break
            if j > i and intensity[j] > intensity[i]:
                is_max = False
                break

        if is_max:
            buffer[i] = 1
            nmax += 1

    return nmax


# ==================== Non-JIT wrappers and preprocessing ====================


def local_maxima_numba(
    intensity: np.ndarray, width: int = 5
) -> tuple[np.ndarray, int]:
    """Find local maxima mask and count using a numba-accelerated C-style rule.

    Args:
        intensity: 1D intensity array.
        width: Neighborhood width; radius is abs(width // 2).

    Returns:
        A tuple (mask, nmax), where mask is bool array with local maxima marks.
    """
    if intensity.ndim != 1 or not isinstance(width, (int, np.integer)) or width < 3:
        raise ValueError("width must be >= 3 for local maxima detection,width must be an integer, and intensity must be 1D")

    xi = np.ascontiguousarray(intensity)

    buffer = np.zeros(xi.size, dtype=np.uint8)
    nmax = _local_maxima_core_numba(xi, int(width), buffer)
    return buffer.astype(np.bool_), int(nmax)


def baseline_snip_numba(
    intensity: np.ndarray,
    lengths: np.ndarray | None = None,
    m: int | None = None,
    decreasing: bool = True,
) -> np.ndarray:
    """Compute the SNIP baseline for a 1D spectrum.

    Args:
        intensity: 1D intensity array.
        m: Maximum window radius; if None, choose based on spectrum length.
        decreasing: Whether to iterate window sizes from large to small.
        numba_max_threads: Maximum number of threads for Numba parallel execution.

    Returns:
        Baseline array (float32).
    """

    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        raise ValueError("intensity must be a 1D numpy array")

    target_dtype = intensity.dtype
    intensity_arr, lengths_arr = prepare_flat_inputs(intensity, lengths)

    if lengths_arr.ndim != 1:
        raise ValueError("lengths must be a 1D array")
    if lengths_arr.size == 0:
        return np.array([], dtype=target_dtype)

    lengths_i64 = np.ascontiguousarray(lengths_arr.astype(np.int64, copy=False))
    if np.any(lengths_i64 < 0):
        raise ValueError("lengths must contain non-negative integers")

    total_len = int(lengths_i64.sum())
    if total_len != int(intensity_arr.size):
        raise ValueError("sum(lengths) must equal intensity.size")

    n = int(intensity_arr.size)
    if n == 0:
        return np.array([], dtype=target_dtype)

    if m is not None and m < 1:
        raise ValueError("m must be a positive integer for SNIP baseline")

    if m is None:
        m_in = min(100, max(10, n // 10))
    else:
        m_in = int(m)
    m_eval = max(1, min(m_in, (n - 1) // 2))

    xi = np.asarray(intensity_arr, dtype=np.float64)
    xi = np.ascontiguousarray(xi)

    baseline = baseline_snip_jit(xi, lengths_i64, m_eval, decreasing)
    return baseline.astype(target_dtype, copy=False)


def baseline_locmin_numba(
    intensity: np.ndarray,
    lengths: np.ndarray | None = None,
    width: int = 5,
    upper: bool = False,
    span: float = 0.1,
    niter: int = 3,
    smooth: str | None = None,
) -> np.ndarray:
    """
    Flat locmin baseline wrapper with basic validation and dtype-preserving output.
    """

    if not isinstance(intensity, np.ndarray) or intensity.ndim != 1:
        raise ValueError("intensity must be a 1D numpy array")

    width_i = int(width)
    apply_lowess = 0
    if width_i < 3:
        raise ValueError("width must be >= 3")
    if smooth is not None:
        smooth_kind = (smooth or "none").strip().lower()
        if smooth_kind not in {"none", "loess"}:
            raise ValueError("locmin_numba only supports smooth='none' or 'loess'")
        apply_lowess = 1 if smooth_kind == "loess" else 0

    if apply_lowess == 1:
        if (not np.isfinite(span)) or (span <= 0.0) or (span > 1.0):
            raise ValueError("span must be finite and in (0, 1] when apply_lowess=1")
        if (not isinstance(niter, (int, np.integer))) or int(niter) < 0:
            raise ValueError("niter must be an integer >= 0")

    iter_i = int(niter)
    target_dtype = intensity.dtype
    intensity_arr, lengths_arr = prepare_flat_inputs(intensity, lengths)

    if lengths_arr.size == 0:
        return np.array([], dtype=target_dtype)

    lengths_i64 = np.ascontiguousarray(lengths_arr.astype(np.int64, copy=False))

    total_len = int(lengths_i64.sum())
    if total_len != int(intensity_arr.size):
        raise ValueError("sum(lengths) must equal intensity.size")

    xi = np.ascontiguousarray(intensity_arr.astype(np.float64, copy=False))
    baseline64 = _baseline_locmin_jit(
        xi,
        lengths_i64,
        width_i,
        bool(upper),
        apply_lowess,
        float(span),
        iter_i,
    )

    return baseline64.astype(target_dtype, copy=False)
