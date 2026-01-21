import numpy as np
from numba import jit, prange, set_num_threads

set_num_threads(4)

@jit(nopython=True, fastmath=True, cache=True, parallel=True)
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
        for i in prange(p, n - p):
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


def snip_baseline_numba(
    intensity: np.ndarray,
    m: int | None = None,
    decreasing: bool = True,
) -> np.ndarray:
    """Compute the SNIP baseline for a 1D spectrum.

    Args:
        intensity: 1D intensity array.
        m: Maximum window radius; if None, choose based on spectrum length.
        decreasing: Whether to iterate window sizes from large to small.

    Returns:
        Baseline array (float32).
    """
    if intensity.ndim != 1:
        raise ValueError("intensity must be a 1D array")

    n = int(intensity.size)
    if n == 0:
        return np.array([], dtype=np.float32)
    if n < 3:
        return np.array(intensity, dtype=np.float32)

    if m is not None and m < 1:
        raise ValueError("m must be a positive integer for SNIP baseline")

    if m is None:
        m_in = min(100, max(10, n // 10))
    else:
        m_in = int(m)
    m_eval = max(1, min(m_in, (n - 1) // 2))

    xi = np.asarray(intensity, dtype=np.float64)
    xi = np.ascontiguousarray(xi)

    baseline = _snip_1d_core(xi, m_eval, decreasing)
    return baseline.astype(np.float32)