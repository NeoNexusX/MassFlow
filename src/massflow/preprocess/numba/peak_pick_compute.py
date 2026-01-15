import numpy as np
from numba import jit
from massflow.tools.logger import get_logger

logger = get_logger("peak_pick_compute")

@jit(nopython=True, fastmath=True, cache=True)
def iterp_jit(arr: np.ndarray, idx_f: float) -> float:
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
    n = int(arr.shape[0])
    m = int(arr.shape[0])
    if n is not None and n != m:
        raise ValueError("Provided n does not match arr length")
    if m == 0:
        return float("nan")

    if idx_f < 0.0:
        idx_f = 0.0
    elif idx_f > m - 1:
        idx_f = float(m - 1)
    i0 = int(np.floor(idx_f))
    i1 = i0 + 1 if i0 + 1 < m else m - 1
    alpha = idx_f - i0

    return float(arr[i0] + alpha * (arr[i1] - arr[i0]))

@jit(nopython=True, fastmath=True, cache=True)
def compute_peak_areas_jit(intensity: np.ndarray,
                           index: np.ndarray,
                           peaks: np.ndarray,
                           left_f: np.ndarray,
                           right_f: np.ndarray,
                           ) -> np.ndarray:
    """
    Compute peak areas using the trapezoidal rule.

    Parameters:
        intensity (np.ndarray): Intensity array.
        index (np.ndarray): Index array (e.g., retention time).
        peaks (np.ndarray): Indices of peak apexes.
        left_bases (np.ndarray): Indices of left bases of peaks.
        right_bases (np.ndarray): Indices of right bases of peaks.

    Returns:
        np.ndarray: Array of computed peak areas.
    """
    areas = np.zeros(len(peaks), dtype=float)
    for k, (lf, rf) in enumerate(zip(left_f, right_f)):
        if rf <= lf:
            areas[k] = 0.0
            continue
        xl = iterp_jit(index, lf)
        xr = iterp_jit(index, rf)
        yl = iterp_jit(intensity, lf)
        yr = iterp_jit(intensity, rf)
        li = int(np.ceil(lf))
        ri = int(np.floor(rf))

        if ri >= li:
            mid_count = ri - li + 1
        else:
            mid_count = 0
        total_count = mid_count + 2

        xs = np.empty(total_count, dtype=float)
        ys = np.empty(total_count, dtype=float)

        xs[0] = xl
        ys[0] = yl

        if mid_count > 0:
            xs[1:1 + mid_count] = index[li:ri + 1]
            ys[1:1 + mid_count] = intensity[li:ri + 1]

        xs[-1] = xr
        ys[-1] = yr
        areas[k] = float(simpson_jit(ys, xs))
    return areas

@jit(nopython=True, cache=True)
def simpson_jit(y, x):
    """
    Numba-compatible implementation of scipy.integrate.simpson for 1D arrays.
    Guarantees the same output as scipy for 1D inputs with variable spacing.
    """
    num = len(x)

    # Edge case: Cannot integrate with fewer than 2 points
    if num < 2:
        return 0.0

    # Edge case: Only 2 points, use the trapezoidal rule
    if num == 2:
        dx = x[1] - x[0]
        return 0.5 * dx * (y[0] + y[1])

    # Determine the stopping point for the basic Simpson rule
    # If num is odd, we can cover all intervals (num-1 intervals, even segments)
    # If num is even, we can only use Simpson's rule up to num-3 (leaving the last segment for separate handling)
    if num % 2 == 1:
        stop = num - 2
    else:
        stop = num - 3

    result = 0.0

    # Basic Simpson Rule Loop (corresponds to SciPy's _basic_simpson)
    # Iterate over triplets of points (i, i+1, i+2) with a step of 2
    for i in range(0, stop, 2):
        h0 = x[i+1] - x[i]
        h1 = x[i+2] - x[i+1]

        # Avoid division by zero
        if h0 == 0 or h1 == 0:
            continue

        # Three-point non-uniform Simpson integration formula
        # Area = (h0 + h1) / 6 * [ term0 * y[i] + term1 * y[i+1] + term2 * y[i+2] ]
        term0 = 2.0 - h1 / h0
        term1 = (h0 + h1)**2 / (h0 * h1)
        term2 = 2.0 - h0 / h1

        area = (h0 + h1) / 6.0 * (term0 * y[i] + term1 * y[i+1] + term2 * y[i+2])
        result += area

    # Handle the last segment for even number of points
    if num % 2 == 0:
        # Use the last three points (num-3, num-2, num-1) to fit a parabola, but only calculate the area of the last segment (num-2 to num-1)
        # Corresponds to the alpha, beta, eta correction logic in SciPy source code

        slice1 = num - 1  # Last point
        slice2 = num - 2  # Second-to-last point
        slice3 = num - 3  # Third-to-last point

        h_last = x[slice1] - x[slice2] # h[1] in scipy
        h_prev = x[slice2] - x[slice3] # h[0] in scipy

        if h_last > 0 and h_prev > 0:
            # Alpha (corresponds to y[slice1])
            num_a = 2 * h_last**2 + 3 * h_prev * h_last
            den_a = 6 * (h_last + h_prev)
            alpha = num_a / den_a

            # Beta (corresponds to y[slice2])
            num_b = h_last**2 + 3 * h_prev * h_last
            den_b = 6 * h_prev
            beta = num_b / den_b

            # Eta (corresponds to y[slice3])
            num_e = h_last**3
            den_e = 6 * h_prev * (h_prev + h_last)
            eta = num_e / den_e

            result += alpha * y[slice1] + beta * y[slice2] - eta * y[slice3]

    return result
