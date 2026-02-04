import numpy as np
from numba import jit
from massflow.tools.logger import get_logger

logger = get_logger("peak_pick_compute")

@jit(nopython=True, fastmath=True, cache=True)
def interp_jit(arr: np.ndarray, idx_f: float) -> float:
    """
    Linear interpolation over an array for a fractional index.

    Parameters:
        arr (np.ndarray): Array to sample from.
        idx_f (float): Fractional index in [0, len(arr)-1].
        n (int | None): Optional length hint; if provided, must match arr length.

    Returns:
        float: Interpolated value at idx_f.
    """
    m = int(arr.shape[0])
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
    Compute peak areas using Simpson's rule for numerical integration.

    Parameters:
        intensity (np.ndarray): Intensity array.
        index (np.ndarray): Index array (e.g., retention time).
        peaks (np.ndarray): Indices of peak apexes.
        left_f (np.ndarray): Indices of left bases of peaks.
        right_f (np.ndarray): Indices of right bases of peaks.

    Returns:
        np.ndarray: Array of computed peak areas.
    """
    areas = np.zeros(len(peaks), dtype=float)
    for k, (lf, rf) in enumerate(zip(left_f, right_f)):
        if rf <= lf:
            areas[k] = 0.0
            continue
        xl = interp_jit(index, lf)
        xr = interp_jit(index, rf)
        yl = interp_jit(intensity, lf)
        yr = interp_jit(intensity, rf)
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
    Numerical integration using Simpson's rule for discrete data.
    Compatible with Numba's nopython mode, supports non-uniformly spaced sampling points,
    and correctly handles boundary cases for odd and even number of points.

    Parameters:
        y: intensity, array of function values y[i] = f(x[i])
        x: index, array of independent variables, must be same length as y and monotonically increasing

    Returns:
        float: Approximate definite integral of the function over the x interval

    Compatibility:
        This function produces results consistent with scipy.integrate.simpson but runs faster
        (Numba JIT compiled, no Python interpreter overhead).
    """
    num = len(x)

    # Edge case: Cannot integrate with fewer than 2 points
    if num < 2:
        return 0.0

    # Edge case: Only 2 points, use the trapezoidal rule
    if num == 2:
        dx = x[1] - x[0]
        return 0.5 * dx * (y[0] + y[1])

    # Simpson's rule requires an even number of intervals (3 points form a group)
    if num % 2 == 1:
        stop = num - 2
    else:
        stop = num - 3

    result = 0.0

    for i in range(0, stop, 2):
        h0 = x[i+1] - x[i]
        h1 = x[i+2] - x[i+1]

        if h0 == 0 or h1 == 0:
            continue

        # Non-uniform Simpson integration formula
        # Area = (h0 + h1) / 6 * [ term0 * y[i] + term1 * y[i+1] + term2 * y[i+2] ]
        term0 = 2.0 - h1 / h0
        term1 = (h0 + h1)**2 / (h0 * h1)
        term2 = 2.0 - h0 / h1

        area = (h0 + h1) / 6.0 * (term0 * y[i] + term1 * y[i+1] + term2 * y[i+2])
        result += area

    # Even number of points: use parabolic fit for the last interval
    # Use the last three points to fit a parabola, then compute only the last interval's area
    if num % 2 == 0:
        slice1 = num - 1
        slice2 = num - 2
        slice3 = num - 3

        h_last = x[slice1] - x[slice2]
        h_prev = x[slice2] - x[slice3]

        if h_last > 0 and h_prev > 0:
            # Alpha: contribution weight of the last point
            num_a = 2 * h_last**2 + 3 * h_prev * h_last
            den_a = 6 * (h_last + h_prev)
            alpha = num_a / den_a

            # Beta: contribution weight of the middle point
            num_b = h_last**2 + 3 * h_prev * h_last
            den_b = 6 * h_prev
            beta = num_b / den_b

            # Eta: correction point weight
            num_e = h_last**3
            den_e = 6 * h_prev * (h_prev + h_last)
            eta = num_e / den_e

            result += alpha * y[slice1] + beta * y[slice2] - eta * y[slice3]

    return result
