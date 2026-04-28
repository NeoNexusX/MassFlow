import numpy as np
from numba import njit
from numpy.typing import NDArray


@njit(cache=True)
def reference_count_filter_mask_jit(
    counts: NDArray[np.int64],
    n_spectra: int,
    use_min_count: bool,
    min_count: int,
    use_min_frequency: bool,
    min_frequency: float,
) -> NDArray[np.bool_]:
    """Build a keep mask from reference support counts."""
    keep = np.ones(counts.size, dtype=np.bool_)

    for i in range(counts.size):
        count = counts[i]
        if use_min_count and count < min_count:
            keep[i] = False
            continue

        if use_min_frequency:
            if n_spectra <= 0:
                keep[i] = False
                continue
            if count / float(n_spectra) < min_frequency:
                keep[i] = False

    return keep
