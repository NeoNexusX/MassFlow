import numpy as np
from module.ms_module import MS, SpectrumBaseModule
from logger import get_logger
from typing import Optional, Tuple, Literal, List
from scipy.signal import find_peaks, lfilter

logger = get_logger("peak alignment")

def get_reference_mz_axis(
        ms_data: MS,
        agg,
        round_digits,
        half_window,
        findlimits,
        snr_threshold,
        noise_method,
        min_height,
        n_sample: int,
        min_distance_da: Optional[float] = None,
        min_distance_factor: float = 0.5,
        ) -> np.ndarray:
    """
    Detect peaks on the global mean spectrum to build the reference m/z axis.

    This function is the core of peak alignment, extracting a unified m/z axis from all pixels.

    Steps:
    a) compute the global mean spectrum (`_compute_mean_spectrum`)
    b) run local maxima detection on the mean spectrum (`_find_peaks_from_array`)
    c) obtain reference peaks and deduplicate/sort (merge peaks that are too close)

    Args:
        ms_data (MS):
            Input MS object containing spectral data for all pixels (usually sparse centroided peaks).
        agg (str):
            Aggregation method for the global mean spectrum: 'mean' or 'sum'.
        round_digits (int):
            Number of decimal places to round m/z when forming the mean spectrum. Helps merge m/z jitter.
        half_window (int):
            Half window size (in indices) used by `_find_peaks_from_array`. Smaller values (e.g. 1) are more sensitive to shoulders.
        findlimits (bool):
            Whether `_find_peaks_from_array` should also return left/right limits for plateaus. Typically False here, we only need peak centers.
        snr_threshold (Optional[float]):
            Signal-to-noise ratio threshold. If provided, `_find_peaks_from_array` discards peaks with SNR below this value.
        noise_method (str):
            Noise estimation method when `snr_threshold` is enabled: 'mad' (robust, default) or 'std' (sensitive to strong peaks).
        min_height (Optional[float]):
            Minimum absolute intensity threshold; peaks below this are discarded by `_find_peaks_from_array`.
        n_sample (int):
            Number of spectra sampled to estimate the deduplication distance (`min_distance_da`) when None.
        min_distance_da (Optional[float]):
            Minimum m/z distance (Da) for deduplication/merging. If None, `_estimate_mass_resolution` with `min_distance_factor` is used to auto-estimate.
        min_distance_factor (float):
            Scaling factor applied when `min_distance_da` is None. Smaller values deduplicate more conservatively, retaining more peaks (default 0.5).

    Returns:
        np.ndarray:
            A 1D strictly increasing NumPy array representing the final reference m/z axis.
    """
    # Compute the global mean spectrum
    bin_centers, mean_intensity = _compute_mean_spectrum(ms_data,agg,round_digits)
    if bin_centers.size == 0:
        logger.warning("Unable to generate reference axis: average mass spectrum is empty.")
        return np.asarray([], dtype=np.float64)

    # Local maxima detection (small window controlled by `half_window`)
    peaks_idx = _find_peaks_from_array(mean_intensity, half_window, findlimits, snr_threshold, noise_method, min_height)
    if peaks_idx.size == 0:
        logger.info("No local maxima detected on the average mass spectrum, returning empty reference axis.")
        return np.asarray([], dtype=np.float64)

    ref = bin_centers[peaks_idx]

    # Deduplicate and sort
    ref = np.sort(np.asarray(ref, dtype=np.float64))
    if ref.size > 1:
        # If `min_distance_da` is not provided, use the estimated Da tolerance (half width)
        if min_distance_da is None:
            da_tol = _estimate_mass_resolution(ms_data=ms_data, units='da', n_sample=n_sample)
            if not np.isfinite(da_tol) or da_tol <= 0:
                da_tol = 0.02  # Default fallback consistent with tolerance estimation
            min_distance_da_value = float(da_tol * float(min_distance_factor))
            logger.info(f"Reference axis deduplication threshold linked: min_distance_da={min_distance_da_value:.6f} Da (factor={min_distance_factor})")
        else:
            min_distance_da_value = float(min_distance_da)

        keep = [0]
        for i in range(1, ref.size):
            if (ref[i] - ref[keep[-1]]) > min_distance_da_value:
                keep.append(i)
        ref = ref[keep]
    return ref

def get_reference_mz_axis2(
    ms_data: MS,
    mz_res: Optional[float] = None,
    px_perc: float = 0.01,
    N_sample: Optional[int] = None,
    smoothing_window: int = 10,
) -> np.ndarray:
    """
    Generate reference m/z axis based on "global m/z histogram + LOESS smoothing + peak detection" (Histogram/Loess method).

    Processing idea:
    - Collect all centroid m/z points from a random sample of pixels, building a global m/z distribution histogram (bin width ~ `mz_res`).
    - Apply LOESS smoothing to the histogram to reduce jitter and discreteness.
    - Detect peaks on the smoothed curve, then filter by original histogram height (>= `px_perc`).

    Args:
        ms_data (MS):
            Input MS object containing spectral data for all pixels (usually sparse centroided peaks).
        mz_res (float):
            Expected m/z resolution (controls histogram bin width).
        px_perc (float):
            Minimum percentage of pixels that must contain a peak (approximate threshold), e.g. 0.01 for ≥1%.
        N_sample (Optional[int]):
            Number of pixels to sample for building the histogram. None for full dataset (consistent with consensus method).
        smoothing_window (int):
            Smoothing window length (odd number best, default 10).

    Returns:
        np.ndarray:
            Reference m/z axis (strictly increasing, float64).
    """
    # Basic checks
    N_pix = len(ms_data)
    if N_pix == 0:
        logger.warning("Empty dataset, returning empty reference axis.")
        return np.asarray([], dtype=np.float64)

    # Sampling strategy: Default to full sampling to align with script behavior.
    if N_sample is None:
        N_sample = N_pix
    else:
        N_sample = int(max(1, min(N_sample, N_pix)))

    # Random sampling without replacement (N_sample == N_pix is equivalent to full sampling)
    rng = np.random.default_rng(0)
    sample_idx = np.arange(N_pix, dtype=int) if N_sample == N_pix else rng.choice(N_pix, size=N_sample, replace=False)

    # Collect all m/z points from the sampled pixels (filter invalid and non-positive intensities)
    all_mz: List[float] = []
    for pid in sample_idx:
        spec = ms_data[int(pid)]
        mz_raw = spec.mz_list
        if mz_raw is None:
            continue
        mz = np.asarray(mz_raw, dtype=np.float64)
        if mz.size == 0:
            continue
        mask = np.isfinite(mz)
        if not mask.all():
            mz = mz[mask]
        if mz.size == 0:
            continue
        all_mz.extend(mz.tolist())

    if not all_mz:
        logger.warning("No valid m/z points collected from the sampled pixels, returning empty reference axis.")
        return np.asarray([], dtype=np.float64)

    all_mz = np.asarray(all_mz, dtype=np.float64)

    if mz_res is None:
        mz_res = _estimate_mass_resolution(ms_data=ms_data, units='da', n_sample=N_sample)
        if not np.isfinite(mz_res) or mz_res <= 0:
            mz_res = 0.005  # Default fallback consistent with tolerance estimation

    # Estimate bin number based on resolution; add margin to avoid edge omission
    mz_min = float(np.min(all_mz)) - 5.0 * mz_res
    mz_max = float(np.max(all_mz)) + 5.0 * mz_res
    n_bins = int(max(1, np.round((mz_max - mz_min) / mz_res) + 1))

    # Weights: Each pixel contributes about 1/N_sample (approximating pixel frequency)
    weights = np.zeros_like(all_mz, dtype=np.float64) + (1.0 / float(N_sample))
    hist, bin_edges = np.histogram(all_mz, bins=n_bins, weights=weights)

    # Smooth the histogram using LOESS
    smoothed = _smooth1d(bin_edges, hist, window=smoothing_window, method='loess', weighting='tri-cubic')

    # Detect peaks on the smoothed curve (local maximas)
    ma, _ = find_peaks(smoothed, height=None, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    if ma.size == 0:
        logger.info("No peaks detected on the smoothed histogram, returning empty reference axis.")
        return np.asarray([], dtype=np.float64)

    # Filter peaks by original histogram height threshold (≥ px_perc)
    keep = ma[hist[ma] >= float(px_perc)]
    if keep.size == 0:
        logger.info("No peaks above the pixel frequency threshold (px_perc), returning empty reference axis.")
        return np.asarray([], dtype=np.float64)

    # Take the left bin edge as reference m/z; strictly ascending
    ref = np.asarray(bin_edges[keep], dtype=np.float64)
    if ref.size > 1:
        ref.sort()

    logger.info(f"Histogram/Loess reference axis generated: {ref.size} features (mz_res={mz_res}, px_perc={px_perc}, N_sample={N_sample})")
    return ref

def peak_matching(
        n_sample, # Used for _estimate_mass_resolution
        gap_to_tol_factor,
        ms_data: MS,
        reference_mz_axis: np.ndarray,
        tolerance: Optional[float] = None,
        units: str = 'ppm',
        combiner: str = 'max',
        unique: bool = False,
        method: str = 'diff',
) -> MS:
    """
    Use an optional alignment method to align the peaks of all pixels to the reference m/z axis.
    - method='diff': Difference matching  
    - method='dp': Dynamic programming global matching

    For each pixel's centroid peak x, call the selected matching method to find the reference peak y within the tolerance.
    Then aggregate the intensity of that pixel's centroid peak x into the corresponding reference peak y, using the combiner method.

    Args:
    - ms_data: MS object, containing sparse centroid peak tables for several pixels
    - reference_mz_axis: Reference m/z axis (strictly increasing array)
    - tolerance: Tolerance value; None to estimate a typical value based on the data (_get_tol)
    - units: 'ppm' or 'da', determining the shape of the tolerance (ppm is dynamic with y)
    - combiner: 'max' | 'sum' | 'mean', aggregation method for peaks falling into the same reference peak
    - unique: If True, each reference peak is only matched once per pixel (diff_matching's unique)

    Returns:
    - MS: Aligned data container, with each spectrum's `mz_list` unified to the reference axis, and `intensity` aggregated.
    """
    # Tolerance and validity of units
    units = (units or 'ppm').lower().strip()
    if units not in ('ppm', 'da'):
        raise ValueError("units must be 'ppm' or 'da'")
    combiner = (combiner or 'max').lower().strip()
    if combiner not in ('max', 'sum', 'mean'):
        raise ValueError("combiner must be 'max', 'sum', or 'mean'")
    method = (method or 'diff').lower().strip()
    if method not in ('diff', 'dp'):
        raise ValueError("method must be 'diff' or 'dp'")

    ref = np.asarray(reference_mz_axis, dtype=np.float64)
    if ref.ndim != 1 or ref.size == 0:
        raise ValueError("reference_mz_axis must be a non-empty 1D array")
    if ref.size > 1 and not np.all(ref[1:] > ref[:-1]):
        raise ValueError("reference_mz_axis must be strictly increasing")

    if tolerance is None:
        tolerance = _estimate_mass_resolution(ms_data=ms_data, units=units, n_sample=n_sample, gap_to_tol_factor=gap_to_tol_factor)
        # Sanity check and fallback to default values.
        if not np.isfinite(tolerance) or tolerance <= 0:
            fallback = 10.0 if units == 'ppm' else 0.02
            logger.warning(f"tolerance estimation is invalid, fallback to default: {fallback} {units}")
            tolerance = fallback

    ms_aligned = MS()
    if ms_data.meta is not None:
        ms_aligned.meta = ms_data.meta

    n_features = ref.size
    N = len(ms_data)

    for i in range(N):
        spec = ms_data[i]
        mz_raw = spec.mz_list
        I_raw = spec.intensity

        # None-guarding and empty spectrum handling.
        if mz_raw is None or I_raw is None:
            new_I = np.zeros((n_features,), dtype=np.float64)
            ms_aligned.add_spectrum(SpectrumBaseModule(
                mz_list=ref,
                intensity=new_I,
                coordinates=getattr(spec, 'coordinates', None)
            ))
            continue

        mz = np.asarray(mz_raw, dtype=np.float64)
        I = np.asarray(I_raw, dtype=np.float64)

        # Truncate to the shortest length when m/z or intensity arrays have different lengths.
        if mz.size != I.size:
            m = int(min(mz.size, I.size))
            if m == 0:
                new_I = np.zeros((n_features,), dtype=np.float64)
                ms_aligned.add_spectrum(SpectrumBaseModule(
                    mz_list=ref,
                    intensity=new_I,
                    coordinates=getattr(spec, 'coordinates', None)
                ))
                continue
            mz = mz[:m]
            I = I[:m]

        if mz.size == 0:
            new_I = np.zeros((n_features,), dtype=np.float64)
            ms_aligned.add_spectrum(SpectrumBaseModule(
                mz_list=ref,
                intensity=new_I,
                coordinates=getattr(spec, 'coordinates', None)
            ))
            continue

        # Filter out invalid m/z or non-positive intensity values.
        valid = np.isfinite(mz) & np.isfinite(I) & (I > 0)
        if not np.any(valid):
            new_I = np.zeros((n_features,), dtype=np.float64)
            ms_aligned.add_spectrum(SpectrumBaseModule(
                mz_list=ref,
                intensity=new_I,
                coordinates=getattr(spec, 'coordinates', None)
            ))
            continue
        mz = mz[valid]
        I = I[valid]

        # Use the selected method to match m/z values in x to reference peaks in y.
        if method == 'diff':
            _, _, _, y_idx = _diff_matching(
                x=mz,
                y=ref,
                tolerance=float(tolerance),
                units=units,
                unique=unique,
                return_indices=True
            )
        # Get the indices of matched reference peaks y for each x.
        else:  # method == 'dp'
            _, _, _, y_idx = _dp_matching(
                x=mz,
                y=ref,
                tolerance=float(tolerance),
                units=units,
                unique=unique,
                return_indices=True,
                gap=0.0
            )

        new_I = np.zeros((n_features,), dtype=np.float64)
        if combiner == 'mean':
            counts = np.zeros((n_features,), dtype=np.int32)

        matched = y_idx >= 0
        if np.any(matched):
            j = y_idx[matched]
            v = I[matched]
            if combiner == 'sum' or combiner == 'mean':
                np.add.at(new_I, j, v)
                if combiner == 'mean':
                    np.add.at(counts, j, 1)
            elif combiner == 'max':
                for jj, vv in zip(j, v):
                    if vv > new_I[jj]:
                        new_I[jj] = vv

        if combiner == 'mean':
            new_I = np.divide(new_I, counts, out=new_I, where=counts > 0)

        ms_aligned.add_spectrum(SpectrumBaseModule(
            mz_list=ref,
            intensity=new_I,
            coordinates=spec.coordinates
        ))

    logger.info(f"Peak alignment using {method} completed: {N} pixels, {n_features} reference features (combiner={combiner}, units={units}, tol={tolerance})")
    return ms_aligned


# Internal helper functions
# 1.Compute global mean spectrum
def _compute_mean_spectrum(ms_data: MS,
                       agg: str = 'mean',
                       round_digits: int = 6) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the global mean spectrum for centroided sparse peak tables:
    Each pixel already has (mz, intensity) as centroided sparse peaks. We merge by rounded m/z and then sum or average.

    Parameters:
    - agg: 'mean' or 'sum' controlling aggregation of global intensity; 'mean' averages by pixel occurrence count.
    - round_digits: round m/z to the given decimal places before merging to reduce floating-point jitter (default 6).

    Returns:
    - (global_mz, global_intensity): both are 1D np.ndarray sorted ascending by m/z.
    """

    N_pix = len(ms_data)
    if N_pix == 0:
        logger.warning("MS object is empty, cannot compute global mean spectrum.")
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)

    agg = (agg or 'mean').lower().strip()
    if agg not in ('mean', 'sum'):
        raise ValueError("agg must be 'mean' or 'sum'")

    # Global accumulation containers: use rounded m/z as key; precision controlled by `round_digits`
    sums = {}
    counts = {}

    for i in range(N_pix):
        spec = ms_data[i]
        if spec.mz_list is None or spec.intensity is None:
            continue
        mz = np.asarray(spec.mz_list, dtype=np.float64)
        I = np.asarray(spec.intensity, dtype=np.float64)
        if mz.size == 0:
            continue

        valid = np.isfinite(mz) & np.isfinite(I) & (I > 0)
        if not valid.any():
            continue
        mz = mz[valid]
        I = I[valid]

        # For the current pixel, merge by rounded m/z first (avoid duplicate keys within the pixel); precision controlled by `round_digits`
        mz_r = np.round(mz, round_digits)
        # Sort to enable groupby
        order = np.argsort(mz_r, kind='mergesort')
        mz_r = mz_r[order]
        I = I[order]

        # Unique m/z keys and their summed intensity for the pixel
        unique_keys, idx = np.unique(mz_r, return_inverse=True)
        pixel_sum = np.zeros(unique_keys.size, dtype=np.float64)
        np.add.at(pixel_sum, idx, I)

        # Write into global containers: each pixel contributes one count per key (for mean)
        for k, intensity_sum in zip(unique_keys, pixel_sum):
            if k in sums:
                sums[k] += float(intensity_sum)
                counts[k] += 1
            else:
                sums[k] = float(intensity_sum)
                counts[k] = 1

    if not sums:
        logger.warning("No valid peaks were aggregated to compute global mean spectrum.")
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)

    # Pack into ascending arrays
    keys = np.asarray(sorted(sums.keys()), dtype=np.float64)
    sums_arr = np.asarray([sums[k] for k in keys], dtype=np.float64)
    counts_arr = np.asarray([counts[k] for k in keys], dtype=np.int32)

    if agg == 'sum':
        global_intensity = sums_arr
    else:
        # Average over the number of pixels the key appears in; pixels without the key are excluded from the denominator
        global_intensity = np.divide(sums_arr,
                                     counts_arr,
                                     out=np.zeros_like(sums_arr),
                                     where=counts_arr > 0)

    logger.info(f"Global mean spectrum computed: {keys.size} features (agg={agg})")
    return keys, global_intensity

# 2.Find local peaks from array
def _find_peaks_from_array(
        x: np.ndarray,
        half_window: float = 3,
        findlimits: bool = False,
        snr_threshold: Optional[float] = None,
        noise_method: str = 'mad',
        min_height: Optional[float] = None):
    """
    First find all local maxima, then filter by peak height/noise ratio (SNR).

    Parameters:
    - x: intensity vector (1D)
    - half_window: half-window size (indices) used for local comparison
    - findlimits: whether to return plateau left/right limits to obtain peak width
    - snr_threshold: if provided, keep peaks with SNR >= threshold (SNR = peak height / noise)
    - noise_method: noise estimation method, 'mad' (default, robust) or 'std'
    - min_height: if provided, keep peaks with intensity >= min_height

    Returns:
    - if findlimits=True: (peaks_idx, left_limits, right_limits)
    - if findlimits=False: peaks_idx
    """

    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if n == 0:
        return (np.asarray([], dtype=int), np.asarray([], dtype=int), np.asarray([], dtype=int)) if findlimits else np.asarray([], dtype=int)

    w = int(max(1, round(float(half_window))))

    # Initial: mark points that reach the maximum within their local window
    is_candidate = np.zeros(n, dtype=bool)
    for i in range(n):
        l = max(0, i - w)
        r = min(n, i + w + 1)
        seg = x[l:r]
        if x[i] >= np.max(seg):
            is_candidate[i] = True

    # Noise estimation (robust global estimate)
    eps = 1e-12
    if snr_threshold is not None:
        nm = (noise_method or 'mad').lower().strip()
        if nm == 'std':
            noise = float(np.std(x))
        else:
            med = float(np.median(x))
            mad = float(np.median(np.abs(x - med)))
            noise = 1.4826 * mad
        noise = max(noise, eps)

    # Compress plateaus: keep only the center for adjacent candidates with equal intensity
    peaks_idx = []
    left_limits = []
    right_limits = []

    i = 0
    while i < n:
        if not is_candidate[i]:
            i += 1
            continue
        # Expand plateau region (same intensity and also candidate)
        j = i
        val = x[i]
        while j + 1 < n and is_candidate[j + 1] and x[j + 1] == val:
            j += 1

        # Use plateau center as peak position
        center = (i + j) // 2

        # Filter by SNR and minimum height
        keep = True
        if min_height is not None and val < float(min_height):
            keep = False
        if keep and snr_threshold is not None:
            snr = float(val) / noise
            if snr < float(snr_threshold):
                keep = False

        if keep:
            peaks_idx.append(center)
            if findlimits:
                # Left/right limits: extend left/right until intensity drops (below peak value)
                L = i
                while L - 1 >= 0 and x[L - 1] == val:
                    L -= 1
                R = j
                while R + 1 < n and x[R + 1] == val:
                    R += 1
                left_limits.append(L)
                right_limits.append(R)

        i = j + 1

    peaks_idx = np.asarray(peaks_idx, dtype=int)
    if findlimits:
        return peaks_idx, np.asarray(left_limits, dtype=int), np.asarray(right_limits, dtype=int)
    else:
        return peaks_idx

# 3.Complete diff_match alignment logic
def _diff_matching(
        x: np.ndarray,
        y: np.ndarray,
        tolerance: float,
        units: str = 'ppm',
        unique: bool = False,
        fill_value: float = np.nan,
        return_indices: bool = False):
    """
    Given pixel peaks `x` and reference peaks `y`, within a difference threshold (dynamic ppm or fixed Da),
    find the nearest `y` for each peak in `x`. The output mapping has the same length as `x`, with unmatched positions filled with NaN.

    Parameters:
    - x: 1D np.ndarray of pixel peak m/z (centroided sparse peaks)
    - y: 1D np.ndarray of reference-axis m/z (strictly increasing)
    - tolerance: threshold size; when units='ppm' it's ppm; when units='da' it's Da
    - units: 'ppm' or 'da'; with ppm the threshold depends on candidate reference peak y[j]
    - unique: if True, each reference peak can be matched at most once; if False, multiple `x` may match the same `y`
    - fill_value: value used to fill `aligned_y` at unmatched positions (default NaN)
    - return_indices: if True, additionally return (x_idx, y_idx) index mapping; unmatched y_idx is -1

    Returns:
    - if return_indices=False: (aligned_x, aligned_y)
    - if return_indices=True: (aligned_x, aligned_y, x_indices, y_indices)
    """

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be 1D arrays")
    if y.size == 0:
        # No reference peaks; everything unmatched
        aligned_x = x.copy()
        aligned_y = np.full_like(x, fill_value, dtype=np.float64)
        if return_indices:
            return aligned_x, aligned_y, np.arange(x.size, dtype=int), np.full(x.size, -1, dtype=int)
        return aligned_x, aligned_y

    # Ensure `y` is strictly increasing; sort ascending if not
    if y.size > 1 and not np.all(y[1:] > y[:-1]):
        order_y = np.argsort(y, kind='mergesort')
        y = y[order_y]
    # Keep `x` in original order; we output in original order

    aligned_x = x.copy()
    aligned_y = np.full(x.shape, fill_value, dtype=np.float64)
    x_indices = np.arange(x.size, dtype=int)
    y_indices = np.full(x.size, -1, dtype=int)

    units = (units or 'ppm').lower().strip()
    if units not in ('ppm', 'da'):
        raise ValueError("units must be 'ppm' or 'da'")

    # Track whether a reference peak has been used (effective when unique=True)
    used_y = np.zeros(y.size, dtype=bool)

    # For each `x` peak, find nearest candidates in `y` (both sides of insertion index)
    for i in range(x.size):
        xi = float(x[i])
        pos = int(np.searchsorted(y, xi, side='left')) # binary search
        candidates = []
        if pos < y.size:
            candidates.append(pos)
        if pos - 1 >= 0:
            candidates.append(pos - 1)

        best_j = -1
        best_diff = np.inf
        for j in candidates:
            yj = float(y[j])
            diff = abs(xi - yj)
            tol_da = tolerance * 1e-6 * yj if units == 'ppm' else float(tolerance)
            if diff <= tol_da:
                # Candidate satisfies tolerance; choose nearest; skip used y when unique=True
                if unique and used_y[j]:
                    continue
                if diff < best_diff:
                    best_diff = diff
                    best_j = j

        if best_j >= 0:
            aligned_y[i] = y[best_j]
            y_indices[i] = best_j
            if unique:
                used_y[best_j] = True

    if return_indices:
        return aligned_x, aligned_y, x_indices, y_indices
    return aligned_x, aligned_y

#4.DP Global optimal alignment under monotonic constraints via dynamic programming.
def _dp_matching(x: np.ndarray,
                 y: np.ndarray,
                 tolerance: float,
                 units: str = 'ppm',
                 unique: bool = False,
                 fill_value: float = np.nan,
                 return_indices: bool = False,
                 gap: float = 0.0):
    """
    Performs global optimal alignment under monotonic constraints using Dynamic Programming (DP).

    - score_function：score(x_i, y_j) = 1 / (1 + |x_i - y_j|)
    - gap penalty: Default 0 (consistent with user description), meaning skipping x or y is not penalized.
    - tolerance: Only matches when |x_i - y_j| <= tol_da; otherwise, considered non-matching.

    - If `return_indices=False`: returns `(aligned_x, aligned_y)`
    - If `return_indices=True`: returns `(aligned_x, aligned_y, x_indices, y_indices)`
    - `y_indices` contains the matched reference index for each x; unmatched positions are marked with `-1`.

    Notes:
    - The DP algorithm inherently ensures one-to-one, monotonic, and non-crossing matches; therefore, the `unique` parameter is always respected (its boolean value is effectively ignored).
    - To ensure consistent behavior with the `diff` method, ppm/Da tolerance conversion and gating are still performed here.
    
    """

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be 1D arrays")

    # reference is empty: all peaks unmatched
    if y.size == 0:
        aligned_x = x.copy()
        aligned_y = np.full_like(x, fill_value, dtype=np.float64)
        if return_indices:
            return aligned_x, aligned_y, np.arange(x.size, dtype=int), np.full(x.size, -1, dtype=int)
        return aligned_x, aligned_y

    # ensure y is strictly increasing
    if y.size > 1 and not np.all(y[1:] > y[:-1]):
        order_y = np.argsort(y, kind='mergesort')
        y = y[order_y]

    # DP assumes x is sorted ascending to ensure monotonic matches; keep original order for output
    order_x = np.argsort(x, kind='mergesort')
    x_sorted = x[order_x]

    n = int(x_sorted.size)
    m = int(y.size)

    units = (units or 'ppm').lower().strip()
    if units not in ('ppm', 'da'):
        raise ValueError("units must be 'ppm' or 'da'")

    # Each reference peak's tolerance (da space); ppm is dynamic with y
    if units == 'ppm':
        tol_vec = (float(tolerance) * 1e-6) * y
    else:
        tol_vec = np.full((m,), float(tolerance), dtype=np.float64)

    # DP table and backtracking pointer
    # Pointer encoding: 1=match(up-left diag), 2=skip x(up), 3=skip y(left)
    dp = np.zeros((n + 1, m + 1), dtype=np.float64)
    ptr = np.zeros((n + 1, m + 1), dtype=np.int8)

    # Optional gap penalty (default 0)
    gap = float(gap or 0.0)

    # Dynamic programming table filling
    for i in range(1, n + 1):
        xi = float(x_sorted[i - 1])
        for j in range(1, m + 1):
            yj = float(y[j - 1])
            diff = abs(xi - yj)
            allowed = diff <= tol_vec[j - 1]
            # Match score: 1 / (1 + diff) if allowed, else -inf
            match_score = (1.0 / (1.0 + diff)) if allowed else -np.inf

            # Three transitions
            up = dp[i - 1, j] - gap          # Skip x_i
            left = dp[i, j - 1] - gap        # Skip y_j
            diag = dp[i - 1, j - 1] + match_score  # Match

            # Choose best
            # Encourage matches by preferring diag in ties
            if diag >= up and diag >= left:
                dp[i, j] = diag
                ptr[i, j] = 1
            elif up >= left:
                dp[i, j] = up
                ptr[i, j] = 2
            else:
                dp[i, j] = left
                ptr[i, j] = 3

    # Backtracking to obtain matched pairs
    y_indices_sorted = np.full((n,), -1, dtype=int)
    i, j = n, m
    while i > 0 and j > 0:
        step = int(ptr[i, j])
        if step == 1:
            # Match (i-1, j-1); only chosen when allowed
            y_indices_sorted[i - 1] = j - 1
            i -= 1
            j -= 1
        elif step == 2:
            # Skip x_i
            i -= 1
        elif step == 3:
            # Skip y_j
            j -= 1
        else:
            # Unset pointer (n=0 or m=0 etc.); prefer reducing the larger index
            if i > j:
                i -= 1
            else:
                j -= 1

    # Restore to original x order
    inv_order = np.empty_like(order_x)
    inv_order[order_x] = np.arange(n, dtype=int)
    y_indices = y_indices_sorted[inv_order]

    aligned_x = x.copy()
    aligned_y = np.full_like(x, fill_value, dtype=np.float64)
    matched = y_indices >= 0
    if np.any(matched):
        aligned_y[matched] = y[y_indices[matched]]

    if return_indices:
        x_indices = np.arange(x.size, dtype=int)
        return aligned_x, aligned_y, x_indices, y_indices
    return aligned_x, aligned_y


# Utility: estimate alignment tolerance; compute resolution; return tol
def _estimate_mass_resolution(
    ms_data: MS,
    units: Literal['ppm', 'da'] = 'ppm',
    method: Literal['median', 'mad'] = 'median',
    n_sample: int = 2000,
    gap_to_tol_factor: float = 0.5,
    random_state: Optional[int] = 42
) -> float:
    """
    Automatically estimate alignment tolerance (half-width).

    When units='ppm', first convert gaps to ppm: ppm_gap = gap / mid_mz * 1e6.
    When units='da', directly use gaps in Da.
    Returns a tolerance constant (half-width) for alignment, not the resolution; if you need resolution, use ~2 * tolerance.

    Parameters:
    - ms_data: MS object
    - units: 'ppm' or 'da'
    - method: 'median' (median) or 'mad' (median absolute deviation)
    - n_sample: number of spectra to sample
    - random_state: RNG seed; default 42 for reproducibility (set None to disable)

    Returns:
    - Estimated tolerance (half-width), unit determined by `units`.
    """
    N = len(ms_data)
    logger.info(
        f"Alignment tolerance estimation parameters: units={units}, method={method}, n_sample={n_sample}, "
        f"gap_to_tol_factor={gap_to_tol_factor}, random_state={random_state}"
    )
    if N == 0:
        logger.warning("MS object is empty, using default tolerance")
        return 5.0 if units == 'ppm' else 0.01

    n_sample = min(max(1, n_sample), N)
    rng = np.random.default_rng(random_state)
    sample_idx = rng.choice(N, size=n_sample, replace=False)

    all_vals = []  # Collect ppm_gap or da_gap

    for idx in sample_idx:
        spec = ms_data[int(idx)]
        mz = spec.mz_list

        if mz is None or np.size(mz) <= 1:
            continue

        mz = np.asarray(mz, dtype=np.float64)
        # Ensure sorted order
        if mz.size > 1 and not np.all(mz[1:] >= mz[:-1]):
            mz = np.sort(mz)

        gaps = np.diff(mz)
        if gaps.size == 0:
            continue

        # Convert gaps to target units for robust statistics
        if units == 'ppm':
            mid_mz = (mz[:-1] + mz[1:]) / 2.0
            vals = gaps / np.maximum(mid_mz, 1e-9) * 1e6
        elif units == 'da':
            vals = gaps
        else:
            raise ValueError("units must be either 'ppm' or 'da'")

        # IQR filter to remove outliers for a robust typical gap
        q1, q3 = np.percentile(vals, [25, 75])
        iqr = q3 - q1
        lower = max(1e-9, q1 - 1.5 * iqr)
        upper = q3 + 1.5 * iqr
        valid = vals[(vals >= lower) & (vals <= upper)]
        if valid.size:
            all_vals.extend(valid)

    if not all_vals:
        logger.warning("No valid gaps found, using default tolerance")
        return 10.0 if units == 'ppm' else 0.02

    all_vals = np.asarray(all_vals, dtype=np.float64)
    if method == 'mad':
        from scipy.stats import median_abs_deviation
        est = float(median_abs_deviation(all_vals, scale='normal'))
    else:
        est = float(np.median(all_vals))

    # Return half-width (tolerance); conventionally resolution ≈ 2 * tolerance
    tol = est * gap_to_tol_factor
    logger.info(f"Estimated alignment tolerance (half-width): {tol:.6f} {units}")
    return tol

#6. Simple smoothing algorithm for 1D data
def _smooth1d(x: np.ndarray,
              y: np.ndarray,
              window: int = 10,
              method: str = "loess",
              weighting: str = "tri-cubic") -> np.ndarray:
    """
    Quick smoothing of evenly spaced data (simplified moving LOESS) for histogram smoothing.

    Parameters:
    - x: Evenly spaced independent variable (pass bin_edges[:-1] for each bin's left boundary)
    - y: Signal array (histogram heights)
    - window: Sliding window length (default 10)
    - method: "loess" | "lowess" | "average" (default "loess")
    - weighting: "tri-cubic" | "gaussian" | "linear" (default "tri-cubic")

    Returns:
    - yhat: Smoothed signal, same length as y
    """
    from scipy import signal, linalg

    y = np.asarray(y, dtype=float)
    leny = y.size
    if leny == 0:
        return y

    halfw = int(np.floor(float(window) / 2.0))
    window = int(2.0 * halfw + 1.0)
    x1 = np.arange(1.0 - halfw, (halfw - 1.0) + 1.0, dtype=float)

    if weighting == 'tri-cubic':
        weight = (1.0 - (np.abs(x1) / halfw) ** 3.0) ** 1.5
    elif weighting == 'gaussian':
        weight = np.exp(- (x1 / halfw * 2.0) ** 2.0)
    elif weighting == 'linear':
        weight = 1.0 - (np.abs(x1) / halfw)
    else:
        weight = (1.0 - (np.abs(x1) / halfw) ** 3.0) ** 1.5

    if method == 'loess':
        V = np.vstack([weight, weight * x1, weight * x1 * x1]).T
        order = 2
    elif method == 'lowess':
        V = np.vstack([weight, weight * x1]).T
        order = 1
    else:
        V = weight.reshape(-1, 1)
        order = 0

    Q, R = linalg.qr(V, mode='economic')
    alpha = Q[halfw - 1] @ Q.T

    yhat = signal.lfilter(alpha * weight, 1.0, y)
    yhat[int(halfw + 1) - 1:-halfw] = yhat[int(window - 1) - 1:-1]

    x2 = np.arange(1.0, (window - 1.0) + 1.0, dtype=float)
    if method == 'loess':
        V2 = np.vstack([np.ones(window - 1), x2, x2 * x2]).T
    elif method == 'lowess':
        V2 = np.vstack([np.ones(window - 1), x2]).T
    else:
        V2 = np.ones((window - 1, 1))

    for j in range(1, halfw + 1):
        if weighting == 'tri-cubic':
            wj = (1.0 - (np.abs(np.arange(1, window) - j) / (window - j)) ** 3.0) ** 1.5
        elif weighting == 'gaussian':
            wj = np.exp(- (np.abs(np.arange(1, window) - j) / (window - j) * 2.0) ** 2.0)
        else:
            wj = 1.0 - (np.abs(np.arange(1, window) - j) / (window - j))

        W = np.kron(np.ones((order + 1, 1)), wj).T
        Q2, R2 = linalg.qr(V2 * W, mode='economic')
        alpha2 = Q2[j - 1] @ Q2.T
        alpha2 = alpha2 * wj
        yhat[int(j) - 1] = alpha2 @ y[:window - 1]
        yhat[int(-j)] = alpha2 @ y[np.arange(leny - 1, leny - window, -1, dtype=int)]

    return yhat