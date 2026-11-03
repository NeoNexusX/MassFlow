from typing import Optional
import numpy as np
from numba import jit, prange
from massflow.preprocess.numba.peak_align_numba import search_nearest_jit
from massflow.tools.funs import prepare_flat_inputs, lengths_to_offsets
from massflow.tools.logger import get_logger

logger = get_logger("massflow.preprocess.normalization")

# -----------------------------------------------------------------------------
# Core Calculation Kernels (Helpers)
# -----------------------------------------------------------------------------


@jit(nopython=True, cache=True)
def _get_tic(arr: np.ndarray) -> float:
    return np.sum(arr)

@jit(nopython=True, cache=True)
def _get_rms(arr: np.ndarray) -> float:
    # RMS = sqrt( sum(x^2) / n )
    s = 0.0
    n = arr.shape[0]
    if n == 0:
        return 0.0
    for i in range(n):
        val = arr[i]
        s += val * val
    return np.sqrt(s / n)

@jit(nopython=True, cache=True)
def _get_ref_idx(
    mz_arr: np.ndarray,
    ref: float,
    ref_tolerance: float,
) -> int:
    """Return the index of the m/z point matching ref within tolerance, or -1."""
    if mz_arr.shape[0] == 0:
        return -1

    query = np.array([ref], dtype=np.float64)
    idx_arr = search_nearest_jit(query, mz_arr, ref_tolerance, 2, -1, False)
    return int(idx_arr[0])


# -----------------------------------------------------------------------------
# flat Logic
# -----------------------------------------------------------------------------

@jit(nopython=True, cache=True)
def _normalize_core_basic(
    in_arr: np.ndarray,
    out_arr: np.ndarray,
    method_enum: int,
    scale_factor: float,
):
    """Normalize one spectrum using TIC/RMS base."""
    if method_enum == 0:
        base = _get_tic(in_arr)
    else:
        base = _get_rms(in_arr)

    if base <= 0.0 or np.isnan(base):
        return

    factor = scale_factor / base
    n = in_arr.shape[0]
    for i in range(n):
        out_arr[i] = np.float64(in_arr[i]) * factor


@jit(nopython=True, cache=True)
def _normalize_core_ref(
    in_arr: np.ndarray,
    out_arr: np.ndarray,
    mz_arr: np.ndarray,
    scale_factor: float,
    ref: float,
    ref_tolerance: float,
) -> int:
    """Normalize one spectrum using ref peak intensity base.

    Returns 0 when normalized, 1 when no m/z matches ref within tolerance,
    2 when the matched ref intensity is non-positive or NaN. On failure the
    input spectrum is copied through unchanged instead of zeroed (BUG-09).
    """
    idx = _get_ref_idx(mz_arr, ref, ref_tolerance)
    if idx < 0:
        for i in range(in_arr.shape[0]):
            out_arr[i] = in_arr[i]
        return 1

    base = np.float64(in_arr[idx])
    if base <= 0.0 or np.isnan(base):
        for i in range(in_arr.shape[0]):
            out_arr[i] = in_arr[i]
        return 2

    factor = scale_factor / base
    n = in_arr.shape[0]
    for i in range(n):
        out_arr[i] = np.float64(in_arr[i]) * factor
    return 0

@jit(nopython=True, parallel=True, cache=True)
def _normalize_flat_basic_jit(
    flat: np.ndarray,
    method_enum: int,
    scale_factors: np.ndarray,
    lengths: np.ndarray,
) -> np.ndarray:
    res = np.zeros(flat.size, dtype=flat.dtype)
    offsets = lengths_to_offsets(lengths)

    for p in prange(lengths.size): # pylint: disable=not-an-iterable
        start = offsets[p]
        end = offsets[p + 1]
        if end > start:
            _normalize_core_basic(
                flat[start:end],
                res[start:end],
                method_enum,
                np.float64(scale_factors[p]),
            )

    return res


@jit(nopython=True, parallel=True, cache=True)
def _normalize_flat_ref_jit(
    flat: np.ndarray,
    mz_flat: np.ndarray,
    scale_factors: np.ndarray,
    ref: float,
    ref_tolerance: float,
    lengths: np.ndarray,
    mz_shared: bool,
):
    """Normalize flat spectra in ref mode; supports per-spectrum or shared mz axis.

    Returns (result, status): status[p] is 0 when spectrum p was normalized,
    1 when no m/z matched ref within tolerance, 2 when the matched ref
    intensity was non-positive or NaN; both failure modes keep the original
    spectrum values.
    """
    res = np.zeros(flat.size, dtype=flat.dtype)
    status = np.zeros(lengths.size, dtype=np.int64)
    offsets = lengths_to_offsets(lengths)

    for p in prange(lengths.size): # pylint: disable=not-an-iterable
        start = offsets[p]
        end = offsets[p + 1]
        if end > start:
            if mz_shared:
                mz_arr = mz_flat
            else:
                mz_arr = mz_flat[start:end]
            status[p] = _normalize_core_ref(
                flat[start:end],
                res[start:end],
                mz_arr,
                np.float64(scale_factors[p]),
                ref,
                ref_tolerance,
            )

    return res, status


# -----------------------------------------------------------------------------
# Main Dispatcher
# -----------------------------------------------------------------------------

def normalizer_numba(
    intensity: np.ndarray,
    method: str = "tic",
    scale: Optional[float] = None,
    mz_flat: Optional[np.ndarray] = None,
    ref: Optional[float] = None,
    ref_tolerance: float = 0.1,
    lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Numba backend for flat-mode normalization.

    Parameters:
        intensity: 1D flat array.
        method: 'tic', 'rms', 'ref'.
        scale: Amplitude scaling factor. If None, default is per-spectrum length for
            tic/rms and 1.0 for ref mode.
        mz_flat: Optional 1D flat mz array; required when method='ref'.
            Also supports shared mz axis (1 spectrum-length array) for
            continuous data where each spectrum uses the same mz grid.
        ref: Optional reference mz value; required when method='ref'.
        ref_tolerance: Tolerance for ref matching. Default 0.1.
        lengths: (Optional) 1D int array of valid lengths per spectrum. If omitted,
            the full input is treated as one spectrum.
    """

    _, intensity_arr, lengths_arr = prepare_flat_inputs(None, intensity, lengths)

    # Map Method
    method_clean = method.strip().lower()
    if method_clean == "tic":
        m_enum = 0
    elif method_clean == "rms":
        m_enum = 1
    elif method_clean == "ref":
        m_enum = 2
    else:
        raise ValueError(f"Unknown method for numba normalizer: {method}")

    if not np.isfinite(ref_tolerance) or float(ref_tolerance) < 0.0:
        raise ValueError("ref_tolerance must be a finite non-negative number")

    n_spec = lengths_arr.size

    if scale is None:
        if m_enum == 2:
            # Ref mode defaults to 1.0 for each spectrum.
            scale_factors = np.ones(n_spec, dtype=np.float64)
        else:
            # TIC/RMS defaults to per-spectrum length, not the flat length.
            scale_factors = lengths_arr.astype(np.float64)
    else:
        if not np.isfinite(scale) or float(scale) < 0.0:
            raise ValueError("scale must be a finite non-negative number")
        scale_factors = np.full(n_spec, float(scale), dtype=np.float64)

    if m_enum == 2:
        if mz_flat is None:
            raise ValueError("mz_flat is required when method='ref'")
        if ref is None:
            raise ValueError("ref is required when method='ref'")
        mz_arr = np.asarray(mz_flat, dtype=np.float64)
        if mz_arr.ndim != 1:
            raise ValueError("mz_flat must be a 1D array")

        mz_shared = mz_arr.size != intensity_arr.size
        if mz_shared:
            if intensity_arr.size % mz_arr.size != 0:
                raise ValueError("mz_flat must be a 1D array with the same length as intensity")
            if np.any(lengths_arr != mz_arr.size):
                raise ValueError(
                    "Shared mz_flat length must match every spectrum length when using ref normalization"
                )

        result, status = _normalize_flat_ref_jit(
            intensity_arr,
            mz_arr,
            scale_factors,
            float(ref),
            float(ref_tolerance),
            lengths_arr,
            mz_shared,
        )

        # Spectra that could not be ref-normalized keep their original
        # values; surface how many so the data loss is visible (BUG-09).
        n_nomatch = int(np.count_nonzero(status == 1))
        n_bad_base = int(np.count_nonzero(status == 2))
        if n_nomatch > 0 or n_bad_base > 0:
            logger.warning(
                f"ref normalization kept {n_nomatch}/{n_spec} spectra unchanged "
                f"(no m/z within {float(ref_tolerance)} of ref={float(ref)}) and "
                f"{n_bad_base}/{n_spec} unchanged (non-positive ref intensity); "
                f"original values preserved"
            )

        return result

    return _normalize_flat_basic_jit(
        intensity_arr,
        m_enum,
        scale_factors,
        lengths_arr,
    )
