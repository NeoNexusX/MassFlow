from __future__ import annotations
import h5py
import os
import gc
from massflow.logger import get_logger
import copy
from dataclasses import dataclass
from typing import List, Optional, Literal, Tuple, Dict
import numpy as np
from numpy.typing import NDArray
from scipy import interpolate, signal, linalg
from massflow.module.ms_module import MS, SpectrumBaseModule

logger = get_logger("peak_alignment")

@dataclass
class PeakBins:
    """Aggregated peak bins

    Fields:
    - peaks: mean peak position per bin (used for subsequent merging)
    - values: mean intensity per bin (used to fill the aligned matrix)
    - counts: number of spectra hitting each bin (for weighting and frequency)
    - tolerance: half-window tolerance used during binning
    - domain: reference domain involved in aggregation (ascending)
    Edge case: bins with counts == 0 have peaks/values set to NaN.
    """
    peaks: NDArray[np.float64]
    values: NDArray[np.float64]
    counts: NDArray[np.int64]
    tolerance: float
    domain: NDArray[np.float64]

@dataclass
class AlignResult:
    """Alignment result object

    Fields:
    - ms_aligned: aligned MS collection with spectra built on the shared reference axis
    - ref: reference peak vector (may come from binpeaks + mergepeaks)
    - count: hit count per reference peak
    - freq: hit frequency (count/nspec)
    - tolerance/units/binfun/binratio: parameter echo for auditing and reproducibility
    """
    ms_aligned: MS
    ref: NDArray[np.float64]
    count: NDArray[np.int64]
    freq: NDArray[np.float64]
    tolerance: float
    units: Literal["relative", "absolute"]
    binfun: Literal["median", "min", "max", "mean"]
    binratio: int

class SpectrumHDF5Lazy(SpectrumBaseModule):
    """
    HDF5-backed spectrum that loads data only when accessed.
    Compatible with SpectrumBaseModule interface.
    """
    def __init__(self, filepath, index, coordinates, mz_ref=None):
        super().__init__(mz_list=None, intensity=None, coordinates=coordinates)
        self._filepath = filepath
        self._index = index
        self._mz_ref = mz_ref

    @property
    def mz_list(self):
        if self._mz_ref is not None:
            return self._mz_ref
        
        with h5py.File(self._filepath, 'r') as f:
            return f['mz'][:]

    @property
    def intensity(self):
        with h5py.File(self._filepath, 'r') as f:
            return f['intensity'][self._index]
            
    def unload(self):
        pass

def _rel_diff(x: NDArray[np.float64], y: Optional[NDArray[np.float64]] = None, ref: Literal["x", "y", "abs"] = "y") -> \
NDArray[np.float64]:
    """Relative/absolute difference computation

    Behavior:
    - If `y` is provided, return element-wise differences using `ref` to choose the denominator.
    - If `y` is None, compute adjacent differences of `x`.

    Parameters:
    - x: float array; acts as the right-hand sequence when `y` is absent
    - y: float array or None; when None, defaults to `x[:-1]`
    - ref: reference for normalization ("x" divide by x, "y" divide by y, "abs" absolute difference)

    Returns:
    - float array of broadcasted length

    Notes:
    - Types are coerced (int→float) consistent with R behavior
    - Division by zero yields `inf`; comparisons `<= tol` treat `inf` as no-hit (consistent with C++)
    """
    if y is None:
        if x.size <= 1:
            return np.array([], dtype=np.float64)
        y = x[:-1]
        x = x[1:]
    # Cast to float64, compatible with NumPy 2.0 (np.float_ removed)
    x = x.astype(np.float64, copy=False)
    y = y.astype(np.float64, copy=False)
    n = max(x.size, y.size)
    xx = np.resize(x, n)
    yy = np.resize(y, n)
    if ref == "x":
        return (xx - yy) / xx
    if ref == "y":
        return (xx - yy) / yy
    return xx - yy

def _udiff_scalar(x: float, y: float, tol_ref: Literal["x", "y", "abs"]) -> float:
    if tol_ref == "x":
        return abs((x - y) / x)
    if tol_ref == "y":
        return abs((x - y) / y)
    return abs(x - y)

def _binary_search(
        x: NDArray[np.float64],
        table: NDArray[np.float64],
        tol: float = 0.0,
        tol_ref: Literal["x", "y", "abs"] = "abs",
        nomatch: int = -1,
        nearest: bool = False,
) -> NDArray[np.int64]:
    """Approximate binary search

    Behavior:
    - In a sorted `table`, for each `x[i]`, compare left/right neighbors by relative/absolute difference.
      If within tolerance or `nearest=True`, return the nearer neighbor (ties choose left).

    Parameters:
    - x: query values
    - table: sorted candidate array (ascending)
    - tol: tolerance half-window
    - tol_ref: difference reference ("x"/"y"/"abs")
    - nomatch: index returned on no-hit (default -1)
    - nearest: when True, allow nearest neighbor even outside tolerance; when False, only return within tolerance

    Returns:
    - int index array pointing into `table`; no-hits are `nomatch`

    Edge cases:
    - `table` must be ascending; otherwise raise (consistent with R)
    - NaN queries are treated as no-hit (`nomatch`)
    Complexity: O(n) traversal + O(log m) search (np.searchsorted), overall O(n log m).

    Raises:
    - ValueError: When `table` is not sorted in ascending order.
    """
    if np.any(np.diff(table) < 0):
        logger.error("_binary_search: 'table' must be sorted")
        raise ValueError("'table' must be sorted")
    pos = np.full(x.shape, nomatch, dtype=np.int64)
    for i, xi in enumerate(x):
        idx = int(np.searchsorted(table, xi))
        i_left = idx - 1
        i_right = idx
        cand_i = None
        if i_left >= 0 and i_right < table.size:
            if xi == table[i_left]:
                pos[i] = i_left
                continue
            if xi == table[i_right]:
                pos[i] = i_right
                continue
            di = _udiff_scalar(xi, table[i_left], tol_ref)
            dj = _udiff_scalar(xi, table[i_right], tol_ref)
            if di <= dj and (nearest or di <= tol):
                cand_i = i_left
            elif dj <= di and (nearest or dj <= tol):
                cand_i = i_right
            else:
                cand_i = nomatch
        elif i_left >= 0:
            di = _udiff_scalar(xi, table[i_left], tol_ref)
            cand_i = i_left if (nearest or di <= tol) else nomatch
        elif i_right < table.size:
            dj = _udiff_scalar(xi, table[i_right], tol_ref)
            cand_i = i_right if (nearest or dj <= tol) else nomatch
        else:
            cand_i = nomatch
        pos[i] = cand_i
    return pos

def _rel_sequence(from_: float, to: float, by: float) -> NDArray[np.float64]:
    """Relative half-width sequence

    Principle:
    - With half-width `half=by/2`, construct geometric ratio `ratio=(1+half)/(1-half)` and compute the length by log scale;
      return the sequence `from * ratio^(i-1)`.

    Parameters: from_, to, by for start, end, and relative step (2×half-width).
    Returns: ascending float array.
    Complexity: O(L) where L is output length; more stable for wide relative grids.
    """
    half = by / 2.0
    ratio = (1.0 + half) / (1.0 - half)
    length_out = int(np.floor(1.0 + (np.log(to) - np.log(from_)) / np.log(ratio)))
    i = np.arange(1, length_out + 1, dtype=np.float64)
    return from_ * np.power(ratio, i - 1.0)

def _mad(x: NDArray[np.float64], constant: float = 1.4826) -> float:
    if x.size == 0:
        return np.nan
    med = float(np.median(x))
    return constant * float(np.median(np.abs(x - med)))

def _es_resolution(x: NDArray[np.float64], tol: Optional[float] = None, ref: Optional[Literal["x", "abs"]] = None) -> float:
    """Resolution estimation

    Behavior:
    - Choose reference (abs or x) by robustness (MAD) between adjacent absolute differences `dx` and relative differences `rx`, and return the minimal difference as resolution.
    - If `tol` is provided, validate grid consistency by checking modulo residuals against `tol`.

    Parameters:
    - x: float array (sorted internally)
    - tol: optional tolerance threshold for grid consistency check
    - ref: reference selector (auto when None)

    Returns: float resolution; NaN when no valid differences are available.
    Numerical stability: filter out values smaller than machine epsilon to avoid divide-by-zero and instability.
    """
    if x.size <= 1:
        return np.nan
    xs = np.sort(x)
    from_ = xs[:-1]
    to = xs[1:]
    rx = 2.0 * ((to / from_) - 1.0) / ((to / from_) + 1.0)
    rx = np.where(np.isnan(rx), np.inf, rx)
    dx = np.diff(xs)
    eps = np.finfo(float).eps
    dx = dx[dx > eps]
    rx = rx[rx > eps]
    chosen = ref
    if chosen is None:
        if dx.size and rx.size:
            lhs = _mad(dx / float(np.max(xs) - np.min(xs)))
            rhs = _mad(rx)
            chosen = "abs" if lhs < rhs else "x"
        elif dx.size:
            chosen = "abs"
        elif rx.size:
            chosen = "x"
        else:
            return np.nan
    if chosen == "abs":
        res = float(np.min(dx)) if dx.size else np.inf
        if tol is None or np.all(np.mod(dx, res) <= tol):
            return res
        return np.nan
    else:
        res = float(np.min(rx)) if rx.size else np.inf
        if tol is None or np.all(np.mod(rx, res) <= tol):
            return res
        return np.nan

def estimate_domain(
        xlist: List[NDArray[np.float64]],
        width: Literal["median", "min", "max", "mean"] = "median",
        units: Literal["relative", "absolute"] = "relative",
) -> Tuple[NDArray[np.float64], float]:
    """Shared reference domain estimation

    Behavior:
    - For each spectrum, estimate resolution and range of peak indices; aggregate resolution by `width`.
    - Generate an ascending reference domain according to units (relative or absolute), applying minimum resolution thresholds.

    Parameters:
    - xlist: list of peak indices per spectrum
    - width: aggregation function ("median"/"min"/"max"/"mean")
    - units: units ("relative"/"absolute")

    Returns:
    - (domain, by): reference domain and resolution step

    Edge cases and numbers:
    - When no peaks are available, return NaN step and empty domain; minimum resolution (relative 5e-7, absolute 1e-4).
    Complexity: O(S·P): per-spectrum O(P) estimation plus a single aggregation.
    """
    ref = "x" if units == "relative" else "abs"
    stats = []
    for x in xlist:
        x = x[np.logical_not(np.isnan(x))]
        res = _es_resolution(x, ref=ref)
        if x.size > 0 and (np.isfinite(res)):
            stats.append((float(np.min(x)), float(np.max(x)), res))
        else:
            stats.append((np.nan, np.nan, res))
    arr = np.array(stats, dtype=np.float64)
    from_ = float(np.floor(np.nanmin(arr[:, 0])))
    to = float(np.ceil(np.nanmax(arr[:, 1])))
    by_values = arr[:, 2]
    if width == "median":
        by = float(np.nanmedian(by_values))
    elif width == "min":
        by = float(np.nanmin(by_values))
    elif width == "max":
        by = float(np.nanmax(by_values))
    else:
        by = float(np.nanmean(by_values))
    min_relative_res = 5e-7
    min_absolute_res = 1e-4
    # ensure 'by' is not too small
    if units == "relative":
        by = max(min_relative_res, round(2.0 * by, 6) * 0.5)
        seq = _rel_sequence(from_, to, by)
    else:
        by = max(min_absolute_res, round(by, 4))
        seq = np.arange(from_, to + by, by, dtype=np.float64)
    return seq, by

def binpeaks(
        peaklist: List[NDArray[np.float64]],
        domain: Optional[NDArray[np.float64]] = None,
        xlist: Optional[List[NDArray[np.float64]]] = None,
        tol: Optional[float] = None,
        tol_ref: Literal["abs", "x", "y"] = "abs",
        merge: bool = False,
        na_drop: bool = True,
) -> PeakBins:
    """Aggregate peaks to the reference domain

    Two-stage behavior:
    1) Initial mapping: per spectrum, use `_binary_search(nearest=False)` to map peak indices into `domain` bins.
    2) De-duplicate conflicts: explicitly detect duplicate hits to the same bin within a spectrum (`while(any(dup))`),
       compute relative/absolute differences to the bin center via `_rel_diff`, keep the closest and mark others as no-hit.
    Finally, mean-aggregate hits to obtain positions and intensities, and record counts.

    Parameters:
    - peaklist: list of peak positions per spectrum
    - domain: reference domain (None → auto-estimate range and step)
    - xlist: list of intensities per spectrum (defaults to positions)
    - tol/tol_ref: tolerance half-window and reference; when tol is None, estimate (median of 0.5×within-spectrum minimal differences)
    - merge: whether to call `mergepeaks` to merge neighboring bins after aggregation
    - na_drop: whether to drop bins with NaN values

    Returns: `PeakBins` containing aggregated positions, intensities, counts, and domain.

    Edge cases:
    - `peaklist[i]` and `xlist[i]` must have equal length; otherwise raise.
    - `domain` must be ascending; otherwise raise.
    Complexity: per spectrum O(P log B) initial mapping + O(P) de-duplication; overall aligns with R implementation.

    Raises:
    - ValueError: When input lengths mismatch or `domain` is not sorted ascending.

    Notes:
    - In relative units (ppm), the effective step is derived dynamically from tolerance semantics to match half-width scaling.
    """
    if xlist is None:
        xlist = peaklist
    if any(len(peaklist[i]) != len(xlist[i]) for i in range(len(peaklist))):
        logger.error("binpeaks: lengths of 'peaklist' and 'xlist' must match")
        raise ValueError("lengths of 'peaklist' and 'xlist' must match")
    if tol is None:
        ref = "abs" if tol_ref == "abs" else "y"
        mins = []
        for peaks in peaklist:
            d = _rel_diff(np.array(peaks, dtype=np.float64), ref=ref)
            if d.size:
                mins.append(float(np.nanmin(d)))
        tol = 0.5 * float(np.median(mins)) if mins else np.nan
    if domain is None:
        lims = []
        for peaks in peaklist:
            if peaks.size:
                lims.append(float(np.min(peaks)))
                lims.append(float(np.max(peaks)))
        if not lims:
            domain = np.array([], dtype=np.float64)
        else:
            lo, hi = float(np.min(lims)), float(np.max(lims))
            if tol_ref == "abs":
                step = tol
                n = int(np.floor((hi - lo) / step)) + 1
                domain = lo + step * np.arange(n, dtype=np.float64)
            else:
                domain = _rel_sequence(lo, hi, tol)
    if domain.size and np.any(np.diff(domain) < 0):
        logger.error("binpeaks: 'domain' must be sorted")
        raise ValueError("'domain' must be sorted")
    peaks_acc = np.zeros(domain.size, dtype=np.float64)
    x_acc = np.zeros(domain.size, dtype=np.float64)
    counts = np.zeros(domain.size, dtype=np.int64)
    for s in range(len(peaklist)):
        p = _binary_search(np.array(peaklist[s], dtype=np.float64), domain, tol, tol_ref, nomatch=-1, nearest=False)
        pos = np.array(peaklist[s], dtype=np.float64)
        vals = np.array(xlist[s], dtype=np.float64)
        # while(any(dup)) explicit loop per R implementation
        # dup marks later occurrences of the same bin index (excluding NA/nomatch)
        while True:
            dup = np.zeros(p.shape, dtype=bool)
            seen_bins: Dict[int, int] = {}
            for i, pi in enumerate(p):
                if pi < 0:
                    continue
                if int(pi) in seen_bins:
                    dup[i] = True
                else:
                    seen_bins[int(pi)] = i
            if not np.any(dup):
                break
            first_dup_idx = int(np.nonzero(dup)[0][0])
            bin_val = int(p[first_dup_idx])
            ids = np.where(p == bin_val)[0]
            diffs = np.abs(_rel_diff(pos[ids], np.full(ids.size, domain[bin_val], dtype=np.float64), ref=tol_ref))
            keep_idx = int(ids[int(np.argmin(diffs))])
            for ii in ids:
                if ii != keep_idx:
                    p[ii] = -1
        # accumulate matches
        matched = np.where(p >= 0)[0]
        for ii in matched:
            b = int(p[ii])
            peaks_acc[b] += pos[ii]
            x_acc[b] += vals[ii]
            counts[b] += 1
    nz = counts != 0
    peaks_acc[nz] = peaks_acc[nz] / counts[nz]
    peaks_acc[~nz] = np.nan
    x_acc[nz] = x_acc[nz] / counts[nz]
    x_acc[~nz] = np.nan
    if merge:
        merged = mergepeaks(peaks_acc.copy(), counts.copy(), x_acc.copy(), tol, tol_ref, na_drop=False)
        peaks_acc = merged.peaks
        x_acc = merged.values
        counts = merged.counts
    if na_drop and np.any(np.isnan(x_acc)):
        keep = ~np.isnan(x_acc)
        peaks_acc = peaks_acc[keep]
        x_acc = x_acc[keep]
        counts = counts[keep]
        domain = domain[keep]
    return PeakBins(peaks=peaks_acc, values=x_acc, counts=counts, tolerance=float(tol), domain=domain)

def mergepeaks(
        peaks: NDArray[np.float64],
        n: NDArray[np.int64],
        x: NDArray[np.float64],
        tol: Optional[float] = None,
        tol_ref: Literal["abs", "x", "y"] = "abs",
        na_drop: bool = True,
) -> PeakBins:
    """Merge neighboring peaks

    Behavior:
    - Within the tolerance half-window, expand left/right to find a neighborhood `[i..j]`,
      compute count-weighted averages of position and intensity, write them to the middle index `p`, and clear others.
    - Gating by counts during expansion: once counts rise again after having fallen past the mode side, treat it as a new peak and stop expansion immediately.

    Parameters:
    - peaks/n/x: aggregated positions, counts, intensities
    - tol/tol_ref: tolerance and reference; when tol is None, estimate as 1% of the average peak spacing
    - na_drop: whether to drop entries that become NaN after merging

    Returns: `PeakBins` (domain left empty, used only for semantic alignment).
    Edge cases: length mismatch or unsorted inputs raise.
    Complexity: linear scan and local merging, O(B).

    Raises:
    - ValueError: When lengths of inputs mismatch or peaks are not sorted ascending.
    """
    if peaks.size != x.size:
        logger.error("mergepeaks: length of 'peaks' and 'x' must match")
        raise ValueError("length of 'peaks' and 'x' must match")
    if np.any(np.diff(peaks[np.logical_not(np.isnan(peaks))]) < 0):
        logger.error("mergepeaks: 'peaks' must be sorted")
        raise ValueError("'peaks' must be sorted")
    if tol is None:
        ref = "abs" if tol_ref == "abs" else "y"
        diffv = _rel_diff(peaks[np.logical_not(np.isnan(peaks))], ref=ref)
        tol = 0.01 * float(np.nanmean(diffv)) if diffv.size else np.nan
    k = int(np.nanmin(np.where(np.logical_not(np.isnan(peaks)), np.arange(peaks.size), np.nan))) if np.any(
        np.logical_not(np.isnan(peaks))) else peaks.size + 1
    while k <= peaks.size - 1:
        i = k
        j = k
        left_of_mode = False
        # left expansion: stop when entering a new peak (counts rise after having fallen)
        while (i - 1) >= 0 and (not np.isnan(peaks[i - 1])) and (_udiff_scalar(peaks[i], peaks[i - 1], tol_ref) <= tol):
            if n[i - 1] < n[i]:
                left_of_mode = True
            if n[i - 1] > n[i] and left_of_mode:
                break
            i -= 1
        right_of_mode = False
        # right expansion: mirror of left, with same gating on counts
        while (j + 1) <= peaks.size - 1 and (not np.isnan(peaks[j + 1])) and (
                _udiff_scalar(peaks[j + 1], peaks[j], tol_ref) <= tol):
            if n[j + 1] < n[j]:
                right_of_mode = True
            if n[j + 1] > n[j] and right_of_mode:
                break
            j += 1
        ij = np.arange(i, j + 1)
        weights = n[ij].astype(np.float64)
        mpeaks = float(np.sum(weights * peaks[ij]) / np.sum(weights))
        mx = float(np.sum(weights * x[ij]) / np.sum(weights))
        p = int(np.floor(np.mean(ij)))
        q = ij[ij != p]
        peaks[p] = mpeaks
        x[p] = mx
        n[p] = int(np.sum(n[ij]))
        peaks[q] = np.nan
        x[q] = np.nan
        n[q] = 0
        k = j + 1
        while k <= peaks.size - 1 and np.isnan(peaks[k]):
            k += 1
    if na_drop and np.any(np.isnan(x)):
        keep = ~np.isnan(x)
        return PeakBins(peaks=peaks[keep], values=x[keep], counts=n[keep], tolerance=float(tol),
                        domain=np.array([], dtype=np.float64))
    return PeakBins(peaks=peaks, values=x, counts=n, tolerance=float(tol), domain=np.array([], dtype=np.float64))

# Single-spectrum alignment
# single spectrum alignment function (extracted from `align_sparse` core logic)
def align_single_spectrum(
        idx: NDArray[np.float64],
        dat: NDArray[np.float64],
        ref: NDArray[np.float64],
        tol: float,
        tol_ref: Literal["x", "y", "abs"]
) -> NDArray[np.float64]:
    """
    Align one spectrum to the shared reference axis.

    Parameters:
        idx: 1D peak index array (m/z), arbitrary order and may contain NaNs.
        dat: 1D intensity array aligned to `idx`.
        ref: Reference axis (ascending, NaNs filtered by caller).
        tol: Half-window tolerance used for gating.
        tol_ref: Difference reference; "x" uses `ref` as denominator (ppm),
                 "abs" uses absolute Da, "y" uses `idx` when needed.

    Returns:
        1D aligned intensity vector of length `len(ref)`.

    Raises:
        ValueError: When index/intensity lengths mismatch.
    """
    nrow = ref.size
    # Preallocate output for a single spectrum (avoid building full matrix)
    aligned = np.zeros(nrow, dtype=np.float64)
    
    if idx.size != dat.size:
        raise ValueError("lengths of index and data must match")
    
    # Filter non-finite values and keep alignment
    valid = np.isfinite(idx)
    idx = idx[valid]
    dat = dat[valid]
    
    if idx.size == 0:
        return aligned
        
    ord = np.argsort(idx)
    idx_sorted = idx[ord]
    dat_sorted = dat[ord]
    
    # Branching: Downsample vs Upsample (consistent with original `align_sparse`)
    if (nrow <= 2 * idx_sorted.size) or (np.any(np.diff(ref) < 0)):
        # Downsample: reference not much larger than spectrum
        pos = _binary_search(ref, idx_sorted, tol, tol_ref, nomatch=-1, nearest=False)
        hits = pos >= 0
        if np.any(hits):
            aligned[np.where(hits)[0]] = dat_sorted[pos[hits]]
    else:
        # Upsample: reference significantly denser than spectrum
        if tol_ref == "x": new_ref = "y"
        elif tol_ref == "y": new_ref = "x"
        else: new_ref = tol_ref
        
        start_pos = _binary_search(idx_sorted, ref, tol, new_ref, nomatch=-1, nearest=False)
        processed = np.zeros(nrow, dtype=bool)
        
        for j in range(idx_sorted.size):
            pj = int(start_pos[j])
            if pj < 0: continue
            
            # Expand to the right from `pj`
            for i in range(pj, nrow):
                if processed[i]: break
                if _udiff_scalar(ref[i], idx_sorted[j], tol_ref) > tol: break
                aligned[i] = dat_sorted[j]
                processed[i] = True
            
            # Expand to the left from `pj - 1`
            for i in range(pj - 1, -1, -1):
                if processed[i]: break
                if _udiff_scalar(ref[i], idx_sorted[j], tol_ref) > tol: break
                aligned[i] = dat_sorted[j]
                processed[i] = True
                
    return aligned

def _normalize_units(units: Literal["ppm", "mz", "relative", "absolute"]) -> Literal["relative", "absolute"]:
    """
    Unit normalization:
    - 'ppm'/'relative' map to relative units; 'mz'/'absolute' map to absolute units.
    Purpose: unify subsequent difference/tolerance reference (`tol_ref='x'|'abs'`) to avoid scattered branches.
    """
    if units in ("ppm", "relative"):
        return "relative"
    return "absolute"

def _estimate_tolerance(
        indexbins: NDArray[np.float64],
        binratio: int,
        units: Literal["relative", "absolute"],
) -> float:
    """
    Auto-estimate tolerance (half-window):
    - Use `_es_resolution(indexbins, ref=('x' if relative else 'abs'))` to estimate resolution
    - Multiply by `binratio` and round by units (relative: 6-digit half-step; absolute: 4 digits)
    Returns: float tolerance half-window used for gating during _binary_search/expansion.
    """
    ref = "x" if units == "relative" else "abs"
    tol = binratio * _es_resolution(indexbins, ref=ref)
    if units == "relative":
        return round(2.0 * tol, 6) * 0.5
    return round(tol, 4)

# [Refactor] `peak_align` now accepts `MS` and returns `AlignResult` with `MS` object
def peak_align(
        ms_data: MS,
        ref: Optional[NDArray[np.float64]] = None,
        binfun: Literal["median", "min", "max", "mean"] = "min",
        binratio: int = 2,
        tolerance: Optional[float] = None,
        units: Literal["ppm", "mz", "relative", "absolute"] = "ppm",
        output_path: str = "aligned_data.h5",
) -> AlignResult:
    """
    Main peak alignment routine (memory optimized, streaming construction).

    Behavior:
    1) Extract peak lists for domain estimation.
    2) Build reference axis (ref) and resolve tolerance.
    3) Stream through `ms_data`, align each spectrum, and build a new `MS`.

    Returns:
        AlignResult: Contains `ms_aligned`, `ref`, `count`, `freq`, `tolerance`, `units`, `binfun`, `binratio`.

    Raises:
        ValueError: If spectrum index/data lengths mismatch or reference construction fails.
    """
    logger.info(f"peak_align: start, spectra={len(ms_data)}, units={units}...")
    
    # 1. Extract Lists (only when needed)
    # Collect `index_list` when building reference peaks or estimating tolerance;
    # collect `data_list` ONLY when building reference peaks (binpeaks needs intensities).
    index_list = []
    data_list = []
    need_ref = ref is None
    need_tol_estimate = tolerance is None
    if need_ref or need_tol_estimate:
        for s in ms_data:
            index_list.append(s.mz_list)
            if need_ref:
                data_list.append(s.intensity)

    units_int = _normalize_units(units)
    tol_ref = "x" if units_int == "relative" else "abs"
    
    # 2. Build Reference Axis (Ref Logic)
    indexbins = None
    min_mz, max_mz = 0.0, 1.0
    
    if need_ref or need_tol_estimate:
        logger.info("peak_align: start domain estimation")
        indexbins, by = estimate_domain(index_list, width=binfun, units=units_int)
        logger.info(f"peak_align: domain estimated, bins={indexbins.size}, by={by}")
        if indexbins.size > 0:
            min_mz, max_mz = indexbins[0], indexbins[-1]

    # Resolve Tolerance
    if need_tol_estimate:
        tol = _estimate_tolerance(indexbins, binratio, units_int)
        logger.info(f"peak_align: tolerance resolved to {tol} (units={units_int})")
    else:
        if units in ("ppm",):
            tol = 1e-6 * float(tolerance)
        elif units in ("mz", "absolute"):
            tol = float(tolerance)
        else:
            tol = float(tolerance)
        logger.info(f"peak_align: tolerance resolved to {tol} (units={units_int})")

    counts = None
    # Build Ref
    if need_ref:
        if indexbins is not None and indexbins.size > 0:
             bins = indexbins
             logger.info(f"peak_align: using estimated domain resolution")
        else:
             forced_res = tol / binratio
             if units_int == "relative":
                 bins = _rel_sequence(min_mz, max_mz, forced_res)
             else:
                 bins = np.arange(min_mz, max_mz + forced_res, forced_res)
             logger.info(f"peak_align: forced bins built, res={forced_res} in {units_int}, range=({min_mz:.6f}, {max_mz:.6f})")

        # Binning & Merging
        peaks = binpeaks(index_list, domain=bins, xlist=data_list, tol=tol, tol_ref=tol_ref, merge=False, na_drop=False)
        merged = mergepeaks(peaks.peaks.copy(), peaks.counts.copy(), peaks.values.copy(), tol, tol_ref, na_drop=False)
        ref_vec = merged.peaks.copy()
        counts = merged.counts.copy()
    else:
        ref_vec = ref
        counts = None
        
    # Filter valid ref
    mask = np.isfinite(ref_vec)
    ref_vec = ref_vec[mask]
    if counts is not None:
        counts = counts[mask]
    else:
        counts = np.zeros(ref_vec.size, dtype=np.int64)
    if ref_vec.size and np.any(np.diff(ref_vec) < 0):
        ord = np.argsort(ref_vec)
        ref_vec = ref_vec[ord]
        if counts is not None:
            counts = counts[ord]
    
    del index_list
    del data_list
    if 'peaks' in locals(): del peaks
    if 'merged' in locals(): del merged
    gc.collect()

    # 3. Stream Alignment and MS Construction (The Optimization)
    logger.info(f"peak_align: ref built ({len(ref_vec)} peaks). Starting stream alignment...")
    
    ms_aligned = MS()
    if ms_data.meta:
        new_meta = copy.copy(ms_data.meta)
        if hasattr(new_meta, 'parser'):
            new_meta.parser = None
        if hasattr(new_meta, '_meta') and isinstance(new_meta._meta, dict):
             new_meta._meta = copy.deepcopy(ms_data.meta._meta)
        ms_aligned.meta = new_meta
        if hasattr(ms_aligned.meta, 'storage_mode'):
            ms_aligned.meta.storage_mode = 'hdf5'

    n_spectra = len(ms_data)
    n_peaks = len(ref_vec)

    with h5py.File(output_path, 'w') as f:
        f.create_dataset("mz", data=ref_vec.astype(np.float64))
        
        dset_int = f.create_dataset("intensity", 
                                   shape=(n_spectra, n_peaks), 
                                   dtype='float32', 
                                   chunks=(1, n_peaks), 
                                   compression="gzip",
                                   compression_opts=4)

        for i, spectrum in enumerate(ms_data):
            aligned_intensity = align_single_spectrum(
                spectrum.mz_list, 
                spectrum.intensity, 
                ref_vec, 
                tol, 
                tol_ref
            )

            dset_int[i, :] = aligned_intensity.astype(np.float32)

            lazy_spec = SpectrumHDF5Lazy(
                    filepath=output_path,
                    index=i,
                    coordinates=spectrum.coordinates,
                    mz_ref=ref_vec 
                )
                
            ms_aligned.add_spectrum(lazy_spec)
            
            del aligned_intensity
            
            if (i + 1) % 5000 == 0:
                logger.info(f"Aligned {i + 1}/{n_spectra} spectra")

    # Final stats
    freq = counts.astype(np.float64) / float(len(ms_data))
    
    logger.info("peak_align: alignment complete.")
    
    return AlignResult(
        ms_aligned=ms_aligned, 
        ref=ref_vec, 
        count=counts, 
        freq=freq, 
        tolerance=tol, 
        units=units_int, 
        binfun=binfun, 
        binratio=binratio
    )

def _smooth_loess_tricubic(y, window=10):
    """
    Specific Loess smoothing implementation (fixed weighting 'tri-cubic').
    Removes unused parameter `x` and fixes boundary alignment logic.
    """
    leny = len(y)
    halfw = int(np.floor(window / 2.))
    window = int(2 * halfw + 1)
    
    # --- 1. Core filter kernel construction ---
    # Local coordinates: [-halfw, ..., 0, ..., halfw]
    x1 = np.arange(1. - halfw, (halfw - 1.) + 1)
    # Tri-cubic weighting
    weight = (1. - np.abs(x1 / halfw) ** 3.) ** 1.5
    
    # Loess design matrix (Order=2)
    # V = [weight, weight*x, weight*x^2]
    V = np.vstack((weight, weight * x1, weight * x1 * x1)).T
    
    # Compute alpha via QR decomposition
    [Q, _] = linalg.qr(V, mode='economic')
    alpha = np.dot(Q[halfw - 1,], Q.T)
    
    # Apply core smoothing via filter
    yhat = signal.lfilter(alpha * weight, 1, y)
    # Correct phase shift (lfilter is causal; shift to center-align)
    yhat[halfw:-halfw] = yhat[window - 2:-1]

    # --- 2. Boundary handling ---
    # Local coordinates at boundary: [1, 2, ..., window-1]
    x1_boundary = np.arange(1., window)
    # Boundary design matrix
    V_boundary = np.vstack((np.ones(window - 1), x1_boundary, x1_boundary * x1_boundary)).T
    
    # Precompute boundary weight shape
    W_ones = np.ones((3, 1))

    # Process both ends (j from 1 to halfw)
    for j in range(1, halfw + 1):
        # Compute weights based on distance j
        # Distance vector: abs([j-1, j-2, ..., 0, ..., window-1-j])
        dist = np.abs(np.arange(1, window) - j)
        w_j = (1. - np.divide(dist, window - j) ** 3.) ** 1.5
        
        # Weighted QR decomposition
        W = (np.kron(W_ones, w_j)).T
        [Q, _] = linalg.qr(V_boundary * W, mode='economic')
        
        # Compute projection coefficients
        alpha = np.dot(Q[j - 1,], Q.T) * w_j
        
        # Left boundary: use the first window-1 points
        yhat[j - 1] = np.dot(alpha, y[:window - 1])
        
        # Right boundary: take the last window-1 points and reverse to match left-bound distance
        # Equivalent to: y[-(window-1):][::-1]
        y_right_edge = y[-(window - 1):][::-1]
        yhat[-j] = np.dot(alpha, y_right_edge)

    return yhat

def median_threshold(X):
    """
    Compute a median-based intensity threshold.

    Args:
        X: intensity data vector or matrix
    """
    md = np.median(X);
    MAD = np.median(np.abs(np.subtract(X, md))) * 1.4826;
    tval = md + 5 * MAD;
    return tval

def _find_peaks(sp, gap=3, int_thr=None):
    """
    Detect local maxima in a signal.
    gap: minimal spacing between peaks (in number of points)
    """
    gap = int(gap)
    ndp = len(sp)

    x = np.zeros(ndp + 2 * gap)
    x[:gap] = sp[0] - 1.e-6
    x[-gap:] = sp[-1] - 1.e-6
    x[gap:gap + ndp] = sp

    peak_candidate = np.ones(ndp, dtype=bool)

    for s in range(gap):
        start = gap - s - 1
        h_s = x[start:start + ndp]
        central = gap
        h_c = x[central:central + ndp]
        end = gap + s + 1
        h_e = x[end:end + ndp]

        peak_candidate &= (h_c > h_s) & (h_c > h_e)

    peakindcs = np.argwhere(peak_candidate).flatten()

    if int_thr is not None:
        peakindcs = peakindcs[sp[peakindcs] > int_thr]

    return peakindcs

def get_reference(mz, mzres, mzmaxshift, mzunits):
    """
    Compute a reference m/z axis (consensus peak locations) via kernel density over a histogram, without relying on an external grid.

    Parameters:
        mz: concatenated m/z array (1D)
        mzres: histogram resolution (bin width), unit dictated by `mzunits`
        mzmaxshift: minimal spacing for peak detection (used as `gap` in `_find_peaks`), same unit as `mzunits`
        mzunits: 'Da' or 'ppm', determines the processing space (operate in ppm space and convert back to m/z)

    Returns:
        refmz: reference m/z axis (ascending list of peak positions)

    Notes:
        1) Histogram construction → Loess/Tricubic smoothing → PCHIP interpolation → adaptive threshold → tie-breaking → peak detection → coordinate back-transform
        2) In 'ppm' mode, process in ppm space to improve numerical stability, and finally convert back to m/z units
    """
    logger.info(
        f"get_reference: start, n={0 if mz is None else np.size(mz)}, mzunits={mzunits}, mzres={mzres}, mzmaxshift={mzmaxshift}")
    if mz is None or np.size(mz) == 0:
        logger.error("get_reference: 'mz' must be a non-empty array")
        raise ValueError("'mz' must be a non-empty array")
    if mzres is None or float(mzres) <= 0:
        logger.error("get_reference: 'mzres' must be > 0")
        raise ValueError("'mzres' must be > 0")
    if mzmaxshift is None or float(mzmaxshift) <= 0:
        logger.error("get_reference: 'mzmaxshift' must be > 0")
        raise ValueError("'mzmaxshift' must be > 0")
    # unit conversion
    if mzunits == 'Da':
        rconst = 0.5
    elif mzunits == 'ppm':
        rconst = 1
        mz = np.log(mz / 1.00794) * 1.0e6
    else:
        rconst = 1

    # determine histogram range (expand by 100×bin as buffer)
    mzmin = np.min(mz) - 100 * mzres
    mzmax = np.max(mz) + 100 * mzres
    nbins = int(np.round((mzmax - mzmin) / mzres) + 1)
    logger.info(f"get_reference: histogram range=({mzmin:.6f}, {mzmax:.6f}), nbins={nbins}")

    # build histogram
    idx = np.round(
        np.divide((nbins - 1) * (mz - mzmin), (mzmax - mzmin), dtype=float) + rconst
    ).astype(int)

    histvals = np.bincount(idx, minlength=nbins)
    nbins = len(histvals)
    histidx = np.arange(1, nbins + 1)
    histmz = np.divide(
        (mzmax - mzmin) * (histidx - rconst),
        nbins - 1,
        dtype=float
    ) + mzmin

    # smoothing with tri-cubic weighting to reduce oscillations before interpolation
    histvals = _smooth_loess_tricubic(histvals)

    # super-resolution interpolation (precision boost; Da uses ×1000, ppm uses its own resolution)
    if mzunits == 'Da':
        iconst = np.max([mzres * 1000, 1])
    else:
        iconst = mzres

    histintidx = np.arange(1, nbins * iconst + 1, 1)
    logger.info(f"get_reference: iconst={iconst}, histint_len={histintidx.size}")
    inthistmz = np.divide(
        (mzmax - mzmin) * (histintidx - rconst),
        len(histintidx) - 1,
        dtype=float
    ) + mzmin

    histintvals = interpolate.pchip_interpolate(histmz, histvals, inthistmz)
    histintvals[histintvals < 0] = 0

    # adaptive threshold (median + 5×MAD) to suppress background/noise
    thrval = median_threshold(histintvals)
    histintvals[histintvals < thrval] = 0

    # break ties among equal peaks (small perturbations) to ensure unique local maxima
    histintvals[histintvals > thrval] = (
            1000 * histintvals[histintvals > thrval] +
            1e-6 * np.cumsum(np.random.uniform(size=np.sum(histintvals > 0)))
    )

    # peak detection (minimum spacing = ceil(mzmaxshift))
    maxidx = _find_peaks(histintvals, gap=np.ceil(mzmaxshift))

    # compute precise centroid positions (invert interpolated coordinates to original space)
    refmz = np.divide(
        (mzmax - mzmin) * (maxidx - rconst),
        len(histintidx) - 1,
        dtype=float
    ) + mzmin

    refmz = refmz.flatten()

    # convert back to m/z units (in ppm mode)
    if mzunits == 'ppm':
        refmz = np.exp(refmz * 1.0e-6) * 1.00794
    logger.info(f"get_reference: done, peaks={refmz.size}, units={mzunits}")
    return refmz