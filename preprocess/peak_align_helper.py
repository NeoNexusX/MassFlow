from __future__ import annotations
from logger import get_logger

logger = get_logger("peak_alignment")

from dataclasses import dataclass
from typing import List, Optional, Literal, Tuple, Dict, Iterable

import numpy as np
from numpy.typing import NDArray

from scipy import interpolate, signal, linalg


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
    - matrix: sparse aligned matrix of shape (len(ref), nspec) (nearest-value fill)
    - ref: reference peak vector (may come from binpeaks + mergepeaks)
    - count: hit count per reference peak
    - freq: hit frequency (count/nspec)
    - tolerance/units/binfun/binratio: parameter echo for auditing and reproducibility
    """
    matrix: NDArray[np.float64]
    ref: NDArray[np.float64]
    count: NDArray[np.int64]
    freq: NDArray[np.float64]
    tolerance: float
    units: Literal["relative", "absolute"]
    binfun: Literal["median", "min", "max", "mean"]
    binratio: int


def reldiff(x: NDArray[np.float64], y: Optional[NDArray[np.float64]] = None, ref: Literal["x", "y", "abs"] = "y") -> \
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
    # 统一为 float64，兼容 NumPy 2.0（移除了 np.float_）
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


def bsearch(
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
    """
    if np.any(np.diff(table) < 0):
        logger.error("bsearch: 'table' must be sorted")
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


def seq_rel(from_: float, to: float, by: float) -> NDArray[np.float64]:
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


def estres(x: NDArray[np.float64], tol: Optional[float] = None, ref: Optional[Literal["x", "abs"]] = None) -> float:
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
    - **kwargs: compatibility placeholders (e.g., verbose, chunkopts, BPPARAM), ignored here

    Returns:
    - (domain, by): reference domain and resolution step

    Edge cases and numbers:
    - When no peaks are available, return NaN step and empty domain; minimum resolution follows Cardinal (relative 5e-7, absolute 1e-4).
    Complexity: O(S·P): per-spectrum O(P) estimation plus a single aggregation.
    """
    ref = "x" if units == "relative" else "abs"
    stats = []
    for x in xlist:
        x = x[np.logical_not(np.isnan(x))]
        res = estres(x, ref=ref)
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
        seq = seq_rel(from_, to, by)
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
    1) Initial mapping: per spectrum, use `bsearch(nearest=False)` to map peak indices into `domain` bins.
    2) De-duplicate conflicts: explicitly detect duplicate hits to the same bin within a spectrum (`while(any(dup))`),
       compute relative/absolute differences to the bin center via `reldiff`, keep the closest and mark others as no-hit.
    Finally, mean-aggregate hits to obtain positions and intensities, and record counts.

    Parameters:
    - peaklist: list of peak positions per spectrum
    - domain: reference domain (None → auto-estimate range and step)
    - xlist: list of intensities per spectrum (defaults to positions)
    - tol/tol_ref: tolerance half-window and reference; when tol is None, estimate (median of 0.5×within-spectrum minimal differences)
    - merge: whether to call `mergepeaks` to merge neighboring bins after aggregation
    - na_drop: whether to drop bins with NaN values
    - **kwargs: compatibility placeholders (ignored)

    Returns: `PeakBins` containing aggregated positions, intensities, counts, and domain.

    Edge cases:
    - `peaklist[i]` and `xlist[i]` must have equal length; otherwise raise.
    - `domain` must be ascending; otherwise raise.
    Complexity: per spectrum O(P log B) initial mapping + O(P) de-duplication; overall aligns with R implementation.
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
            d = reldiff(np.array(peaks, dtype=np.float64), ref=ref)
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
                domain = seq_rel(lo, hi, tol)
    if domain.size and np.any(np.diff(domain) < 0):
        logger.error("binpeaks: 'domain' must be sorted")
        raise ValueError("'domain' must be sorted")
    peaks_acc = np.zeros(domain.size, dtype=np.float64)
    x_acc = np.zeros(domain.size, dtype=np.float64)
    counts = np.zeros(domain.size, dtype=np.int64)
    for s in range(len(peaklist)):
        p = bsearch(np.array(peaklist[s], dtype=np.float64), domain, tol, tol_ref, nomatch=-1, nearest=False)
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
            diffs = np.abs(reldiff(pos[ids], np.full(ids.size, domain[bin_val], dtype=np.float64), ref=tol_ref))
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
    - **kwargs: compatibility placeholders (ignored)

    Returns: `PeakBins` (domain left empty, used only for semantic alignment).
    Edge cases: length mismatch or unsorted inputs raise.
    Complexity: linear scan and local merging, O(B).
    """
    if peaks.size != x.size:
        logger.error("mergepeaks: length of 'peaks' and 'x' must match")
        raise ValueError("length of 'peaks' and 'x' must match")
    if np.any(np.diff(peaks[np.logical_not(np.isnan(peaks))]) < 0):
        logger.error("mergepeaks: 'peaks' must be sorted")
        raise ValueError("'peaks' must be sorted")
    if tol is None:
        ref = "abs" if tol_ref == "abs" else "y"
        diffv = reldiff(peaks[np.logical_not(np.isnan(peaks))], ref=ref)
        tol = 0.01 * float(np.nanmean(diffv)) if diffv.size else np.nan
    k = int(np.nanmin(np.where(np.logical_not(np.isnan(peaks)), np.arange(peaks.size), np.nan))) if np.any(
        np.logical_not(np.isnan(peaks))) else peaks.size + 1
    while k <= peaks.size - 1:
        i = k
        j = k
        left_of_mode = False
        # left expansion: stop when entering a new peak (counts rise after having fallen)
        while (i - 1) > 0 and (not np.isnan(peaks[i - 1])) and (_udiff_scalar(peaks[i], peaks[i - 1], tol_ref) <= tol):
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


def align_sparse(
        index_list: List[NDArray[np.float64]],
        data_list: List[NDArray[np.float64]],
        ref: NDArray[np.float64],
        tol: float,
        tol_ref: Literal["x", "y", "abs"],
) -> NDArray[np.float64]:
    """Sparse matrix domain alignment

    Behavior:
    - fast path: when a spectrum's m/z grid matches the reference axis within tolerance point-by-point, copy the entire intensity column and skip search.
    - downsample: when `len(ref) <= 2*len(idx)` or the reference domain is not ascending, for each `ref[b]` use `bsearch` over sorted `idx` to find nearest hits.
    - upsample: otherwise, center on spectrum peaks, locate a start position on the reference axis and expand left/right, filling the first unprocessed bin within the tolerance window (consistent with C++ `do_approx1`).

    Parameters:
    - index_list/data_list: peak indices and intensities per spectrum (lengths must match)
    - ref: reference peak domain (ascending; NaNs are filtered at entry)
    - tol/tol_ref: tolerance half-window and reference

    Returns: aligned 2D matrix (rows align to `ref`, columns are spectra).
    Edge cases: empty input returns a zero matrix; length mismatch raises; non-finite values are removed while keeping index-intensity order.
    Complexity: downsample O(B log P); upsample O(P log B + expansion); fast path O(B) direct copy.
    """
    ref = ref[np.isfinite(ref)]
    nrow = ref.size
    ncol = len(index_list)
    if nrow == 0 or ncol == 0:
        return np.zeros((nrow, ncol), dtype=np.float64)
    mat = np.zeros((nrow, ncol), dtype=np.float64)
    for s in range(ncol):
        idx = np.array(index_list[s], dtype=np.float64)
        dat = np.array(data_list[s], dtype=np.float64)
        if idx.size != dat.size:
            logger.error("align_sparse: lengths of index and data must match (spectrum %d)" % s)
            raise ValueError("lengths of index and data must match")
        valid = np.isfinite(idx)
        idx = idx[valid]
        dat = dat[valid]
        if idx.size == 0:
            continue
        ord = np.argsort(idx)
        idx_sorted = idx[ord]
        dat_sorted = dat[ord]
        if (nrow <= 2 * idx_sorted.size) or (np.any(np.diff(ref) < 0)):
            pos = bsearch(ref, idx_sorted, tol, tol_ref, nomatch=-1, nearest=False)
            hits = pos >= 0
            if np.any(hits):
                mat[np.where(hits)[0], s] = dat_sorted[pos[hits]]
        else:
            if tol_ref == "x":
                new_ref = "y"
            elif tol_ref == "y":
                new_ref = "x"
            else:
                new_ref = tol_ref
            start_pos = bsearch(idx_sorted, ref, tol, new_ref, nomatch=-1, nearest=False)
            processed = np.zeros(nrow, dtype=bool)
            for j in range(idx_sorted.size):
                pj = int(start_pos[j])
                if pj < 0:
                    continue
                for i in range(pj, nrow):
                    if processed[i]:
                        break
                    if _udiff_scalar(ref[i], idx_sorted[j], tol_ref) > tol:
                        break
                    mat[i, s] = dat_sorted[j]
                    processed[i] = True
                for i in range(pj - 1, -1, -1):
                    if processed[i]:
                        break
                    if _udiff_scalar(ref[i], idx_sorted[j], tol_ref) > tol:
                        break
                    mat[i, s] = dat_sorted[j]
                    processed[i] = True
    return mat


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
    - Use `estres(indexbins, ref=('x' if relative else 'abs'))` to estimate resolution
    - Multiply by `binratio` and round by units (relative: 6-digit half-step; absolute: 4 digits)
    Returns: float tolerance half-window used for gating during bsearch/expansion.
    """
    ref = "x" if units == "relative" else "abs"
    tol = binratio * estres(indexbins, ref=ref)
    if units == "relative":
        return round(2.0 * tol, 6) * 0.5
    return round(tol, 4)


def peak_align(
        index_list: List[NDArray[np.float64]],
        data_list: List[NDArray[np.float64]],
        ref: Optional[NDArray[np.float64]] = None,
        spectra: str = "intensity",
        index: Optional[str] = "mz",
        binfun: Literal["median", "min", "max", "mean"] = "min",
        binratio: int = 2,
        tolerance: Optional[float] = None,
        units: Literal["ppm", "mz", "relative", "absolute"] = "ppm",
) -> AlignResult:
    """Main peak alignment routine

    Behavior:
    - Normalize units and estimate tolerance (relative uses 6-digit half-step; absolute uses 4 digits).
    - Build reference domain and peaks when no external `ref` is provided (`estimate_domain`→`binpeaks`→`mergepeaks`).
    - Filter NaNs from the reference axis before alignment to improve stability.
    - Perform sparse alignment (`align_sparse`), including fast path and two branches.
    - Compute counts and frequencies.

    Parameters:
    - index_list/data_list: peak indices and intensities per spectrum
    - ref: external reference peak vector; None means internal construction
    - spectra/index: name placeholders (kept for Cardinal API consistency)
    - binfun/binratio/tolerance/units: parameter set consistent with Cardinal
    - **kwargs: forwarded to `estimate_domain`

    Returns: `AlignResult` (matrix, reference axis, count, frequency, tolerance, units).
    Edge cases: skip reference peak construction when `ref` is provided; skip tolerance estimation when `tolerance` is provided.
    Performance: main cost lies in domain estimation and alignment; fast path significantly accelerates.
    """
    logger.info(
        f"peak_align: start, spectra={len(index_list)}, ref={'provided' if ref is not None else 'None'}, units={units}, tolerance={tolerance}, binfun={binfun}, binratio={binratio}")
    units_int = _normalize_units(units)
    tol_ref = "x" if units_int == "relative" else "abs"
    indexbins: Optional[NDArray[np.float64]] = None
    if ref is None or tolerance is None:
        indexbins, by = estimate_domain(index_list, width=binfun, units=units_int)
        logger.info(f"peak_align: domain estimated, bins={indexbins.size}, by={by}")
    if tolerance is None:
        tol = _estimate_tolerance(indexbins, binratio, units_int)
    else:
        # Cardinal 包装层规则：当用户以 ppm 提供 tol 时需乘以 1e-6；mz/absolute 直接使用
        if units in ("ppm",):
            tol = 1e-6 * float(tolerance)
        elif units in ("mz", "absolute"):
            tol = float(tolerance)
        else:
            # 用户以相对单位提供（relative），直接使用
            tol = float(tolerance)
    logger.info(f"peak_align: tolerance resolved to {tol} (units={units_int})")
    if ref is None:
        res = estres(indexbins, ref=tol_ref)
        if units_int == "relative":
            res = round(2.0 * res, 6) * 0.5
        else:
            res = round(res, 4)
        bins = indexbins
        peaks = binpeaks(index_list, domain=bins, xlist=data_list, tol=tol, tol_ref=tol_ref, merge=False, na_drop=False)
        merged = mergepeaks(peaks.peaks.copy(), peaks.counts.copy(), peaks.values.copy(), tol, tol_ref, na_drop=False)
        ref_vec = merged.peaks.copy()
        counts = merged.counts.copy()
    else:
        ref_vec = ref
        counts = None
    # remove NaN entries from ref_vec
    mask = np.isfinite(ref_vec)
    ref_vec = ref_vec[mask]
    if counts is not None:
        counts = counts[mask]
    mat = align_sparse(index_list, data_list, ref_vec, tol, tol_ref)
    # counts None indicates an external ref was provided; initialize zeros
    if counts is None:
        counts = np.zeros(ref_vec.size, dtype=np.int64)
    freq = counts.astype(np.float64) / float(len(index_list))
    logger.info(f"peak_align: done, B={ref_vec.size}, S={len(index_list)}, tol={tol}, units={units_int}")
    return AlignResult(matrix=mat, ref=ref_vec, count=counts, freq=freq, tolerance=tol, units=units_int, binfun=binfun,
                       binratio=binratio)


# reference construction helpers
def smooth1D(x, y, window=10, method='loess', weighting='tri-cubic'):
    leny = len(y)
    halfw = np.floor((window / 2.))
    window = int(2. * halfw + 1.)
    x1 = np.arange(1. - halfw, (halfw - 1.) + 1)

    if weighting == 'tri-cubic':
        weight = (1. - np.divide(np.abs(x1), halfw) ** 3.) ** 1.5
    elif weighting == 'gaussian':
        weight = np.exp(-(np.divide(x1, halfw) * 2.) ** 2.)
    elif weighting == 'linear':
        weight = 1. - np.divide(np.abs(x1), halfw)

    if method == 'loess':
        V = (np.vstack((np.hstack(weight), np.hstack(weight * x1), np.hstack(weight * x1 * x1)))).transpose()
        order = 2
    elif method == 'lowess':
        V = (np.vstack((np.hstack((weight)), np.hstack((weight * x1))))).transpose()
        order = 1
    elif method == 'average':
        V = weight.transpose()
        order = 0

        # % Do QR decomposition
    [Q, R] = linalg.qr(V, mode='economic')

    halfw = halfw.astype(int)
    alpha = np.dot(Q[halfw - 1,], Q.transpose())

    yhat = signal.lfilter(alpha * weight, 1, y)
    yhat[int(halfw + 1) - 1:-halfw] = yhat[int(window - 1) - 1:-1]

    x1 = np.arange(1., (window - 1.) + 1)
    if method == 'loess':
        V = (np.vstack((np.hstack(np.ones([1, window - 1])), np.hstack(x1), np.hstack(x1 * x1)))).transpose()
    elif method == 'lowess':
        V = (np.vstack((np.hstack(np.ones([1, window - 1])), np.hstack(x1)))).transpose()
    elif method == 'average':
        V = np.ones([window - 1, 1])

    for j in np.arange(1, (halfw) + 1):
        # % Compute weights based on deviations from the jth point,
        if weighting == 'tri-cubic':
            weight = (1. - np.divide(np.abs((np.arange(1, window) - j)), window - j) ** 3.) ** 1.5
        elif weighting == 'gaussian':
            weight = np.exp(-(np.divide(np.abs((np.arange(1, window) - j)), window - j) * 2.) ** 2.)
        elif method == 'linear':
            weight = 1. - np.divide(np.abs(np.arange(1, window) - j), window - j)

        W = (np.kron(np.ones((order + 1, 1)), weight)).transpose();
        [Q, R] = linalg.qr(V * W, mode='economic')

        alpha = np.dot(Q[j - 1,], Q.transpose())
        alpha = alpha * weight
        yhat[int(j) - 1] = np.dot(alpha, y[:int(window) - 1])
        yhat[int(-j)] = np.dot(alpha, y[np.arange(leny - 1, leny - window, -1, dtype=int)])

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


def to_mz(value):
    """Convert from ppm space back to m/z"""
    return np.exp(value * 1.0e-6) * 1.00794


def to_ppm(value):
    """Convert from m/z to ppm space"""
    return np.log(value / 1.00794) * 1.0e6


def __findpeaks(sp, gap=3, int_thr=None):
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
        mzmaxshift: minimal spacing for peak detection (used as `gap` in `__findpeaks`), same unit as `mzunits`
        mzunits: 'Da' or 'ppm', determines the processing space (operate in ppm space and convert back to m/z)

    Returns:
        refmz: reference m/z axis (ascending list of peak positions)

    Notes:
        1) Histogram construction → tri-cubic smoothing → PCHIP interpolation → adaptive threshold → tie-breaking → peak detection → coordinate back-transform
        2) In 'ppm' mode, process in ppm space and finally convert back to m/z for numerical stability
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
        mz = to_ppm(mz)
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
    histvals = smooth1D(histmz, histvals)

    # super-resolution interpolation (precision boost; Da uses ×1000, ppm uses its own resolution)
    if mzunits == 'Da':
        iconst = np.max([mzres * 1000, 1])
    else:
        iconst = mzres

    histintidx = np.arange(1, nbins * iconst + 1, 1)
    logger.debug(f"get_reference: iconst={iconst}, histint_len={histintidx.size}")
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
    maxidx = __findpeaks(histintvals, gap=np.ceil(mzmaxshift))

    # compute precise centroid positions (invert interpolated coordinates to original space)
    refmz = np.divide(
        (mzmax - mzmin) * (maxidx - rconst),
        len(histintidx) - 1,
        dtype=float
    ) + mzmin

    refmz = refmz.flatten()

    # convert back to m/z units (in ppm mode)
    if mzunits == 'ppm':
        refmz = to_mz(refmz)

    logger.info(f"get_reference: done, peaks={refmz.size}")
    return refmz