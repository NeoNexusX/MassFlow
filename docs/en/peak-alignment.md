# MassFlow Peak Alignment

This document systematically describes the peak alignment workflow and implementation in MassFlow, focusing on reference m/z axis construction, per-spectrum peak matching and aggregation, tolerance estimation, and the unified API. The core implementation resides in:

- `preprocess/peak_align_helper.py`
- `preprocess/ms_preprocess.py`

The intent is to present implementation-centric details: function signatures, parameter meaning, data dimensions and flow, and reproducible usage examples.

## Overview

“Reference axis construction” and “peak matching and aggregation” are layered:

- Reference domain and reference peaks
  - Automatically estimate a shared reference domain and resolution step: `estimate_domain`
  - Aggregate peaks from each spectrum into the domain under tolerance and take means: `binpeaks`
  - Merge neighboring reference peaks within tolerance: `mergepeaks`
- Density method for reference peaks
  - When `ref_method='density'`, concatenate all spectra m/z and call `get_reference(mz_all, ref_mzres, ref_mzmaxshift, ref_mzunits)` to compute the reference axis
- Sparse matrix alignment
  - On the unified reference axis, choose nearest hits for each spectrum and fill intensities: `align_sparse`
  - Auto tolerance estimation and unit normalization: `_estimate_tolerance`, `_normalize_units`
- Unified entry
  - `MSIPreprocessor.peak_alignment` collects data, chooses the reference axis source (external or internal), calls `peak_align`, and constructs a new `MS` collection

### Pipeline

```mermaid
graph TD
    A["MS object (ms_data)"] --> B("Preprocess entry MSIPreprocessor.peak_alignment");

    subgraph "Reference axis source"
        B -- ref argument --> C{External ref?};
        C -- Yes --> F["Reference m/z axis (np.ndarray)"];
        C -- No --> D{ref_method?};
        D -- 'density' --> E[get_reference];
        D -- Other/None --> G["Internal build estimate_domain→binpeaks→mergepeaks"];
        E --> F;
        G --> F;
    end

    subgraph "Peak matching and aggregation"
        A --> H["peak_align_helper.align_sparse"];
        F --> H;
        H --> I["Aligned matrix (B×S)"];
        I --> J[Count/Frequency stats];
    end

    J --> K["Build new MS (unified reference axis)"];
```


## Core Functions

### MSIPreprocessor.peak_alignment

```python
def preprocess.ms_preprocess.MSIPreprocessor.peak_alignment(
    ms_data, ref=None, units='ppm', tolerance=None,
    binfun='min', binratio=2,
    ref_method=None, ref_mzres=None, ref_mzmaxshift=None, ref_mzunits=None
) -> MS
```

- Behavior
  - Collect each spectrum’s `mz_list` and `intensity`; validate equal lengths and non-emptiness
  - Reference axis source:
    - If `ref` is provided, use it directly
    - If `ref_method=='density'`, concatenate all spectra m/z and call `get_reference(mz_all, ref_mzres, ref_mzmaxshift, ref_mzunits)` to compute the reference axis
    - Otherwise, internally estimate the reference domain and generate reference peaks in the lower-level logic
  - Call `peak_align` to obtain `AlignResult` (matrix, reference axis, count, frequency, tolerance and units)
  - Build a new `MS` where all spectra share `res.ref` and intensities come from the corresponding column
- Parameters
  - `ms_data` — Unaligned `MS` collection
  - `ref` — External reference axis; `None` means auto estimation or computing via `ref_method`
  - `units` — Tolerance units, `'ppm'` or `'mz'/'absolute'`
  - `tolerance` — Half-window tolerance; `None` means auto estimation
  - `binfun` — Resolution aggregation strategy, `'median'|'min'|'max'|'mean'`
  - `binratio` — Tolerance scaling factor, default 2
  - `ref_method` — Reference axis source; `'density'` uses kernel-density histogram `get_reference`
  - `ref_mzres` — Histogram resolution; defaults map to ppm→`100.0`, Da→`0.01`
  - `ref_mzmaxshift` — Minimal peak distance for detection (as `gap` for `__findpeaks`); defaults map to ppm→`100.0`, Da→`0.05`
  - `ref_mzunits` — Units for reference axis computation; defaults derived from `units` (ppm→`'ppm'`; mz/absolute→`'Da'`)
- Returns
  - New `MS` object whose spectra share the same reference m/z axis
- Code reference
  - `preprocess/ms_preprocess.py:159`

### peak_align

```python
def preprocess.peak_align_helper.peak_align(
    index_list, data_list, ref=None,
    spectra='intensity', index='mz',
    binfun='min', binratio=2,
    tolerance=None, units='ppm'
) -> AlignResult
```

- Function
  - Unit normalization (`'ppm'/'relative'`→`'relative'`, `'mz'/'absolute'`→`'absolute'`)
  - Reference domain estimation and reference peak construction (when `ref` is absent: `estimate_domain`→`binpeaks`→`mergepeaks`)
  - Auto tolerance estimation (when `tolerance` is absent): `_estimate_tolerance`
  - Sparse matrix alignment: `align_sparse` (fast path plus two branches)
  - Output count and frequency stats
- Parameters
  - `index_list`/`data_list` — list of peak indices and intensities per spectrum
  - `ref` — external reference peak vector; `None` means internal build
  - `binfun`/`binratio`/`tolerance`/`units` — consistent with Cardinal’s interface
- Returns
  - `AlignResult`: matrix, reference axis, count, frequency, tolerance and units
- Code reference
  - `preprocess/peak_align_helper.py:602`

### AlignResult dataclass

```python
class preprocess.peak_align_helper.AlignResult(matrix, ref, count, freq, tolerance, units, binfun, binratio)
```

- Fields
  - `matrix` — aligned matrix of shape `(len(ref), nspec)` (rows on reference axis, columns are spectra)
  - `ref` — reference peak vector (ascending)
  - `count` — per-reference peak hit count
  - `freq` — hit frequency (`count/nspec`)
  - `tolerance/units/binfun/binratio` — runtime parameter echo
- Code reference
  - `preprocess/peak_align_helper.py:32`

### estimate_domain

```python
def preprocess.peak_align_helper.estimate_domain(xlist, width='median', units='relative') -> (domain, by)
```

- Function
  - Estimate resolution and range for each spectrum’s peak indices, aggregate resolution by `width`
  - Generate an ascending reference domain according to units (relative sequence or linear sequence), and apply minimum resolution thresholds
- Returns
  - `(domain, by)`: reference domain and resolution step
- Code reference
  - `preprocess/peak_align_helper.py:241`

### binpeaks

```python
def preprocess.peak_align_helper.binpeaks(peaklist, domain=None, xlist=None, tol=None, tol_ref='abs', merge=False, na_drop=True) -> PeakBins
```

- Function
  - Initial mapping: for each spectrum, use `bsearch(nearest=False)` to map peak indices to bins in `domain`
  - Duplicate resolution: when multiple peaks from the same spectrum hit the same bin, keep the one closest to the bin center by relative/absolute difference; mark the rest as no-hit
  - Mean-aggregate hits to obtain peak positions and intensities; record counts
- Returns
  - `PeakBins`: `peaks/values/counts/tolerance/domain`
- Code reference
  - `preprocess/peak_align_helper.py:298`

### mergepeaks

```python
def preprocess.peak_align_helper.mergepeaks(peaks, n, x, tol=None, tol_ref='abs', na_drop=True) -> PeakBins
```

- Function
  - Expand left/right within tolerance to find neighborhood `[i..j]`, compute count-weighted averages for position and intensity, write to the middle index, and clear others
  - Expansion gating on counts: once counts rise again after passing the mode side, treat as a new peak and stop expansion
- Code reference
  - `preprocess/peak_align_helper.py:418`

### align_sparse

```python
def preprocess.peak_align_helper.align_sparse(index_list, data_list, ref, tol, tol_ref) -> np.ndarray
```

- Function
  - Fast path: when a spectrum’s m/z grid matches the reference within tolerance point-by-point, copy the entire intensity column directly
  - Downsample: when `len(ref) <= 2*len(idx)` or the reference domain is not ascending, search nearest hits of the reference axis in sorted indices
  - Upsample: otherwise, center on each spectrum peak, locate a start position on the reference axis, extend left/right, and fill the first unprocessed bins within tolerance
- Returns
  - Aligned 2D matrix (rows aligned to `ref`, columns are spectra)
- Code reference
  - `preprocess/peak_align_helper.py:494`

### bsearch

```python
def preprocess.peak_align_helper.bsearch(x, table, tol=0.0, tol_ref='abs', nomatch=-1, nearest=False) -> np.ndarray
```

- Function
  - In a sorted `table`, for each `x[i]`, compare left and right neighbors by relative/absolute difference; return the nearer one within tolerance (tie breaks to left)
- Code reference
  - `preprocess/peak_align_helper.py:97`

### reldiff

```python
def preprocess.peak_align_helper.reldiff(x, y=None, ref='y') -> np.ndarray
```

- Function
  - Compute relative/absolute differences; when `y` is absent, compute adjacent differences of `x`
- Code reference
  - `preprocess/peak_align_helper.py:53`

### estres

```python
def preprocess.peak_align_helper.estres(x, tol=None, ref=None) -> float
```

- Function
  - Choose the reference by robustness (MAD) between adjacent absolute vs relative differences, and return the minimal difference as resolution
  - Optional tolerance checks grid consistency
- Code reference
  - `preprocess/peak_align_helper.py:191`

### get_reference (density method)

```python
def preprocess.peak_align_helper.get_reference(mz, mzres, mzmaxshift, mzunits) -> np.ndarray
```

- Function
  - Based on kernel density estimation over a histogram, identify high-frequency peaks and return an ascending consensus m/z reference axis
  - Includes smoothing (tri-cubic weighting), PCHIP interpolation, MAD thresholding, tie-breaking, and local maxima detection
- Parameters
  - `mz` — concatenated m/z array (1D)
  - `mzres` — histogram resolution (bin width)
  - `mzmaxshift` — minimal distance for peak detection (used as `gap`)
  - `mzunits` — `'Da'` or `'ppm'` (in ppm mode, processing occurs in ppm space and converts back to m/z)
- Code reference
  - `preprocess/peak_align_helper.py:793`


## Examples

### Example 1: Use internal reference axis with auto tolerance

```python
>>> from preprocess.ms_preprocess import MSIPreprocessor
>>> ms_ali = MSIPreprocessor.peak_alignment(ms_data=ms,tolerance=10,units='ppm')
output:
INFO:     25-11-17 05:02 197 ms_preprocess - peak_alignment: ref=None, units=ppm, tolerance=10, binfun=min, binratio=2, spectra=17176, ref_method=None
INFO:     25-11-17 05:02 657 peak_alignment - peak_align: start, spectra=17176, ref=None, units=ppm, tolerance=10, binfun=min, binratio=2
INFO:     25-11-17 05:02 663 peak_alignment - peak_align: domain estimated, bins=114537, by=8e-06
INFO:     25-11-17 05:02 675 peak_alignment - peak_align: tolerance resolved to 9.999999999999999e-06 (units=relative)
INFO:     25-11-17 05:46 700 peak_alignment - peak_align: done, B=53061, S=17176, tol=9.999999999999999e-06, units=relative
INFO:     25-11-17 05:46 239 ms_preprocess - peak_alignment: ref_len=53061, tolerance=9.999999999999999e-06, units=relative, spectra=17176
INFO:     25-11-17 05:46 249 ms_preprocess - ms_aligned built: B=53061, S=17176
>>> spectrum = ms_ali[0]
>>> plot_spectrum(spectrum,
            mz_range=[420,430],
            intensity_range=[0,5],
            plot_mode="stem")
output:
  
```

![py-400-406-0-5](/Users/dre/Desktop/未命名文件夹/py-400-406-0-5.png)

```R
>>> from preprocess.ms_preprocess import MSIPreprocessor
>>> ms_peakalign <- peakAlign(ms_peakpick_realize,tolerance=10,units="ppm")
output:
ms_peakalign <- peakAlign(ms_peakpick_realize,tolerance=10,units="ppm")
detected ~5616.3 peaks per spectrum
using bin ratio of 2 to create peak bins (per tolerance half-window)
using peak bins with relative resolution of 5e-06
binning peaks to create shared reference
# processing chunk 1/20 (859 items | 38.32 MB)
# processing chunk 2/20 (859 items | 39.3 MB)
# processing chunk 3/20 (859 items | 37.61 MB)
# processing chunk 4/20 (859 items | 38.93 MB)
# processing chunk 5/20 (858 items | 40.93 MB)
# processing chunk 6/20 (859 items | 39.4 MB)
# processing chunk 7/20 (859 items | 38.69 MB)
# processing chunk 8/20 (859 items | 35.85 MB)
# processing chunk 9/20 (858 items | 38.6 MB)
# processing chunk 10/20 (859 items | 39.05 MB)
# processing chunk 11/20 (859 items | 38.79 MB)
# processing chunk 12/20 (859 items | 39.55 MB)
# processing chunk 13/20 (858 items | 39.29 MB)
# processing chunk 14/20 (859 items | 36.37 MB)
# processing chunk 15/20 (859 items | 37.66 MB)
# processing chunk 16/20 (859 items | 39.32 MB)
# processing chunk 17/20 (858 items | 38.65 MB)
# processing chunk 18/20 (859 items | 39.49 MB)
# processing chunk 19/20 (859 items | 38.93 MB)
# processing chunk 20/20 (859 items | 37.93 MB)
# collecting 17176 results from 20 chunks
merging peak bins with relative centroid differences <= 1e-05
aligned to 60979 reference peaks with relative tolerance 1e-05 (10 ppm)
>>> plot(ms_peakalign,pixel=1,xlim=c(400,406),ylim=c(0,5),type="h")
output:
```

![car-400-406-0-5](/Users/dre/Desktop/未命名文件夹/car-400-406-0-5.png)

### Example 2: Density method with explicit parameters

```python
>>> from preprocess.ms_preprocess import MSIPreprocessor
>>> ms_aligned = MSIPreprocessor.peak_alignment(
...     ms_data=ms,
...     ref_method='density',      # use kernel-density histogram method
...     ref_mzres=3,               # histogram resolution
...     ref_mzmaxshift=3,          # minimal distance for peak detection
...     units='ppm',               # tolerance units
...     tolerance=None,            # auto tolerance estimation
...     binfun='min',              # resolution aggregation strategy
...     binratio=2                 # tolerance scaling factor
... )
>>> print(f"ref_len={len(ms_aligned[0])}")
output:
ref_len=19523
```


## Parameter and Unit Conventions

- Unit normalization:
  - `'ppm'/'relative'` normalized to relative units; `'mz'/'absolute'` to absolute units
  - Difference reference in alignment/search: relative→`tol_ref='x'`; absolute→`tol_ref='abs'`
- Tolerance estimation:
  - When `tolerance=None`, first estimate resolution `estres`, then multiply by `binratio`
  - Rounding rules: relative uses a 6-digit half-step (`round(2*res, 6)/2`), absolute uses 4 digits
- Reference axis:
  - When an external `ref` is provided, skip internal construction
  - `'density'` requires or auto-maps `ref_mzres/ref_mzmaxshift/ref_mzunits`


## Data Flow Summary

1. Collect all spectra’s `mz_list` and `intensity`
2. Choose the reference axis source (external `ref`, `'density'` method, or internal estimation of domain and peaks)
3. Normalize units and parse/estimate tolerance
4. Align on the reference axis to obtain a `(B, S)` matrix
5. Build a new `MS` with unified m/z axis and filled intensities


## Extensibility Notes

- External grid alignment: if you already have a high-quality external grid (e.g., from third-party tools), pass it as `ref` to skip domain construction
- Tolerance strategy: tune `binratio` according to instrument drift and resolution; or specify `tolerance` manually
- Scale optimization: for sparse peak lists, prefer absolute units and a shorter reference axis to reduce search and extension costs


## Glossary

- Reference Axis: row index of the aligned matrix; ascending m/z list
- Relative tolerance (ppm): half-window proportional to m/z
- Absolute tolerance (Da): fixed m/z half-window
- downsample/upsample: alignment strategies when the reference axis is shorter/longer
- Count-weighted mean: aggregation strategy using hit count as weights