# MassFlow

This document introduces the baseline correction module in MassFlow, focusing on the spectrum-level entry `SpectrumPreprocess.baseline_correction_spectrum`, the data-manager-level entry `Preprocess.baseline_correction`, and helper functions defined in `preprocess/baseline_correction_helper.py`.

## Overview
- Input and output
  - Spectrum level: Input `massflow.module.spectrum.Spectrum` (or `SpectrumImzML`) with 1D `intensity` and `mz_list`; output is a new `SpectrumImzML` with the same `mz_list` and coordinates and baseline-corrected `intensity`.
  - Data-manager level: Input `massflow.module.ms_data_manager.MSDataManager` (typically `MSDataManagerImzML`); output is a new `MSDataManagerImzML` whose spectra have been baseline-corrected via `Preprocess.baseline_correction`, using batch processing and disk swapping to control memory usage.
- Methods
  - LocMin (Local Minima Interpolation): baseline from local extrema anchors with optional smoothing.
  - SNIP (Statistics-Sensitive Non-linear Iterative Peak-clipping): iterative clipping with adaptive early stop.
  - SNIP (Numba-accelerated variant): same algorithm as SNIP, implemented in `snip_baseline_numba` for faster processing on large spectra.
  - ASLS (Asymmetric Least Squares): robust baseline estimation with peak preservation. Due to high computational cost, it is slower and best suited for precise correction of single or small numbers of spectra.
- Baseline scaling
  - `baseline_scale` scales the estimated baseline by `(0,1]` (default `1.0`) to prevent over-subtraction.
  - The scaled baseline is subtracted from the original signal; set `baseline_scale=1.0` to keep native algorithm behavior.

### Function Relationship Diagram

```mermaid
graph LR
  A[Baseline Correction] --> B{Method Selection}
  B --> C[LocMin Method]
  B --> D[SNIP Method]
  B --> E[ASLS Method]
  C --> F[Estimate Baseline]
  D --> F
  E --> F
  F --> G[Scale Baseline]
  G --> H[Subtract & Clip]
  H --> I[Corrected Data]
```

## Core API

### Preprocess.baseline_correction (Data manager level)
```python
massflow.preprocess.dm_pre_fun.Preprocess.baseline_correction(
  data_manager: MSDataManager,
  method: str = "asls",
  smooth: str = "none",
  span: float = 0.1,
  s: float | None = 0.0,
  upper: bool = False,
  width: int = 5,
  lam: float = 1e7,
  p: float = 0.01,
  niter: int = 15,
  baseline_scale: float = 1.0,
  m: int | None = None,
  decreasing: bool = True,
  batch_size: int = 256,
  temp_dir: str = "./temp_baseline_data",
) -> MSDataManagerImzML
```
- Description: Data-manager-level entry for baseline correction. Processes all spectra in an `MSDataManager` using the same algorithms as the spectrum-level API, streaming baseline-corrected spectra to disk.
- Input: `MSDataManagerImzML` (or subclass) containing spectra to be baseline-corrected.
- Output: A new `MSDataManagerImzML` with the same `mz_list` and coordinates as the original, and baseline-corrected `intensity` values.
- Notes:
  - Internally uses `BatchPreprocess.baseline_correction_batch` to apply `SpectrumPreprocess.baseline_correction_spectrum` over batches of spectra.
  - Batch-wise processing clears in-memory spectra (`MSDataManager.clear_batch_data_memory`) and writes corrected data out to disk, allowing large datasets to be processed within limited memory.
  - Supported methods include classic Python implementations (`"locmin"`, `"snip"`, `"asls"`) and the Numba-accelerated SNIP variant (`"snip_numba"`).

Example (data manager level):

```python
from massflow.module.mass_spectrum_set import MassSpectrumSet
from massflow.module.ms_data_manager_imzml import MSDataManagerImzML
from massflow.preprocess.dm_pre_fun import Preprocess
from massflow.tools.plot import plot_spectrum

FILE_PATH = "data/example.imzML"
ms = MassSpectrumSet()
dm = MSDataManagerImzML(ms, filepath=FILE_PATH)
dm.load_full_data_from_file()

dm_corrected = Preprocess.baseline_correction(
    data_manager=dm,
    method="snip_numba",  # use Numba-accelerated SNIP for speed
    m=50,
    baseline_scale=1.0,
    batch_size=256,
)

sp_orig = dm.ms[0]
sp_corr = dm_corrected.ms[0]

plot_spectrum(
    base=sp_orig,
    target=sp_corr,
    mz_range=(400.0, 450.0),
    intensity_range=(0.0, 2.0),
    metrics_box=True,
    title_suffix="DM_SNIPNumba",
)

dm.close()
dm_corrected.close()
```

### SpectrumPreprocess.baseline_correction_spectrum
```python
massflow.preprocess.spectrum_pre_fun.SpectrumPreprocess.baseline_correction_spectrum(
  data: Spectrum | SpectrumImzML,
  method: str = "asls",
  smooth: str = "none",
  span: float = 0.1,
  s: float | None = 0.0,
  upper: bool = False,
  width: int = 5,
  lam: float = 1e7,
  p: float = 0.01,
  niter: int = 15,
  baseline_scale: float = 1.0,
  m: int | None = None,
  decreasing: bool = True,
) -> SpectrumImzML
```
- Description: Unified entry for baseline correction of a single spectrum. Dispatches to LocMin, SNIP, or ASLS and returns a new `SpectrumImzML` with baseline-corrected `intensity` and preserved `mz_list`/coordinates.
- Supported methods: `"locmin"`, `"snip"`, `"snip_numba"`, `"asls"`.
- Notes:
  - This spectrum-level API does not perform automatic memory cleanup; it simply returns a new spectrum instance. For large datasets, prefer `Preprocess.baseline_correction` at the data-manager level.
  - To obtain the baseline explicitly (for analysis or plotting), use `baseline_corrector` directly.
- Exceptions: `ValueError` (unsupported method), `TypeError` (invalid input type)

### baseline_corrector
```python
preprocess.baseline_correction_helper.baseline_corrector(
  intensity: np.ndarray,
  index: np.ndarray | None = None,
  method: str = "asls",
  smooth: str = "none",
  span: float = 0.1,
  s: float | None = 0.0,
  upper: bool = False,
  width: int = 5,
  lam: float = 1e7,
  p: float = 0.01,
  niter: int = 15,
  baseline_scale: float = 1.0,
  m: int | None = None,
  decreasing: bool = True,
) -> tuple[np.ndarray, np.ndarray]
```
- Parameters:
  - `intensity`: 1D intensity array.
  - `index`: optional 1D coordinate array; validated for alignment and optionally used for plotting.
  - `method`: `'locmin' | 'snip' | 'snip_numba' | 'asls'`.
  - `smooth`, `span`, `s`, `upper`, `width`: LocMin options (smoothing and anchor detection).
  - `lam`, `p`, `niter`: ASLS options.
  - `baseline_scale`: scale factor in `(0,1]` applied to the estimated baseline (default `1.0`).
  - `m`, `decreasing`: SNIP options (for both pure Python and Numba versions).
- Returns:
  - `(corrected, scaled_baseline)`: baseline-corrected intensity and the scaled baseline.
- Exceptions:  
  - `ValueError`: unsupported `method`.
  - Input validation errors for non-1D or empty arrays; invalid `baseline_scale`.

### locmin_baseline
```python
preprocess.baseline_correction_helper.locmin_baseline(
  intensity: np.ndarray,
  smooth: str = "none",
  span: float = 0.1,
  s: float | None = 0.0,
  upper: bool = False,
  width: int = 5
) -> np.ndarray
```
- Description: Baseline estimation by interpolation from local extrema anchors. Detect local minima (or maxima if `upper=True`) using a windowed rule, force endpoints as anchors, then linearly interpolate. Optionally smooth the baseline via Loess (`smooth='loess', span`) or spline (`smooth='spline', s`).
- Parameters:
  - `smooth`: `'none' | 'loess' | 'spline'`
  - `span`: Loess span proportion (0 < span ≤ 1), default 0.1
  - `s`: Spline smoothing target residual sum of squares; `0.0` means interpolation
  - `upper`: Use local maxima as anchors when `True` (otherwise minima)
  - `width`: Neighborhood width for local extrema detection (default 5)
- Notes:
  - `s=0.0` yields exact interpolation through all anchor points; `s>0.0` smooths the spline toward a residual sum of squares ≈ `s`.
  - Endpoints are always included as anchors to ensure full coverage.
  - Minima and maxima detection are symmetric; maxima produce an upper envelope.
- Exceptions:
  - `ValueError`: invalid `smooth` value; `span` not in `(0,1]`; `s < 0`; `width < 3`.
  - `TypeError`: `width` is not an integer, intensity is not a 1D array.

Example:

```python
import numpy as np
from massflow.module.spectrum import Spectrum
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.tools.plot import plot_spectrum

sp = Spectrum(mz_list=mz_data, intensity=intensity_original, coordinate=[0, 0, 0])

corrected_sp = SpectrumPreprocess.baseline_correction_spectrum(
    data=sp,
    method="locmin",
    upper=False,
    width=11,
    smooth="none",
    baseline_scale=1.0,
)

plot_spectrum(
    base=sp,
    target=corrected_sp,
    mz_range=(400, 450),
    intensity_range=(0.0, 2.0),
    metrics_box=True,
    title_suffix="LocMin",
    overlay=True,
)
```
![Locmin example before/after (image to be added)]()

### snip_baseline
```python
preprocess.baseline_correction_helper.snip_baseline(
  intensity: np.ndarray,
  m: int | None = None,
  decreasing: bool = True,
) -> np.ndarray
```
- Description: Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP) baseline estimation algorithm (pure Python implementation). Performs multi-scale iterative processing to effectively estimate and remove baseline while preserving real mass spectral peaks.
- Parameters:
  - `m`: Window half-size; if None, auto-selects based on spectrum length: `min(100, max(10, n//10))`.
  - `decreasing`: Iterate from large window to small (`True`) or small to large (`False`).
- Returns:
  - `baseline`: Estimated baseline as 1D numpy array.
- Exceptions:
  - `TypeError`: `m` is not an integer.
  - `ValueError`: `m <= 0`; `intensity` is not a 1D array.

Example:

```python
from massflow.module.spectrum import Spectrum
from massflow.preprocess.spectrum_pre_fun import SpectrumPreprocess
from massflow.tools.plot import plot_spectrum

sp = Spectrum(mz_list=mz_data, intensity=intensity_original, coordinate=[0, 0, 0])

corrected_sp = SpectrumPreprocess.baseline_correction_spectrum(
    data=sp,
    method="snip",
    m=50,
    decreasing=True,
    baseline_scale=1.0,
)

plot_spectrum(
    base=sp,
    target=corrected_sp,
    mz_range=(400, 450),
    intensity_range=(0.0, 2.0),
    metrics_box=True,
    title_suffix="SNIP",
    overlay=True,
)
```
![SNIP example before/after (image to be added)]()

### snip_baseline_numba
```python
preprocess.numba.baseline_correction_numba.snip_baseline_numba(
  intensity: np.ndarray,
  m: int | None = None,
  decreasing: bool = True,
) -> np.ndarray
```
- Description: Numba-accelerated implementation of the SNIP baseline estimation algorithm. Uses parallel loops (`prange`) and JIT compilation to speed up baseline estimation on large spectra while preserving the same algorithmic behavior as `snip_baseline`.
- Parameters:
  - `intensity`: 1D input spectrum; must be non-empty.
  - `m`: Window half-size; if None, auto-selects based on spectrum length: `min(100, max(10, n//10))` (internally clamped to valid range).
  - `decreasing`: Iterate from large window to small (`True`) or small to large (`False`).
- Returns:
  - `baseline`: Estimated baseline as a float32 numpy array.
- Exceptions:
  - `ValueError`: `intensity` is not a 1D array, spectrum too short, or `m < 1` when provided.

Notes:
- When calling `baseline_corrector` with `method="snip_numba"`, this Numba implementation is used internally.


### asls_baseline
```python
preprocess.baseline_correction_helper.asls_baseline(
  intensity: np.ndarray,
  lam: float = 1e5,
  p: float = 0.001,
  niter: int = 15,
) -> np.ndarray
```
- Description: Asymmetric Least Squares (ASLS) baseline estimation algorithm. Uses weighted least squares with asymmetric penalties to estimate baseline while preserving peak signals through iterative reweighting.
- Parameters:
  - `lam`: Smoothness control parameter (positive float; larger values produce smoother baselines, typical range 1e4–1e8)
  - `p`: Asymmetry parameter (0–1; smaller values provide better peak preservation, typical range 0.001–0.1)
  - `niter`: Iteration count (positive integer; typical range 5–30)
- Returns:
  - `baseline`: Estimated baseline as 1D numpy array
 - Exceptions:
  - `ValueError`: `lam <= 0`; `p` not in `(0,1)`; `niter` not a positive integer.
  - `TypeError`:  `intensity` is not a 1D array

## Parameters and Tuning
- General
  - `baseline_scale` (`(0,1]`): smaller values reduce over-subtraction; `1.0` preserves native behavior.
- LocMin
  - `width`: Larger window reduces anchor density; start at 5–9.
    - Must be an integer and `>= 3`.
  - `upper`: Use maxima for upper envelope; keep `False` for baseline.
  - `smooth`: `'none'` for raw interpolation; `'loess'` for local smoothing; `'spline'` for global smoothing.
  - `span`: Loess span 0.05–0.3; larger values smooth more aggressively.
  - `s`: Spline smoothing target residual sum of squares; `0.0` or `None` for interpolation.
- ASLS
  - `lam` (smoothness): `1e4–1e8`, larger → smoother baseline.
  - `p` (asymmetry): `0.001–0.1`, smaller → stronger peak preservation.
  - `niter` (iterations): `5–30`, more iterations stabilize weights.
- SNIP
  - `m` (half-window): default `min(100, max(10, n//10))`; overly large windows may over-subtract.
  - `decreasing`: `True` (coarse→fine) usually more robust; `False` (fine→coarse) emphasizes local-first.

## Tips
- Ensure `mz` and `intensity` lengths match when reading NPY files.
- Use `metrics_box=True` to visualize SNR, TIC ratio and related metrics inline.
- Overlay mode: set `overlay=True` to plot original and corrected spectra on the same axis. Omit or set `overlay=False` for stacked subplots.

## References
- `massflow.preprocess.spectrum_pre_fun.py` (spectrum-level unified entry, parameter defaults)
- `massflow.preprocess.dm_pre_fun.py` (data-manager-level batch baseline correction)
- `preprocess/baseline_correction_helper.py` (LocMin, SNIP, Numba-SNIP, ASLS implementations)
- `massflow.module.spectrum.py` and `massflow.module.spectrum_imzml.py` (Spectrum and SpectrumImzML data structures)
- `massflow.tools.plot.py` (Plotting utilities for Spectrum/SpectrumImzML)

## Error Handling and Logging

- All input validation errors log via `logger.error` before raising `TypeError` or `ValueError`.
- Baseline scaling: `baseline_scale` must be a finite number in `(0,1]`; violations log and raise `ValueError`.
- LocMin:
  - `smooth` must be one of `'none' | 'loess' | 'spline'`.
  - `span` must be in `(0,1]`; `s` must be a finite number and `>= 0`.
  - `width` must be an integer and `>= 3`.
  - Loess/Spline smoothing depends on external functions (`_smooth1d`, `UnivariateSpline`); if they fail (e.g., missing dependency, invalid parameters), exceptions will propagate.
- SNIP/SNIPNumba:
  - `m` must be an integer (type violation raises `TypeError`); if provided, it must be `>= 1`.
  - `intensity` must be a non-empty 1D array; otherwise `ValueError` is raised.
- ASLS:
  - `lam` must be a positive finite number; `p` must be in `(0,1)`; `niter` must be a positive integer.
