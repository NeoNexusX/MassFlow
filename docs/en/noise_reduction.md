# MassFlow Noise Suppression and Filtering (filter.py)

This document introduces the noise suppression and filtering module in MassFlow, focusing on the functions in `preprocess/filter.py` and their coordinated use with `MSIPreprocessor.noise_reduction`. The content includes API descriptions, example code, parameters and tuning suggestions, application scenarios, and common issues.
- Related Modules
  - `preprocess/filter.py`
  - `preprocess/ms_preprocess.py` (provides a unified entry point for `MSIPreprocessor.noise_reduction`)

## Overview

- Input and Output
  - Input: `module.ms_module.SpectrumBaseModule` (1D `intensity` must exist; `mz_list` is optional)
  - Output: 1D smoothed `intensity` or a new `SpectrumBaseModule` (returned via `MSIPreprocessor.noise_reduction`)
- Algorithm Categories
  - Time-domain convolution smoothing: `smooth_signal_ma` (moving average/custom kernel), `smooth_signal_gaussian` (discrete Gaussian)
  - Polynomial fitting smoothing: `smooth_signal_savgol` (Savitzky-Golay)
  - Wavelet denoising: `smooth_signal_wavelet` (thresholding based on PyWavelets)
  - Smoothing based on m/z neighborhood search: `smooth_ns_signal_ma`, `smooth_ns_signal_gaussian`, `smooth_ns_signal_bi` (bilateral)
- Preprocessing
  - `smooth_preprocess`: Sets negative intensities to zero, cleans up data references to avoid subsequent algorithm exceptions.
  - Neighborhood search preprocessing (internal use): Fills NaN intensities, falls back to equidistant indices if `mz_list` is missing.

### Function Relationship Diagram

```mermaid
graph LR
  A[MSIPreprocessor.noise_reduction(data, method,...)] --> P[smooth_preprocess]
  A --> MA[smooth_signal_ma]
  A --> G[smooth_signal_gaussian]
  A --> SG[smooth_signal_savgol]
  A --> WT[smooth_signal_wavelet]
  A --> MA_NS[smooth_ns_signal_ma]
  A --> G_NS[smooth_ns_signal_gaussian]
  A --> BI_NS[smooth_ns_signal_bi]
```

## Core API

### MSIPreprocessor.noise_reduction

```python
preprocess.ms_preprocess.MSIPreprocessor.noise_reduction(
  data: SpectrumBaseModule | SpectrumImzML,
  method: str = "ma",
  window: int = 2,
  sd: float = 2,
  sd_intensity: float | None = None,
  p: int = 2,
  coef: np.ndarray | None = None,
  polyorder: int = 2,
  wavelet: str = "db4",
  threshold_mode: str = "soft"
) -> SpectrumBaseModule | SpectrumImzML
```

- Description: Unified entry point for denoising. Dispatches to the specific filter implementation based on `method`, and returns a spectrum object with the same coordinates as the input, where `intensity` is the smoothed result.
- Supported `method`s:
  - `"ma"`, `"gaussian"`, `"savgol"`, `"wavelet"`
  - `"ma_ns"`, `"gaussian_ns"`, `"bi_ns"`
- Returns: A new `SpectrumBaseModule` (or `SpectrumImzML`) instance, with `mz_list` and coordinates preserved, and `intensity` replaced with the smoothed result.
- Exceptions: `ValueError` (unsupported method), `TypeError` (invalid input type).

### smooth_signal_ma

```python
preprocess.filter.smooth_signal_ma(
  x: SpectrumBaseModule,
  coef: np.ndarray | None = None,
  window: int = 5
) -> np.ndarray
```

- Description: Moving average (or custom convolution kernel) smoothing. Uses `edge` padding for 1D boundaries and automatically normalizes weights.
- Parameters:
  - `coef`: Convolution kernel; if `None`, a uniform kernel of length `window` is used.
  - `window`: Window length (positive integer; if even, it is automatically adjusted to an odd number).
- Returns: A 1D `intensity` array of the same length as the input.
- Exceptions: `ValueError` (`window <= 0` when `coef` is not provided), `TypeError` (`intensity` is not 1D).

### smooth_signal_gaussian

```python
preprocess.filter.smooth_signal_gaussian(
  x: SpectrumBaseModule,
  sd: float | None = None,
  window: int = 5
) -> np.ndarray
```

- Description: Discrete Gaussian kernel smoothing; defaults to `sd = window / 4`.
- Dependencies: `scipy.stats.norm` (if unavailable, falls back to a `numpy` implementation).
- Parameters: `sd` (Gaussian standard deviation), `window` (odd length).
- Exceptions: `ValueError` (window is not positive), `TypeError` (`intensity` is not 1D).

### smooth_signal_savgol

```python
preprocess.filter.smooth_signal_savgol(
  x: SpectrumBaseModule,
  window: int = 5,
  polyorder: int = 2
) -> np.ndarray
```

- Description: Savitzky-Golay polynomial fitting smoothing; automatically ensures `window` is odd and `window >= 3`.
- Dependencies: `scipy.signal.savgol_filter`.
- Parameters: `polyorder < window` (if not satisfied, it is automatically lowered).
- Returns: A 1D `intensity` array of the same length as the input.

### smooth_signal_wavelet

```python
preprocess.filter.smooth_signal_wavelet(
  x: SpectrumBaseModule,
  wavelet: str = "db4",
  threshold_mode: str = "soft"
) -> np.ndarray
```

- Description: Wavelet threshold denoising; automatically estimates noise and processes it using the Donoho-Johnstone threshold; the reconstructed length is strictly aligned with the input (truncated/padded if necessary).
- Dependencies: `pywt` (PyWavelets).
- Parameters: `wavelet` (e.g., `db4`, `db8`, `haar`, `coif2`), `threshold_mode` (`soft`/`hard`).

### Neighborhood Search Preprocessing (Internal)

```python
preprocess.filter.smooth_ns_signal_pre(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 1
) -> tuple[np.ndarray, np.ndarray, np.ndarray]
```

- Description: Performs kNN search on the m/z axis (`cKDTree`) and returns neighbor intensities, distances, and indices.
- Behavior:
  - If `mz_list` is missing or its length does not match the intensity, it falls back to `np.arange(N)`.
  - NaN values in `intensity` are filled with the neighborhood mean (if the entire neighborhood is NaN, the global median is used; if still NaN, it is set to 0).
  - If `k > N`, `k` is automatically adjusted to `N`.
- Parameters: `k >= 1`, `p >= 1` (Minkowski distance parameter; `p=1` for Manhattan, `p=2` for Euclidean).

### smooth_ns_signal_ma

```python
preprocess.filter.smooth_ns_signal_ma(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 1
) -> np.ndarray
```

- Description: Performs equal-weight average smoothing over the m/z neighborhood (row-normalized).
- Parameters: `k` (number of neighbors), `p` (distance metric).

### smooth_ns_signal_gaussian

```python
preprocess.filter.smooth_ns_signal_gaussian(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 1,
  sd: float | None = None
) -> np.ndarray
```

- Description: Applies Gaussian weighting to neighborhood distances; the exponent is `clip`ped before numerical underflow, and row normalization avoids division by zero.
- Default: `sd = median(max_row_distance) / 2`.
- Parameters: `k`, `p`, `sd`.

### smooth_ns_signal_bi (Bilateral)

```python
preprocess.filter.smooth_ns_signal_bi(
  x: SpectrumBaseModule,
  k: int = 5,
  p: int = 2,
  sd_dist: float | None = None,
  sd_intensity: float | None = None
) -> np.ndarray
```

- Description: Bilateral weighting that considers both m/z distance and intensity difference; suitable for edge-preserving (peak position) smoothing.
- Defaults:
  - `sd_dist = median(max_row_distance) / 2`
  - `sd_intensity = scipy.stats.median_abs_deviation(intensity, scale="normal")`
- Numerical Stability: `clip` before exponential underflow, row normalization.

### smooth_preprocess

```python
preprocess.filter.smooth_preprocess(data: SpectrumBaseModule) -> SpectrumBaseModule
```

- Description: General-purpose preprocessing before smoothing, sets negative intensities to zero and cleans up references, returning the updated object.

## Quick Start Examples

> The following examples are simplified demonstrations. Image placeholders are left blank.

### 1. Moving Average (MA)
```python
import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from module.ms_module import MS, SpectrumBaseModule
from module.ms_data_manager_imzml import MSDataManagerImzML
from preprocess.ms_preprocess import MSIPreprocessor
# Direct parameter settings
window_size = 7  # Window size for moving average
output_filename = "ma_denoised_plot.png"  # Output plot filename
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem

# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing (using directly set window size)
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="ma",
    window=window_size
)
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='MA'
)
# Save results
mz_den = np.asarray(denoised.mz_list, dtype=float)
inten_den = np.asarray(denoised.intensity, dtype=float)
out_base = f"{input_file_stem}_ma_denoised"
np.save(f"{out_base}_intensity.npy", inten_den)
np.save(f"{out_base}_mz.npy", mz_den)

print(f"Processing completed! Window size: {window_size}")
print(f"Plot saved as: {output_filename}")
print(f"Data saved as: {out_base}_intensity.npy and {out_base}_mz.npy")

```
![Example before and after MA filtering (image to be added)]()
## For subsequent examples, the import and save sections are consistent with MA and will be omitted.
### 2. Gaussian Smoothing

```python
# Parameter settings
window_size = 7  # Window size for Gaussian filter
standard_deviation = 2.0  # Standard deviation for Gaussian filter
output_filename = "gaussian_denoised_plot.png"  # Output plot filename
# Check if window size is odd
if window_size % 2 == 0:
    print("Warning: Window size should be an odd number. Incrementing by 1.")
    window_size += 1
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Apply Gaussian noise reduction
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="gaussian",
    window=window_size,
    sd=standard_deviation
)
print(f"Denoising complete. Denoised spectrum length: {len(denoised.intensity)}")
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='seagreen',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Gaussian'
)
```

![Example before and after Gaussian filtering (image to be added)]()

### 3. Savitzky-Golay Smoothing (Peak-preserving)

```python
# Parameter settings
window_size = 11  # Window size for Savitzky-Golay filter (odd number, >=3)
polyorder = 3  # Polynomial order for Savitzky-Golay filter (must be less than window)
output_filename = "savgol_denoised_plot.png"  # Output plot filename
# Check window parameters
if window_size % 2 == 0:
    window_size += 1
    print(f"Warning: Window size must be odd. Adjusted to {window_size}.")
if window_size < 3:
    window_size = 3
    print(f"Warning: Window size must be at least 3. Adjusted to {window_size}.")
if polyorder >= window_size:
    polyorder = window_size - 1
    print(f"Warning: Polyorder must be less than window. Adjusted to {polyorder}.")
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Apply Savitzky-Golay noise reduction
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="savgol",
    window=window_size,
    polyorder=polyorder
)
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='darkorange',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Savitzky-Golay'
)

```

![Example before and after SG filtering (image to be added)]()

### 4. Wavelet Denoising

```python
# Parameter settings
wavelet_type = "db4"  # Wavelet type (e.g., 'db4', 'sym8', 'coif5')
threshold_mode = "hard"  # Thresholding mode ('soft' or 'hard')
output_filename = "wavelet_denoised_plot.png"  # Output plot filename
# Data path
FILE_PATH = "data\\neg-gz4.imzML"
input_file_stem = Path(FILE_PATH).stem
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]

print(f"Original spectrum loaded. Length: {len(sp.intensity)}")
print(f"Processing with Wavelet filter: wavelet={wavelet_type}, mode={threshold_mode}")
# Apply Wavelet noise reduction
denoised = MSIPreprocessor.noise_reduction(
    data=sp,
    method="wavelet",
    wavelet=wavelet_type,
    threshold_mode=threshold_mode
)
# Plotting
denoised.plot(
    save_path=output_filename,
    figsize=(12, 8),
    dpi=300,
    color='purple',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Wavelet'
)
```

![Example before and after wavelet denoising (image to be added)]()

### 5. Neighborhood Search Moving Average (ma_ns)

```python
# Parameters
K = 7  # Number of neighbors
P = 1  # Minkowski distance parameter
OUTPUT_FILENAME = "ma_ns_denoised_plot.png"
FILE_PATH = "data\\neg-gz4.imzML"
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing
denoised_sp = MSIPreprocessor.noise_reduction(
    data=sp,
    method="ma_ns",
    window=K,  # For ns methods, window is used as k
    p=P
)
# Plotting
denoised_sp.plot(
    save_path=OUTPUT_FILENAME,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='MA-NS'
)
```

![Example before and after MA_NS filtering (image to be added)]()

### 6. Neighborhood Search Gaussian (gaussian_ns)

```python
# Parameters
K = 7
P = 2
SD = None  # Let the function determine sd automatically
OUTPUT_FILENAME = "gaussian_ns_denoised_plot.png"
FILE_PATH = "data\\neg-gz4.imzML"
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing
denoised_sp = MSIPreprocessor.noise_reduction(
    data=sp,
    method="gaussian_ns",
    window=K,
    p=P,
    sd=SD
)
# Plotting
denoised_sp.plot(
    save_path=OUTPUT_FILENAME,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='Gaussian-NS'
)
```
![Example before and after Gaussian_NS filtering (image to be added)]()

### 7. Neighborhood Search Bilateral (bi_ns)

```python
# Parameters
K = 7
P = 2
SD_DIST = None
SD_INTENSITY = None
OUTPUT_FILENAME = "bi_ns_denoised_plot.png"
FILE_PATH = "data\\neg-gz4.imzML"
# Load data
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoising processing
denoised_sp = MSIPreprocessor.noise_reduction(
    data=sp,
    method="bi_ns",
    window=K,
    p=P,
    sd=SD_DIST,
    sd_intensity=SD_INTENSITY
)
# Plotting
denoised_sp.plot(
    save_path=OUTPUT_FILENAME,
    figsize=(12, 8),
    dpi=300,
    color='steelblue',
    plot_mode='line',
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    original=sp,
    metrics_box=True,
    title_suffix='BI-NS'
)
```

![Example before and after Bilateral_NS filtering (image to be added)]()

### 8. Unified Entry Point (Recommended): MSIPreprocessor.noise_reduction

```python
from preprocess.ms_preprocess import MSIPreprocessor

pre = MSIPreprocessor()
spec_ma = pre.noise_reduction(data=spec, method='ma', window=7)
spec_gauss = pre.noise_reduction(data=spec, method='gaussian', window=9, sd=2.5)
spec_sg = pre.noise_reduction(data=spec, method='savgol', window=9, polyorder=3)
spec_wt = pre.noise_reduction(data=spec, method='wavelet', wavelet='db4', threshold_mode='soft')

spec_ma_ns = pre.noise_reduction(data=spec, method='ma_ns', window=7, p=2)
spec_gauss_ns = pre.noise_reduction(data=spec, method='gaussian_ns', window=7, sd=None, p=2)
spec_bi_ns = pre.noise_reduction(data=spec, method='bi_ns', window=7, sd=None, sd_intensity=None, p=2)
```

## Parameter Descriptions and Tuning Suggestions

- `window` (Convolution/SG window)
  - Positive integer; even numbers are automatically converted to the nearest odd number (for center alignment).
  - Recommendation: 5–15 (adjust based on peak width and noise intensity).
- `coef` (Custom convolution kernel)
  - Automatically normalized when passed; its length determines the effective window.
- `sd` (Standard deviation of Gaussian kernel)
  - Classic convolution defaults to `sd = window / 4`.
  - Neighborhood Gaussian defaults to `sd = median(max_row_distance) / 2` (data-adaptive).
- `polyorder` (Polynomial order for SG)
  - Must be less than `window`; recommended 2–4.
- `wavelet` and `threshold_mode`
  - Recommendation: `db4` + `soft`; for strong noise, try `hard`.
- `k` (Number of neighbors for neighborhood search)
  - `k >= 1`; too large can lead to over-smoothing; generally 5–11.
- `p` (Minkowski distance parameter)
  - `p=1` for Manhattan, `p=2` for Euclidean; higher `p` emphasizes distance differences more.
- `sd_intensity` (Scale of intensity difference weight for bilateral filter)
  - Defaults to the `median_abs_deviation` of the intensity (robust).

## Use Cases

- Basic denoising and quick smoothing: `ma`, `gaussian`
- Peak shape preservation and quantitative analysis: `savgol` (suitable for sharp peaks where peak position should not shift).
- Strong/mixed noise: `wavelet` (multi-scale robust thresholding).
- Irregular m/z sampling or smoothing based on m/z distance: `ma_ns`, `gaussian_ns`.
- Edge-preserving (peak position) denoising: `bi_ns` (bilateral filter is more faithful around sharp peaks).

## Common Problems and Troubleshooting

- Incorrect input shape
  - Error: `TypeError: x.intensity must be a 1D numpy array`
  - Solution: Ensure `intensity` is a 1D `np.ndarray`.
- Invalid window parameter
  - Error: `ValueError: window must be a positive integer`
  - Solution: `window >= 3`, and it will be converted to an odd number.
- SG constraint not met
  - Situation: `polyorder >= window`; the library will raise an error or overfit.
  - Solution: Lower `polyorder` (the implementation already automatically corrects it to `window-1`).
- Missing dependencies
  - Dependencies: `numpy`, `scipy` (`stats`, `signal`, `spatial`), `pywt`.
  - Solution: Install according to `requirements.txt`; or if `scipy` is missing for `gaussian`, a `numpy` implementation is used.
- Neighborhood search boundaries
  - Situation: `mz_list` is missing or has a mismatched length.
  - Handling: Automatically falls back to `mz_list = np.arange(N)` and fills `NaN` intensities.
- Numerical stability
  - Exponential weight underflow: The implementation `clip`s to a safe range.
  - Row normalization: Avoids division by zero; if a row sum is 0, it is set to 1.
- Invalid method value for the unified entry point
  - Error: `Unsupported smoothing method`
  - Solution: Use one of `ma, gaussian, savgol, wavelet, ma_ns, gaussian_ns, bi_ns`.

## Practical Suggestions

- Preprocessing: First, use `smooth_preprocess` to clean up negative intensities.
- Window selection: Use the target peak width as a reference; a large window will flatten sharp peaks.
- Adaptive scales: Prefer using the default `sd`/`sd_intensity` (based on robust data estimation).
- Evaluation metrics: Compare changes in peak height/width/noise floor to avoid over-smoothing.

## References

- `preprocess/filter.py` (Core implementation of filtering/denoising)
- `preprocess/ms_preprocess.py` (Unified entry point and parameter dispatch)
- `module/ms_module.py` (Data structure and visualization for `SpectrumBaseModule`)