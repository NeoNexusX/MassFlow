# MassFlow

This document introduces the normalization module in MassFlow, focusing on the unified entry `MSIPreprocessor.normalization_spectrum` and helper functions defined in `preprocess/normalizer_helper.py`. It provides an overview, API specification, example code, parameter notes, and troubleshooting tips.

## Overview
- Input and output
  - Input: `module.ms_module.SpectrumBaseModule` (with 1D `intensity` and optional `mz_list`)
  - Output: A new `SpectrumBaseModule` with the same `mz_list` and coordinates; `intensity` is replaced by normalized (and optionally scaled) values.
- Methods
  - TIC (Total Ion Current) normalization: scales intensities such that the sum equals 1.
  - RMS (Root Mean Square) normalization: scales intensities such that the RMS equals 1.
  - Median normalization: scales intensities such that the median equals 1.

### Function Relationship Diagram

```mermaid
graph LR

  A[MSIPreprocessor.normalization_spectrum] --> B{Primary Method}
  B --> C[TIC]
  B --> D[RMS]
  B --> E[Median]
  C --> F[Apply Scaling 'none'/'unit' ]
  D --> F
  E --> F
  F --> G[Return Spectrum with normalized intensity]
```

## Core API

### MSIPreprocessor.normalization_spectrum

```python
preprocess.ms_preprocess.MSIPreprocessor.normalization_spectrum(
  data: SpectrumBaseModule | SpectrumImzML,
  scale_method: str = "none",
  method: str = "tic",
  scale: float = 1.0
) -> SpectrumBaseModule | SpectrumImzML
```
- Description: Unified entry for spectrum normalization. Dispatches to TIC, RMS, or Median normalization and optionally applies post-scaling. Returns a spectrum object preserving `mz_list` and coordinates, with normalized `intensity`.
- Methods:
  - 'tic' (sum equals 1)
  - 'rms' (RMS equals 1)
  - 'median' (median equals 1)
- Parameters:
  - `scale`: Amplitude scaling factor applied after normalization (default 1.0).


### normalizer
```python
preprocess.normalizer_helper.normalizer(
  intensity: np.ndarray,
  scale_method: str = "none",
  method: str = "tic",
  scale: float = 1.0
) -> np.ndarray
```
- Description: Unified normalization dispatcher. Validates input and routes to `'tic'`, `'rms'`, or `'median'` normalization, then applies optional amplitude `scale` and `'unit'` min-max scaling.
- Parameters:
  - `intensity`: 1D numpy array to normalize.
  - `scale_method`: `'none' | 'unit'` (applied after primary normalization).
  - `method`: `'tic' | 'rms' | 'median'` (primary normalization).
  - `scale`: Amplitude scaling factor applied after normalization (default 1.0).
    - Must be a finite non-negative number.
  - Notes:
    - `scale_method` is case-insensitive and supports only `'none'` and `'unit'`.
- Returns:
  - `intensity`: normalized (and optionally scaled) 1D numpy array.
- Raises:
  - `ValueError`: unsupported `method`, `scale_method`; `scale` must be finite and non-negative.
  - `TypeError`: input not a non-empty 1D array.

### tic_normalize
```python
preprocess.normalizer_helper.tic_normalize(
  intensity: np.ndarray,
  scale_method: str = "none",
  scale: float = 1.0
) -> np.ndarray
```
- Description: Normalize by total ion current (sum). Divides by sum (TIC) when TIC > 0, then applies amplitude scaling `scale`, followed by optional `'unit'` min-max scaling.
- Parameters:
  - `scale_method`: `'none' | 'unit'` (min-max scaling to [0, 1] applied after normalization)
  - `scale`: Amplitude scaling factor applied after normalization (default 1.0)
    - Must be a finite non-negative number.
- Returns:
  - `intensity`: 1D numpy array; sum equals 1 before amplitude/optional unit scaling
- Exceptions: 
  - `ValueError`: TIC ≤ 0; unsupported `scale_method`; 
  - `TypeError`: input not a non-empty 1D array

Example: (after Savitzky-Golay denoising)

```python
import sys
import os
import numpy as np
from module.ms_module import MS
from module.ms_data_manager_imzml import MSDataManagerImzML
from preprocess.ms_preprocess import MSIPreprocessor
from tools.plot import plot_spectrum

FILE_PATH = "data\\neg-gz4.imzML"
ms = MS()
ms_md = MSDataManagerImzML(ms, filepath=FILE_PATH, coordinates_zero_based=False)
ms_md.load_full_data_from_file()
sp = ms[0]
# Denoise then normalize
denoised = MSIPreprocessor.noise_reduction_spectrum(
    data=sp,
    method="savgol",
    window=11,
    polyorder=3
)
normalized_tic = MSIPreprocessor.normalization_spectrum(
    data=denoised,
    method="tic",
    scale_method="none"
)

tic_origin = float(np.sum(denoised.intensity))
tic_after = float(np.sum(normalized_tic.intensity))
print(f"TIC normalized sum={tic_after:.6f}")

plot_spectrum(
    base = denoised,
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5),
    title_suffix='Savgol_denoised'
)

plot_spectrum(
    base=normalized_tic,
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5 / tic_origin),
    title_suffix='TIC_normalized_none'
)
```
![Original Spectrum after denoised]()
![TIC example]()

### rms_normalize
```python
preprocess.normalizer_helper.rms_normalize(
  intensity: np.ndarray,
  scale_method: str = "none",
  scale: float = 1.0
) -> np.ndarray
```
 - Description: Normalize by root mean square (RMS). Divides by RMS when RMS > 0; if RMS ≤ 0 (or NaN), raises `ValueError`. Applies amplitude `scale`, then optional `'unit'` scaling.
- Parameters:
  - `scale_method`: `'none' | 'unit'` 
  - `scale`: Amplitude scaling factor applied after normalization (default 1.0)
    - Must be a finite non-negative number.
- Returns:
  - `intensity`: 1D numpy array; RMS equals 1 (before amplitude/optional unit scaling)
- Notes:
  - Implementation matches R style: `b = sqrt(mean(x^2))`; if `b > 0` then `y = scale * x / b`, otherwise an error is raised.
- Exceptions:
  - `ValueError`: RMS ≤ 0; unsupported `scale_method`.
  - `TypeError`: input not a non-empty 1D array.

Example: (after Savitzky-Golay denoising)

```python
# Denoise then normalize (same as TIC example)
denoised = MSIPreprocessor.noise_reduction_spectrum(
    data=sp,
    method="savgol",
    window=11,
    polyorder=3
)

normalized_rms = MSIPreprocessor.normalization_spectrum(
    data=denoised,
    method="rms",
    scale_method="none"
)

# RMS before and after (RMS equals 1 after normalization when input RMS > 0)
rms_origin = float(np.sqrt(np.nanmean(np.square(denoised.intensity))))
rms_after = float(np.sqrt(np.nanmean(np.square(normalized_rms.intensity))))
print(f"RMS normalized value={rms_after:.6f}")

# Plot normalized spectrum; scale the y-range using original RMS for visibility
plot_spectrum(
    base=normalized_rms,
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5 / max(rms_origin, 1e-12)),
    title_suffix='RMS_normalized_none'
)
```
![RMS example]()

### median_normalize
```python
preprocess.normalizer_helper.median_normalize(
  intensity: np.ndarray,
  scale_method: str = "none",
  scale: float = 1.0
) -> np.ndarray
```
- Description: Normalize by the median. Divides by the median when median > 0, then applies amplitude scaling `scale`, followed by optional `'unit'` min-max scaling.
- Parameters:
  - `scale_method`: `'none' | 'unit'` (min-max scaling to [0, 1] applied after normalization)
  - `scale`: Amplitude scaling factor applied after normalization (default 1.0)
    - Must be a finite non-negative number.
- Returns:
  - `intensity`: 1D numpy array; median equals 1 before amplitude/optional unit scaling
- Exceptions:
  - `ValueError`: median ≤ 0; unsupported `scale_method`.
  - `TypeError`: input not a non-empty 1D array.

Example: (after Savitzky-Golay denoising)

```python
denoised = MSIPreprocessor.noise_reduction_spectrum(
    data=sp,
    method="savgol",
    window=11,
    polyorder=3
)
med_origin = float(np.median(denoised.intensity))
normalized_med = MSIPreprocessor.normalization_spectrum(
    data=denoised,
    method="median",
    scale_method="none"
)
med_after = float(np.median(normalized_med.intensity))
print(f"Median_value_after={med_after:.6f}")

plot_spectrum(
    base=normalized_med,
    mz_range=(500.0, 510.0),
    intensity_range=(0.0, 1.5 / med_origin),
    title_suffix='Median_normalized_none'
)
```
![Median example]()

### apply_scaling
```python
preprocess.normalizer_helper.apply_scaling(
  intensity: np.ndarray,
  scale_method: str
) -> np.ndarray
```
- Description: Apply scaling transformation after primary normalization.
  - `'none'`: return original values
  - `'unit'`: min-max scale to `[0, 1]`
- Parameters:
  - `scale_method`: `'none' | 'unit'` (case-insensitive)
- Returns:
  - `intensity`: scaled 1D numpy array
- Exceptions:
  - `ValueError`: unsupported `scale_method`.
  - `TypeError`: input not a non-empty 1D array.

### Optional 0–1 scaling (unit scaling)
```python
normalized_unit = MSIPreprocessor.normalization_spectrum(
    data=denoised,
    method="tic",          # or "median"
    scale_method="unit"    # min-max scaling to [0, 1]
)
plot_spectrum(
        base = normalized_unit,
        mz_range=(500.0, 510.0),
        intensity_range=(0.0, 0.1),
        title_suffix='TIC_normalized_unit'     
    )
```
![TIC example unit scaling]()

- Use `'unit'` to scale the normalized intensity to `[0, 1]` for consistent visualization/comparison across spectra.

## Parameters and Tuning
- General
  - Choose `'tic'` for consistent total intensity across pixels; suitable for visualization, relative quantitation.
  - Choose `'median'` to reduce the influence of extreme values; robust for noisy spectra.
  - Use `'unit'` scaling for plotting or UI normalization; avoid if you need numeric properties like sum=1 or median=1 to remain interpretable.
- TIC
  - Sensitive to large peaks; consider prior denoising/baseline correction.
- Median
  - More robust against heavy-tailed distributions; preserves rank structure.

## Troubleshooting
- `ValueError: TIC value is not greater than 0`  
  Ensure the spectrum has non-zero sum after preprocessing; avoid using entirely empty or clipped spectra.
- `ValueError: Median value is not greater than 0`  
  Ensure intensities are not all zero; consider baseline correction or avoiding aggressive clipping.
 - `ValueError: RMS value is not greater than 0`  
   Check for empty/constant spectra or all-NaN values; ensure preprocessing retains signal.
- Unsupported method or scale_method  
  Check spelling: `method='tic'|'median'`, `scale_method='none'|'unit'`.

## References
- `preprocess/normalizer_helper.py` (TIC/median implementations and scaling)
- `preprocess/ms_preprocess.py` (Unified entry point and parameter dispatch)
- `module/ms_module.py` (Data structure for `SpectrumBaseModule`)
- `tools/plot.py` (Plotting utilities for `SpectrumBaseModule`)