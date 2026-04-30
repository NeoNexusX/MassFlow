# MassFlow

This document introduces the normalization module in MassFlow, focusing on:
- Pipeline-level entry point `Preprocessor.normalization`
- Core functions in `preprocess/helper/normalizer_helper.py`

## Overview

- Input and Output
  - Input is in flat array form (`mz_data`, `intensity`, `lengths`), representing m/z and intensity data for a batch of spectra.
  - Output is a normalized flat array; m/z and coordinates remain unchanged, intensity is replaced with normalized values.
- Methods
  - TIC (Total Ion Current) normalization (`tic_numba`): Scales each spectrum's intensity so the sum equals a specified value.
  - RMS (Root Mean Square) normalization (`rms_numba`): Scales each spectrum's intensity so the RMS equals a specified value.
  - Reference peak normalization (`ref_numba`): Scales each spectrum's intensity so the signal near a specified m/z is normalized to a target value.
- All methods execute per-spectrum within Numba JIT kernels, leveraging `prange` parallel loops.

### Function Relationship Diagram

```mermaid
graph LR
  A["Preprocessor.normalization"] --> B["FlatPreprocess.normalization_flat"]
  B --> C["normalizer"]
  C --> D["normalizer_numba<br/>(Numba JIT)"]
  D --> E["Normalized Data"]
```

## Core API

### Preprocessor.normalization (Pipeline Level)

```python
massflow.preprocess.api.PreprocessorAPI.normalization(
    *,
    method: str = "tic_numba",
    scale: float | None = None,
    mz_flat: np.ndarray | None = None,
    ref: float | None = None,
    ref_tolerance: float = 0.1,
) -> Self
```

- Description: Pipeline-level normalization entry point. Performs normalization on all spectra in batches via `Preprocessor` chainable calls.
- Input: Bind a data manager via `Preprocessor(dm)`.
- Output: After calling `.start()`, returns a new `MSDataManagerImzML` whose spectral intensities are normalized.
- Supported methods: `"tic_numba"`, `"rms_numba"`, `"ref_numba"`.
- Parameters:
  - `method`: Normalization method, default `"tic_numba"`.
  - `scale`: Post-normalization amplitude scaling factor (optional). If provided, intensities are multiplied by this value.
  - `mz_flat`: Flat m/z array (required for `ref_numba` method).
  - `ref`: Reference m/z value (required for `ref_numba` method).
  - `ref_tolerance`: Reference peak matching tolerance (used by `ref_numba` method), default 0.1.
- Exceptions:
  - `ValueError`: Unsupported method; invalid scale; ref_numba missing required parameters.

Example:

```python
from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess.preprocessor import Preprocessor
from massflow.tools.plot import plot_spectrum

FILE_PATH = "data/example.imzML"
dm = MSDataManagerImzML(filepath=FILE_PATH)
dm.load_head_data()

normalized = (
    Preprocessor(dm)
    .noise_reduction(method="gaussian_numba", window=7)
    .normalization(method="tic_numba")
    .start()
)

plot_spectrum(normalized.ms[0])
normalized.close()
dm.close()
```

![image-20260430140609289](https://neonexus-picture.oss-ap-southeast-1.aliyuncs.com/test/image-20260430140609289.png)

### normalizer (Helper Function Level)

```python
preprocess.helper.normalizer_helper.normalizer(
    intensity: np.ndarray,
    method: str = "tic_numba",
    scale: float | None = None,
    mz_flat: np.ndarray | None = None,
    ref: float | None = None,
    ref_tolerance: float = 0.1,
    lengths: np.ndarray | None = None,
) -> np.ndarray
```

- Description: Flat-array normalization dispatcher. Calls the corresponding Numba backend based on `method`.
- Parameters:
  - `intensity`: One-dimensional flat intensity array.
  - `lengths`: Spectrum length array, used to split processing per-spectrum within Numba.
  - Other parameters are consistent with the pipeline API.
- Returns: Normalized intensity array.
- Supported methods are limited to `"tic_numba"`, `"rms_numba"`, and `"ref_numba"`.

### FlatPreprocess.normalization_flat (Flat Batch Level)

```python
massflow.preprocess.flat_pre_fun.FlatPreprocess.normalization_flat(
    mz_data: np.ndarray,
    intensity: np.ndarray,
    lengths: np.ndarray,
    method: str = "tic_numba",
    scale: float | None = None,
    mz_flat: np.ndarray | None = None,
    ref: float | None = None,
    ref_tolerance: float = 0.1,
) -> FlatBatchResult
```

- Description: The underlying implementation of `Preprocessor.normalization`. Returns a `FlatBatchResult` containing mz_data, normalized intensity, and lengths.

## Method Details

### TIC Normalization (tic_numba)

Each spectrum's intensity is divided by its total intensity sum, making sum = 1. If `scale` is provided, the result is further multiplied by scale.

Use case: When inter-pixel total intensity comparability is needed; suitable for visualization and relative quantitative analysis.

### RMS Normalization (rms_numba)

Each spectrum's intensity is divided by its root mean square (`sqrt(mean(x^2))`), making RMS = 1. If `scale` is provided, the result is further multiplied by scale.

Use case: Less sensitive to strong peaks; more stable in data with high dynamic range signals.

### Reference Peak Normalization (ref_numba)

For each spectrum, finds the peak nearest to `ref` within the tolerance window, and uses that peak's intensity as the normalization factor. If `scale` is provided, the result is further multiplied by scale.

Use case: When an internal standard m/z is known and normalization relative to that standard is desired. Requires both `mz_flat` and `ref` parameters.

## Parameter and Tuning Recommendations

- Method Selection
  - When inter-pixel total intensity comparability is needed, prefer `"tic_numba"`.
  - When robustness to extreme values is desired, choose `"rms_numba"`.
  - When an internal standard peak is known, choose `"ref_numba"`.
- `scale` (Amplitude scaling factor)
  - Additional scaling after normalization. Common use: set to 1000 or 10000 after TIC normalization for convenient numerical display.
  - Must be finite and non-negative.
- `ref` + `mz_flat` (Reference peak normalization)
  - `ref`: Target internal standard peak m/z value.
  - `mz_flat`: Flat m/z array of the same length as intensity.
  - `ref_tolerance`: Internal standard peak search tolerance (absolute difference), default 0.1 Da.

## References

- `preprocess/helper/normalizer_helper.py`: Normalization dispatcher (normalizer)
- `preprocess/numba/normalization_numba.py`: Numba-accelerated implementation (normalizer_numba)
- `preprocess/flat_pre_fun.py`: Flat preprocessing entry point (FlatPreprocess.normalization_flat)
- `preprocess/api.py`: Chainable API (PreprocessorAPI.normalization)

## Error Handling and Logging

- Input validation: intensity must be a non-empty 1D numpy array, otherwise `TypeError` is raised.
- scale validation: Must be finite and non-negative.
- ref_tolerance validation: Must be finite and non-negative.
- method validation: Only `"tic_numba"`, `"rms_numba"`, `"ref_numba"` are supported.
- ref_numba missing parameter validation: Raises `ValueError` if `mz_flat` or `ref` is missing.
