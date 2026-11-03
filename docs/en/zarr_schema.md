# MassFlow MSI Zarr Schema

Schema version: `1.1`. This is the MassFlow schema version, independently of the Zarr storage-format version. The current writer uses Zarr's default storage format; group attributes are accessed through `.attrs` (in Zarr v3 they are stored in `zarr.json`).

## Supported layouts

| Layout | `/ion_image` | `/spectra` | Encoding |
| --- | --- | --- | --- |
| Ion-only | Present | Absent | continuous |
| Spectra-only | Absent | Present | continuous or processed |
| Combined | Present | Present | continuous for both views |

`MSDataManagerZarr` writes spectra-only data. `IonImageDataManagerZarr` defaults to combined output (`include_spectra=True`); use `include_spectra=False` for ion-only output. The low-level `MSIZarrWriter` defaults `include_spectra` to `False`.

All stores contain `axes/coordinates` and `metadata`. Continuous stores also contain `axes/mz`. Version 1.0 and the legacy `/data` group are not supported.

```text
sample.zarr/                  # combined layout
├── metadata/                # scalar attributes; optional mask array
├── axes/
│   ├── coordinates          # (n_pixels, 3), <u4
│   └── mz                   # (n_ions,), <f4 or <f8
├── ion_image/               # attrs: row_axis="ion", encoding="continuous"
│   ├── intensity            # flat, ion-major
│   └── offsets              # (n_ions + 1,), <i8
├── spectra/                 # attrs: row_axis="pixel", encoding="continuous"
│   ├── intensity            # flat, pixel-major
│   └── offsets              # (n_pixels + 1,), <i8
└── stats/
    ├── tic                  # required for combined
    ├── mean_spectrum        # required for combined
    └── rms                  # optional
```

## Required root attributes

```json
{
  "format": "massflow.msi_zarr",
  "format_version": "1.1",
  "write_state": "complete",
  "coordinate_order": ["x", "y", "z"],
  "coordinate_base": 1,
  "spatial_shape": [93, 227]
}
```

`spatial_shape` is `[height, width]`. A writer uses `write_state="writing"` until finalization; readers require `"complete"`.

## Axes and coordinates

`axes/coordinates` has shape `(n_pixels, 3)` and dtype `<u4`, with one-based `[x, y, z]` coordinates. Valid bounds are `1 <= x <= width` and `1 <= y <= height`.

Ion-major finalization normalizes z to 1 and keeps the last occurrence of duplicate `(x, y)` positions, retaining the surviving input-row order. Both views and per-pixel statistics use the same retained rows in combined stores. Pixel-major writing preserves all input rows and their z values.

`axes/mz` is required for continuous data, with dtype `<f4` or `<f8` and `unit="m/z"`. Processed spectra instead store per-point m/z in `spectra/mz`.

## Row arrays

Every layout group contains one-dimensional `intensity` and `offsets` arrays. The writer defaults intensity to float32; readers accept floating-point intensity arrays. Offsets have dtype `<i8` and satisfy:

```text
offsets[0] = 0
offsets[i + 1] >= offsets[i]
offsets[-1] = len(intensity)
```

| View | Rows | Flat intensity length | Offsets |
| --- | --- | --- | --- |
| Ion-major | Ions | `n_ions * n_pixels` | `[0, n_pixels, 2*n_pixels, ...]` |
| Continuous spectra | Pixels | `n_pixels * n_ions` | `[0, n_ions, 2*n_ions, ...]` |
| Processed spectra | Pixels | `sum(lengths)` | Cumulative spectrum lengths |

A row is `intensity[offsets[i]:offsets[i+1]]`. For processed data, the same slice selects `spectra/mz`; its length equals the intensity length. Processed writers accept empty spectra, represented by repeated offsets. Continuous rows must match the shared axis length.

An ion row contains sampled pixel values, not a flattened dense image. To display it, scatter the row using `axes/coordinates` into a zero-initialized `(height, width)` array.

## Metadata

`metadata` attributes store scalar fields such as `name`, numeric `version`, `spectrum_count_num`, image dimensions, pixel sizes, spatial offsets, instrument model, acquisition fields, encoding flags, and spectrum-type flags. The optional `metadata/mask` array has image shape and marks occupied locations. Field names omit the `meta_` prefix used in printed output; see [Metadata](./metadata.md).

## Statistics

| Array | Shape | Definition | Attributes |
| --- | --- | --- | --- |
| `stats/mean_spectrum` | `(n_ions,)` | Mean over retained sampled pixels for each ion | `source="/ion_image"`, `reduced_axis="pixel"`, `denominator="sampled_pixels"` |
| `stats/tic` | `(n_pixels,)` | Sum of each spectrum's intensities | `source="/spectra"`, `reduced_axis="mz"` |
| `stats/rms` | `(n_pixels,)` | `sqrt(sum(intensity**2) / spectrum_length)` | `source="/spectra"`, `reduced_axis="mz"`, `denominator="spectrum_points"` |

All three statistics use float64. Empty processed spectra have TIC=RMS=0. Mean-spectrum normalization uses the number of sampled pixels, not the dense image area.

Combined writing always enables TIC and mean spectrum, with RMS optional. For ion-only writing mean spectrum is optional; spectra-only writing can enable TIC and RMS independently.

## Reader behavior

`MSIZarrReader` uses `/ion_image` as the main layout when available; a combined store additionally exposes pixel spectra through `read_spectrum()` and `read_flat_batch()`. Spectra-only stores use `/spectra` as their main layout. An ion-only store does not provide pixel-major spectrum reads.

The shared m/z axis is cached per reader and reused by single-spectrum and batch reads. Copy a returned shared axis before modifying it. Non-contiguous flat-batch selections are gathered by chunk while preserving the requested row order.

For the detailed Chinese reference, see [Zarr 格式规范](/zh/zarr_schema).
