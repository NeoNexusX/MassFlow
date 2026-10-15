# QC metrics, review triggers and actions

Triggers are prompts to look closer, not pass/fail limits: a value past a
trigger needs an explanation (from the image, the acquisition design, or the
user) before the verdict.

## File and acquisition (`load_imzml`)

| Field | Meaning | Review trigger | Action |
|---|---|---|---|
| `warnings` | Automatic repairs while opening | Any | Report; geometry taken from coordinates means the header was wrong or incomplete |
| `acquisition.polarity` + `_source` | Ion polarity and its evidence | `null`, `mixed`, or `file name` only | Ask the user if adduct steps or the reference-ion check depend on it; `set_acquisition` |
| `acquisition.analyzer_class` | tof / orbitrap / ftms / quadrupole / unknown | `unknown` | Ask if tolerances matter; defaults then assume ~10 ppm |
| `processing_history` | Steps that produced this file | Present | Evaluate as a processed output: compare with its input when available |
| `files.ibd_mb` | Binary size | Far smaller than spectra × points × bytes suggests | Possible truncated copy: re-copy before analysis |

## Spatial layout (`profile_spatial`)

| Field | Meaning | Review trigger | Action |
|---|---|---|---|
| `occupancy` | Acquired / declared grid pixels | < 0.3 | If `bounding_box_fill` is also low: tissue-shaped raster (normal). If fill is high: declared grid larger than the raster |
| `bounding_box_fill` | Acquired / bounding-box pixels | < 0.5 | Irregular acquisition region; expect no off-tissue frame |
| `duplicate_coordinates` | Spectra reusing an (x, y) | > 0 | Hold unless `z_levels` > 1 explains it; images overwrite these spectra |
| `z_levels` | Distinct z values | > 1 | Multi-section/3D: evaluate per plane |
| `pixel_size_um` | Pixel pitch | Missing | No physical scale; ask the user |
| `physical_extent_mm` | Section size | < 0.5 mm or > 50 mm per side | Likely test run or wrong pixel size |
| `n_pixels_acquired` | Spectra | < 100 | Likely test/debug acquisition: confirm before processing |

## Mass axis (`profile_mass_axis`)

| Field | Meaning | Review trigger | Action |
|---|---|---|---|
| `mz_sorted` | m/z ascending within spectra | `false` | Hold: MassFlow's pipeline rejects it; re-export |
| `representation` vs `spectrum_type` | Metadata flag vs what the spectra look like | Disagree / `null` | Trust the spectra; report the wrong flag |
| `mz_min`–`mz_max` | Acquired range | Narrower than the intended analysis (e.g. lipids need ~600–1000) | Caveat: target classes may be outside the range |
| `empty_spectra_sampled` | Spectra with no positive signal | > 5% of the sample | Off-tissue/dropout pixels: check the TIC image |
| `sampling.bin_width_da` + `uniform_bins` | Uniform bin grid width | Uniform and ≥ 0.02 Da | Re-binned data: accuracy ≈ ± half a bin |
| `sampling.points_per_fwhm` | Samples across one peak (profile) | < 3 | Sampling-limited: skip smoothing; tolerances bin-sized |
| `sampling.narrow_points_per_fwhm` (+ `narrow_peak_mz`) | Samples across the narrowest peaks (10th percentile, profile) | < ~4.3 | Note it: preprocessing narrows MassFlow's default windows (peak-pick width, smoothing window 5 → 3) so these peaks are not lowered or hidden |
| `sampling.peak_fwhm_ppm` / `resolving_power` | Measured peak width (profile) | Much worse than the analyzer class typical | Either re-binning or instrument issue |
| `n_mz_bins` (processed, aligned) | Aligned features | < 100 or > 20,000 | Too strict / too loose peak picking or alignment |

## TIC image (`render_tic_image`)

| Field | Meaning | Review trigger | Action |
|---|---|---|---|
| `row_striping`, `column_striping` | Robust line-to-line variation of line means / median | ≥ 0.08 and ≥ 1.5× the other direction | Raster artifact along the scan direction; caveat for intensity maps; normalization reduces but rarely removes it |
| `edge_to_interior_ratio` | Boundary vs interior mean TIC | ≥ 1.5 or ≤ 0.5 | Edge effect: caveat for edge pixels; distinguish from real edge biology with ion images |
| `dropout_fraction` | Pixels < 5% of median TIC | ≥ 5% | Holes, tears, missed shots or bare slide inside the raster |
| `hot_pixels` | Pixels > median + 10 robust SD | Many | Use `clip_percentile` 99 (or lower) when viewing images |
| `cv` | Coefficient of variation of TIC | > 1.0 | Strong heterogeneity or hot spots; check the image |
| `cv` with `dynamic_range_p99_p01` | Uniformity of the TIC | CV < 0.1 and p99/p01 < 1.6 | Data already TIC-normalized upstream: no second normalization; no absolute-intensity claims |
| `dynamic_range_p99_p01` | Spread of pixel TIC | > 100 | Very uneven signal; normalization choice matters |

## Processing outputs (post-processing review)

| Check | Expected | Finding if not |
|---|---|---|
| Acquired pixels, bounding box | Unchanged (or = mask preview's kept pixels) | Lost spectra |
| `representation`, `is_continuous` | centroid, true after alignment | Pipeline did not run as intended |
| Feature drop after a filter pass | Equals the preview's removal count; with both filters in one pass, the union of the two previews' removals (an `intensity_quantile` cut is recomputed after the adduct filter, so it shifts slightly) | Different data or parameters than previewed |
| Spectrum overlay at matched x/y | Stems on raw apexes; major peaks present | Peak picking too strict/loose or misalignment |
| Reference-ion offset | Same as (or better than) raw | Alignment shifted masses |
