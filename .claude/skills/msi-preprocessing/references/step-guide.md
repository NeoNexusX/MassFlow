# Preprocessing step guide

MassFlow applies steps in this order whatever order they are registered in:
baseline_correction → noise_reduction → normalization → peak_pick →
filter_peaks_by_isotopes → peak_align → filter_peaks_by_adducts →
filter_peaks_by_stats. Profile-only steps (baseline, noise, peak_pick) are
dropped on centroid input. All defaults below are MassFlow's; "plan" means the
value `plan_preprocessing` derives from the data.

## Baseline correction (profile data)

- **What it fixes**: a slowly varying background under the peaks — chemical
  noise and matrix background (strongest at low m/z in MALDI-TOF), detector
  offset. Corrected intensities below 0 are set to 0; m/z and pixels unchanged.
- **Methods**: `locmin_numba` (default; anchors on local minima; `width` 5–9,
  `niter` 15) and `snip_numba` (iterative; better for dense peaks or complex
  backgrounds; `m` window auto). `baseline_scale` in (0, 1] subtracts only part
  of the estimate to avoid over-subtraction.
- **Skip when**: the floor between peaks already sits at zero (ratio ≪ 0.1),
  the data are centroided, or the vendor already subtracted a baseline (common
  in re-binned exports).
- **Pitfall**: over-subtraction erases small peaks and changes relative
  intensities; compare the floor, not the peak tops, in the preview.

## Noise reduction (profile data)

- **What it fixes**: high-frequency jitter that peak picking would turn into
  false peaks.
- **Methods**: `gaussian_numba` (default; `sd = window/4`), `ma_numba` (moving
  average; fast, broadens narrow peaks), `savgol_numba` (Savitzky–Golay;
  preserves peak shape; needs `polyorder ≤ window − 2` to actually smooth).
- **Plan rule**: `window ≈ samples per FWHM` (odd, 3–15). With fewer than 3
  samples per FWHM (re-binned data) every window merges neighbouring bins, so
  the plan skips smoothing — a Gaussian window 5 on 1-bin peaks kept only ~35%
  of peak height in testing.
- **Acceptance**: numeric (`try_step`: SNR gain must beat attenuation for most
  pixels; median peak retention ≥ 0.5) plus the visual veto at 10/20/50 m/z.

## Normalization

Per-pixel scaling that removes pixel-to-pixel variation in total signal
(section thickness, matrix density, ionization efficiency). It changes how
every intensity is read, which is why it is the user's choice.

| Method | Use when | Avoid when |
|---|---|---|
| `tic_numba` (total ion current) | Default convention for MALDI imaging, segmentation and visual comparison | A few very intense or very localized ions dominate the TIC (heme in vessels, a dosed drug, strong matrix peaks on bare slide): TIC normalization then imprints their pattern onto every ion |
| `rms_numba` | Spectra with many weak peaks; less sensitive to a broad noise floor than TIC | Intense outlier peaks (RMS weights them more) |
| `ref_numba` (reference peak) | An internal standard was applied (e.g. sprayed with the matrix) or a known stable ion exists; set `ref` and `ref_tolerance` (about half a bin, or 2–3× the peak width) | The reference is absent in some regions (division by noise) |
| none | Single-section exploration of absolute patterns; later quantification with calibration | Comparing intensities across regions with different total ion yield |

`scale` multiplies the result (e.g. 1000 for readable numbers). Normalization
runs before peak picking (SNR-based picking is scale-invariant), so choosing it
does not change which peaks are found. It can be applied later as its own pass
on the aligned output if deferred.

## Peak picking (profile → centroid)

- **Knobs**: `snr` (default 2; 1.5–5), `method` (noise estimate: `mad` default,
  robust to intense peaks; `quantile`, `sd`, `diff`), `width` (local-maximum
  window in samples: plan = samples per FWHM, 3 for re-binned data),
  `return_type` (`height` default; `area` only for well-sampled peaks),
  `nbins`/`overlap` (chunked noise for very long FT-ICR spectra).
- **Reading the preview**: stems on apexes, no missed visible peaks, no picks on
  baseline grass. The SNR sweep tells how many picks sit just above the noise:
  more than ~1.6× as many peaks at SNR 2 as at SNR 3 means raise `snr` to 3.
- **Typical counts**: MALDI-TOF profile spectra give hundreds to a few thousand
  peaks per spectrum before alignment.

## Isotope filter (per spectrum, between peak picking and alignment)

- Keeps the monoisotopic (lowest-m/z) peak of each isotope cluster and removes
  M+1, M+2, …, using the ¹³C spacing 1.00335/z and a decreasing-intensity check
  from M+2 (`start_intensity_check=2`).
- **Plan rules**: tolerance = a third of the peak FWHM (ppm), or one bin for
  re-binned data (in Da, capped at MassFlow's 0.1 Da / 100 ppm limits);
  `max_charge=1` for MALDI, which produces almost only singly charged ions —
  higher charges only add chances of grouping unrelated peaks.
- `min_isotope_count=3` needs M, M+1 and M+2: weak ions whose M+2 is below the
  detection limit keep their M+1 peak (conservative by design).
- Raw data keep the isotope patterns needed later for formula checks.

## Peak alignment (centroid → shared axis)

- Builds one reference m/z axis from all pixels and snaps each spectrum's peaks
  to it; unmatched positions become 0.
- `tolerance=None` estimates the matching tolerance from the data; `binfun`
  `min` (default: finest bins, longest axis), `median`, `max`; `binratio`
  scales the tolerance; `min_frequency=0.01` (default) drops reference peaks
  found in < 1% of spectra; `min_count` sets an absolute floor.
- **Signs of a bad tolerance**: too wide merges neighbouring ions (feature
  count unexpectedly low); too narrow splits one ion into near-duplicate
  features (adduct preview points at Δm/z ≈ 0, many features spaced a few ppm
  apart).
- TOF m/z drift across a section (topography, temperature) is absorbed when it
  is smaller than the tolerance.

## Adduct filter (aligned data)

- Groups features whose neutral masses agree under the CAMERA primary rules —
  positive: [M+H]⁺, [M+Na]⁺, [M+K]⁺, [M+NH₄]⁺; negative: [M−H]⁻, [M+Cl]⁻,
  [M−2H+Na]⁻, [M−2H+K]⁻ — and whose ion images correlate
  (`min_spatial_correlation`, default 0.7); keeps the most intense member.
- **Tolerance** = 2·ppm·M + mzabs. Plan: from the peak width for well-sampled
  data; for re-binned data `mzabs` = one bin (two binned peaks can differ by a
  bin) with ppm 5.
- **Polarity must match the acquisition**; the wrong rule table removes the
  wrong peaks.
- Na⁺/K⁺ ratios vary across tissue, so true adduct pairs can correlate below
  0.7 — the filter is conservative and keeps them. Removed adducts lose
  adduct-ratio maps, which can be biology (ion homeostasis, damage): keep the
  aligned file.
- Read `removal_fraction` (typically a few percent), `removed_by_rule`, and the
  borderline examples; the audit CSV lists every removal with its evidence.

## Statistical peak filter (aligned data)

- Keeps features detected in at least `min_frequency` of spectra (or
  `min_count` spectra) AND with aggregate intensity above `min_intensity` or the
  `intensity_quantile` (strictly greater than the quantile value).
- Alignment already removed features below 1% frequency; start at 0.02.
- Never use a cut-off above the smallest region of interest's share of pixels:
  a 3% structure disappears at `min_frequency=0.05`.

## Tissue mask (aligned data)

- UMAP + KMeans clusters pixels; off-tissue background (bare slide, matrix) is
  a cluster framing the tissue (`exterior_fraction` near 1) whose spectra carry
  little signal above the middle of the m/z range (`high_mass_ratio` ≤ 0.7,
  matrix clusters dominate). Its TIC is often *higher* than tissue — matrix
  ionizes without tissue suppression — and TIC-normalized data have no TIC
  contrast at all, so `tic_ratio` is only a fallback. In testing, a k35 matrix
  margin had 2.2× the tissue TIC but a 0.49 high-mass ratio. A frame-like
  cluster with normal high-mass signal is tissue (edge effect or adjacent
  tissue).
- Modes: `remove_clusters` (explicit labels; most predictable),
  `remove_background` (peels outer layers), `keep_clusters` (keeps labels plus
  everything they enclose).
- Mask before statistics or segmentation of the tissue; per-pixel
  normalizations (TIC/RMS) are unaffected by masking.
- Tissue-shaped rasters (low bounding-box fill) usually have no background to
  remove.

## Backends

`backend="cardinal"` (R Cardinal) exists for baseline, noise, normalization,
peak picking and alignment, in `finalize_preprocessing` only, and needs R with
Cardinal installed (MassFlow `docs/en/r_environment_setup.md`). Previews and
`try_step` verify the python backend; say so if the user picks Cardinal.

## Export formats

| Format | Layout | Use |
|---|---|---|
| `ion_image_zarr` | Ion-major images + stats (mean spectrum, RMS) | Fast ion images, `massflow.morch` deep-learning datasets |
| `ion_image_h5` / `ion_image_binary` | Ion-major | Same, for HDF5/binary consumers |
| `spectra_zarr` | Pixel-major spectra + TIC/RMS | Spectrum-wise access; works for unaligned data |
