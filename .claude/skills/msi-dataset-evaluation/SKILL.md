---
name: msi-dataset-evaluation
description: "Evaluate an imzML mass spectrometry imaging (MSI) dataset with the mzorbit tools: file integrity, acquisition facts (ion polarity, analyzer, ion source, m/z sampling and binning), spatial layout, mass axis, TIC artifacts (striping, edge effects, hot pixels, dropouts), biological metadata and a reference-ion calibration check, ending in one Proceed / Proceed with caveats / Hold verdict. Also reviews a processed output against its raw input. Use when a new .imzML arrives, or the user asks to inspect, QC, sanity-check, assess or compare MSI data (数据评估, 质控, 质量检查, 数据检查, 处理结果检查). Read-only: hand any processing to the msi-preprocessing skill, and segmentation, marker ions or compound annotation to msi-spatial-analysis."
---

# MSI dataset evaluation

Answer one question with measured evidence: **what is this dataset, and is it fit
for the intended analysis?** The mzorbit tools measure and flag; this skill
decides what the measurements mean. Evaluation never changes data — processing
belongs to the `msi-preprocessing` skill, which calls back here for its final
review.

Domain background (instruments, matrices, reference ions, artifact catalogue)
is in `references/msi-domain-knowledge.md`; every metric with its review
trigger is in `references/qc-metrics.md`. Read them when a finding needs
interpretation beyond the rules below.

## Principles

- A successful tool call is not a passed check. Read the numbers, the flags and
  the image itself.
- Report evidence first, then the judgment, with its uncertainty.
- Use the tool's `flags` as the starting list of findings, not the whole list:
  look at the images for what metrics cannot see (tissue folds, smears, missing
  regions).
- `warnings` on any result are automatic repairs (tolerant loading, image size
  taken from coordinates, dot-free links). Report them; they are not errors.
- Never invent facts. Polarity, organism, tissue and condition come from the
  file, the user, or explicit evidence; label every inference as such.
- Describe artifacts as possible causes ("consistent with laser-energy drift"),
  not certainties, unless independent evidence exists.

## Intake workflow

Run these in order; stop early only if a tool fails (then follow
`details.recovery`).

1. **`load_imzml`** — record storage (continuous/processed), representation
   (profile/centroid), analyzer class and ion source, polarity **with its
   source** (`imzML cvParam`, `user`, `inherited from …`, `file name`), file
   sizes, and `processing_history` (a non-empty history means this is a
   processed output: note the steps and where facts were inherited from).
   - Polarity unknown or only from the file name, and adduct-related work is
     planned → ask the user once, then `set_acquisition`.
2. **`profile_spatial`** — image shape, acquired pixels, occupancy,
   `bounding_box_fill`, pixel size, physical extent, duplicate coordinates,
   z levels, flags. Check the extent is plausible for a tissue section
   (typically 2–25 mm per side) given the pixel size.
3. **`profile_mass_axis`** — m/z range, sortedness, empty spectra, and
   `sampling`: bin width, peak FWHM, samples per FWHM, `sampling_limited`.
   Read `recommended_parameters` only as context; the processing skill acts
   on them.
4. **`render_tic_image`** — look at the PNG, then `quality` and `flags`:
   striping (row/column), edge-to-interior ratio, dropouts, hot pixels, CV.
   Use `clip_percentile` (default 99) if hot pixels dominate the colour scale.
5. **`get_biometa`** — for missing fields, infer the organ from TIC morphology
   only as a labelled guess (e.g. "outline and cortex/medulla-like banding are
   consistent with kidney"); never infer organism or condition. Ask one
   consolidated question, then `set_biometa`. Missing biometa is a flag, not a
   blocker.
6. **Targeted checks** — cheap and decisive when something is unclear:
   - `render_spectrum` at a high-TIC and a low-TIC pixel (coordinates from the
     TIC image) to see real spectra; centroid data are drawn as stems.
   - `render_ion_images` for one expected matrix ion and one expected tissue
     ion (see the reference-ion table) to see on/off-tissue contrast and the
     acquisition footprint.
   - **Calibration and polarity check** (centroid data only): run
     `annotate_adducts` with a 3–6 compound reference panel suited to the
     tissue and polarity (`references/msi-domain-knowledge.md`), passing
     `imzml_path`, the `ion_source`, and a `ppm_tol` that reflects the
     data's real accuracy — about 2–3× the analyzer's typical error (TOF
     ~50, Orbitrap ~10, FT-ICR ~5, unknown ~25 ppm), or half a bin in ppm
     at the lowest panel m/z for re-binned data. The tool always returns
     the nearest peak and `ppm_tol` only sets the `within_tol` flag (its
     default, 100 ppm, means little on TOF or Orbitrap data), so judge the
     `best_ppm` values together: several panel ions sharing one offset confirm
     polarity and calibration (report the offset); a shared large offset
     means calibration drift; scattered values mean no real match — wrong
     polarity, a different m/z range, or unusual chemistry. When the
     polarity comes only from the file name or is unknown, pass
     `polarity="both"` so positive and negative adducts compete; with
     `auto` the check can only test the polarity already assumed. For
     profile data, defer this check to the post-processing review.

## Interpretation rules

- **Representation**: metadata and spectra must agree. A flag "metadata says X
  but the spectra look Y" means the flag is wrong: evaluate as Y and say so.
- **m/z sortedness false** → Hold: MassFlow's pipeline rejects the data; it must
  be re-exported.
- **Sampling-limited / re-binned profile** (≈1 sample per FWHM on a uniform bin
  grid): not a defect of the sample, but it caps mass accuracy at about one
  bin, makes smoothing harmful, and means isotope/adduct tolerances must be
  bin-sized. Report it as a caveat with the bin width.
- **Occupancy**: low occupancy with low `bounding_box_fill` is a tissue-shaped
  (irregular) acquisition region — normal for many instruments. Low occupancy
  with a high fill means the declared grid is larger than the raster (often
  harmless). Duplicate coordinates are a real problem (overlapping rasters or a
  lost z dimension).
- **TIC artifacts**: striping aligned with the scan direction is an acquisition
  artifact (laser energy/raster timing); edge ratios far from 1 suggest matrix
  accumulation or ablation differences at section edges — or an off-tissue
  matrix margin inside the raster; dropouts suggest missed shots, holes or bare
  slide inside the raster; hot pixels need colour clipping. None of these are
  chemistry. Preprocessing has no de-striping step (normalization reduces
  striping but rarely removes it), so strongly striped data can only be
  processed as an exploratory run: say so in the verdict.
- **Nearly uniform TIC** (CV < 0.1, p99/p01 < 1.6, flagged): the data were
  almost certainly TIC-normalized before export. Record it: normalizing again
  is pointless, absolute intensities are not comparable, and TIC cannot locate
  background.
- **Test acquisitions**: fewer than ~100 spectra or a physical extent under
  ~0.5 mm usually means a test or debug run — ask before investing in it.
- **Off-tissue pixels**: rasters often include a margin of bare slide around
  the section, even tissue-shaped ones. Matrix there ionizes without tissue
  suppression, so the margin is frequently *brighter* in the TIC than the
  tissue; it is recognised by its spectra (matrix ions, little signal above the
  middle of the m/z range), which `run_umap_clustering` reports per cluster as
  `high_mass_ratio`. A bright rim is therefore not proof of an edge effect.

## Verdict

End intake with exactly one verdict:

- **Proceed** — no hold condition; every flag explained and harmless for the
  intended analysis.
- **Proceed with caveats** — usable, with limitations that must travel with the
  data (e.g. re-binned spectra, striping, edge effect, polarity only from the
  file name, missing biometa).
- **Hold — needs …** — any of: unsorted m/z; a tool failure that recovery
  could not resolve; duplicate coordinates without an explanation; no
  tissue-like signal in the TIC; unresolved representation conflict; a test-size
  acquisition the user has not confirmed.

## Post-processing review

Used by `msi-preprocessing` after `finalize_preprocessing`, or on request.
Compare the processed output with its input; expected changes are fine,
unexplained ones are findings.

1. `load_imzml(output)` — `processing_history` lists the intended steps;
   polarity and biometa were inherited.
2. `profile_spatial(output)` — acquired pixels and bounding box unchanged
   (unless a tissue mask was applied: then the kept fraction must match its
   preview).
3. `profile_mass_axis(output)` — representation `centroid`, `is_continuous`
   true after alignment, feature count (`n_mz_bins`) plausible (typically a
   few hundred to a few thousand for lipid/metabolite MALDI data). A filter
   pass must remove exactly what its preview reported.
4. `render_tic_image(output)` — same morphology; no new striping, holes or edge
   effects introduced by processing. Absolute TIC values change after peak
   picking (sum of peak heights) — compare patterns, not numbers.
5. `render_spectrum(raw, x=…, y=…, compare_imzml_path=output, mz_min=…,
   mz_max=…)` at 2–3 pixels (high, typical and low TIC) over a 50–200 m/z window
   with major peaks: stems should sit on raw apexes, no major peak missing or
   shifted by more than about half a bin / the alignment tolerance.
6. Reference-ion check (as in intake step 6) on the centroid output.
7. Verdict for downstream use, as above.

## Report

Write in the user's language. Keep it short and evidence-led:

1. **Dataset identity** — file, storage/representation, instrument, polarity
   (with source), pixel grid and size, m/z range, biometa.
2. **Measurements** — the key numbers from each profile (table).
3. **Findings** — each flag or visual finding: evidence, likely cause, impact.
4. **Verdict** — one of the three, with the conditions or caveats.
5. **Next step** — usually "preprocess with msi-preprocessing" with the
   parameters the facts imply (e.g. "skip smoothing: re-binned 0.12 Da data").

Do not interpret clusters, assign molecular identities, or draw biological
conclusions here; those are downstream analyses the user requests separately
(the msi-spatial-analysis skill, on processed data).
