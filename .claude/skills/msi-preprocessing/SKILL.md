---
name: msi-preprocessing
description: "Preprocess an imzML mass spectrometry imaging dataset with MassFlow through the mzorbit tools: plan with MassFlow's default parameters (changed only on evidence), decide baseline/smoothing/normalization with previews, pick and align peaks, optionally filter isotopes, adducts and rare peaks, remove off-tissue background with a tissue mask, write the processed imzML with provenance, and export Zarr/ion-image formats. Works in guided mode (asks at real decisions) or hands-off autonomous mode, and recovers from errors automatically. Use when the user asks to preprocess, peak-pick, align, normalize, denoise, filter peaks, remove background, or prepare MSI data for clustering, statistics or annotation (预处理, 峰提取, 峰对齐, 归一化, 去噪, 去背景, 数据处理). Evaluate the dataset first with msi-dataset-evaluation if that has not been done; hand interpretation of the processed data (regions, marker ions, compounds) to msi-spatial-analysis."
---

# MSI preprocessing

Turn a raw imzML dataset into analysis-ready data, with every decision backed by
evidence and every output traceable. The mzorbit tools measure, preview and
execute; this skill decides. Raw data are never modified: each pass writes a
new imzML whose sidecar records the steps, parameters, inherited polarity and
biometa.

Step-by-step domain guidance (what each step does to MSI data, its knobs, and
its pitfalls) is in `references/step-guide.md`; failure handling is in
`references/error-recovery.md`.

## Before you start

1. **Evaluation**: if this dataset has not been evaluated in the session, run
   the intake of the `msi-dataset-evaluation` skill (at least `load_imzml`,
   `profile_mass_axis`, `render_tic_image`). Do not process a dataset whose
   verdict is Hold.
2. **Goal**: know what the data are for — exploratory segmentation, targeted
   ions, feature tables for statistics, or export for deep learning. The goal
   decides which optional steps matter. Ask in one line if it is unclear.
3. **Mode**:
   - **Guided** (default): decide technical parameters yourself; ask the user
     at the semantic decisions marked *user* below, after stating your
     recommendation and evidence.
   - **Autonomous**: when the user asks for hands-off processing ("全自动",
     "直接处理", "don't ask me"), decide everything with the default policies
     below, never block, and report each decision with its evidence and how to
     change it. The tools ask for the user's choice by default; a hands-off
     request delegates those choices to you, so mark them `delegated`.
4. **Outputs**: write next to the input unless told otherwise:
   `<stem>_aligned.imzML` (pass 1), `<stem>_filtered.imzML` (pass 2),
   `<stem>_tissue.imzML` (mask). Never reuse the input path; avoid dots in
   names (the tool rewrites them anyway).

## Decision ownership

| Step | Guided | Autonomous default | Evidence |
|---|---|---|---|
| baseline_correction | model proposes, **user** confirms | apply only if `baseline_to_peak_ratio` ≥ 0.1 **and** the floor gap is visible | `preview_baseline`, then `try_step` |
| noise_reduction | plan + `try_step` (window 3–9) | same | plan reason, `try_step(auto_tune)`, `preview_noise_reduction` |
| normalization | **user** chooses | none; `tic_numba` only when the goal is clustering/segmentation | `compare_normalization` |
| peak_pick | model: plan params (MassFlow defaults), changed only for a preview problem | same | `preview_peak_pick` |
| filter_peaks_by_isotopes | offer when redundancy matters; **user** | apply when the goal needs a feature table/clustering | plan params |
| peak_align | auto (skipped if already aligned) | same | plan |
| filter_peaks_by_adducts | **user** after preview | apply if polarity comes from the file/user and the preview removes ≤ 25% with evidence above threshold | `preview_adduct_filter` |
| filter_peaks_by_stats | **user** after preview | apply the smallest cut-off ≥ 0.02 that removes ≤ 30% | `preview_peak_filter` |
| tissue mask | **user** after preview | apply only the suggested background, keeping ≥ 40% of pixels | `run_umap_clustering`, `preview_tissue_mask` |
| export | on request | on request | — |

The user's explicit instruction always overrides a default. A plan step with
`decision: "auto"` and `execute: false` (e.g. smoothing skipped on re-binned
data) is a finding, not a question: report it, don't ask.

## Parameter policy

Preprocessing should keep as many real features as possible; removing noise
is the job of the later steps (alignment's `min_frequency`, the peak filters,
the tissue mask).

- **MassFlow's defaults first.** The plan uses them (`peak_pick` width 5,
  snr 2; Gaussian smoothing window 5) and only ever narrows a window, for
  narrow peaks: under-sampled peaks (< 3 samples per FWHM) get width 3 and no
  smoothing; when the narrowest peaks span fewer than ~4.3 samples per FWHM,
  width and window drop to 3. It never widens them — a wide window swallows
  small peaks next to taller ones.
- **Picking is for recall.** Many small or near-noise picks are expected.
  Never raise `snr`, `width` or the smoothing just to get fewer peaks.
- **Change a parameter only on evidence**: a problem a preview reports, a
  `try_step` failure, or something clearly wrong in the images (a real peak
  without a stem, peaks flattened or merged). Say what the evidence was. Keep
  smoothing windows within 3–9; when in doubt, keep the default.

## Workflow

### 1. Plan

Call `plan_preprocessing`. Read, in this order:

- `context.acquisition` — polarity and its source, analyzer, ion source,
  `points_per_fwhm` and the narrowest peaks' `narrow_points_per_fwhm`,
  `bin_width_da`, `sampling_limited`. These decide whether a default window is
  narrowed and set the filter tolerances.
- `context.notes` and `warnings` — automatic repairs.
- `auto_steps`, `verify_steps`, `user_decisions` — the work list.
- Each step's `params`, `reason` and `evidence`. Use the params as given unless
  a preview shows a reason to change them; when you change one, say why.

### 2. Profile-only steps (profile data; centroid input skips this section)

1. **Baseline** — `preview_baseline`. Read the floor gap in the images and
   `baseline_to_peak_ratio`. State your call ("apply"/"skip") with both pieces
   of evidence, then follow the mode. When the decision is "apply", verify it
   with `try_step(step="baseline_correction", params=<plan params>,
   auto_tune=true)`: the preview judges the floor, this checks the peaks — on
   a sample the floor must drop by at least half while peak heights keep
   ≥ 95%, and gentler settings are tried if peaks suffer (over-subtraction
   erases small peaks). Carry its accepted params forward; if it escalates,
   skip baseline (autonomous) or show the user the attempt history (guided).
2. **Noise reduction** — if the plan skipped it because peaks are
   under-sampled, report that and move on. Otherwise run
   `try_step(step="noise_reduction", params=<plan params>,
   accepted_steps=[<accepted baseline step, if any>], auto_tune=true)` (fixed
   seed), so the sample is smoothed after the same baseline finalize will
   apply. Smoothing is accepted only if the typical peak keeps ≥ 85% of its
   height; when peaks suffer the retries narrow the window (down to 3) before
   switching to Savitzky–Golay, and an under-smoothed sample widens only up to
   9 and only as far as the narrowest peaks allow. If it returns `accepted`,
   confirm with `preview_noise_reduction(params=<accepted params>)` at the
   default 10/20/50 m/z windows; veto if peaks are clipped, shifted, broadened
   or merged at any scale. If it escalates, skip denoising (autonomous) or ask
   (guided).
3. **Normalization** — `compare_normalization`; show the three images and the
   `uniformity` numbers; explain TIC vs RMS vs reference-peak (see the step
   guide); follow the mode. TIC-normalized TIC images are flat by construction —
   that is not over-normalization. If evaluation flagged the raw TIC as already
   nearly uniform, the data were normalized upstream: recommend none.
4. **Peak picking** — `preview_peak_pick(accepted_steps=<accepted baseline/noise
   steps>)` with the plan params. The aim is to keep every candidate peak, so
   use `finalize_step` unless `stats.problems` lists one of these:
   - `window_check` `too_wide` — the window is wider than MassFlow's default
     (fewer candidates kept) or can hide resolved neighbours of the narrowest
     peaks: use `suggested_params` (the widest safe width; the basis compares
     picks per spectrum); the `narrow_width…` close-ups show the narrowest
     peaks picked with it.
   - `noise_check` `collapsed` — the noise estimate is zero (often
     zero-filled bins), so every positive point passes the SNR test: preview
     again with a `prominence` a few times the typical grass height and
     compare peaks per spectrum.

   Also look at the images: a clearly visible peak without a stem means
   picking is too strict (lower `width` or `snr` and preview again); stems on
   shoulders instead of apexes mean the peaks were distorted (revisit the
   smoothing). Many picks near the noise (`snr_sweep`) or far below the base
   peak (`weak_pick_fraction`) are normal in MSI and not a reason to change
   anything: alignment screens them.

### 3. Pass 1 — core chain

`finalize_preprocessing(imzml_path=<raw>, output_path=<stem>_aligned.imzML,
accepted_steps=[...])` with, in any order (the tool applies MassFlow's order):
accepted baseline/noise/normalization steps, `peak_pick`, optional
`filter_peaks_by_isotopes` (plan params), and `peak_align` (plan params).
Check `output_summary` (pixel count unchanged, `centroid`, aligned feature
count) and every warning. Alignment keeps MassFlow's defaults: its
`min_frequency=0.01` drops picks that do not recur across pixels, which is
where most near-noise picks go. Set its `tolerance` only when the output
shows a problem (one ion split into features a few ppm apart, or neighbouring
ions merged; see the step guide).

Input that is already aligned (plan `context.is_aligned`) skips pass 1.

### 4. Pass 2 — peak-list filters (optional, on the aligned output)

Both filters need aligned data, so they are previewed on the pass-1 output:

- **Adducts** — `preview_adduct_filter(imzml_path=<aligned>)`; polarity and
  tolerances come from the inherited facts (`param_sources` says where).
  Judge `removal_fraction`, `removed_by_rule`, and how close
  `removal_correlation` sits to the threshold. Re-preview with other
  `min_spatial_correlation`/tolerances if needed.
- **Rare/weak peaks** — `preview_peak_filter(imzml_path=<aligned>,
  min_frequency=…)`; use `removed_by_frequency_cutoff` to choose. `peak_align`
  already drops peaks seen in < 1% of spectra (its default `min_frequency=0.01`),
  so a stats cut-off of 0.01 removes nothing; meaningful cut-offs start
  around 0.02. Never set a frequency cut-off above the smallest region of
  interest's share of pixels.

Then `finalize_preprocessing(imzml_path=<aligned>,
output_path=<stem>_filtered.imzML, accepted_steps=[<finalize_step of each
accepted preview>])`. The feature drop must match the previews: one filter
removes exactly its preview's count; both filters remove the union of the two
previews' sets (an `intensity_quantile` cut is recomputed after the adduct
filter, so it shifts slightly). Any other difference means the data or
parameters differ from what was previewed.

### 5. Tissue mask (optional)

Useful when a rectangular raster includes off-tissue area and the goal involves
normalization, statistics or clustering (background pixels distort all three).
Skip it for tissue-shaped rasters (low `bounding_box_fill`, no background
frame).

1. `run_umap_clustering(imzml_path=<latest aligned output>, random_state=0)`.
   Read `clusters` (`exterior_fraction`, `high_mass_ratio`, `tic_ratio`) and
   `suggested_background`. Background is a frame-like cluster with little
   signal above the middle of the m/z range (matrix-dominated spectra,
   `high_mass_ratio` ≤ 0.7) — often *brighter* in the TIC than tissue, so never
   decide by TIC alone. A frame-like cluster with normal `high_mass_ratio` is
   tissue (an edge effect or adjacent tissue): keep it. These clusters only
   serve the mask; they include the background, so do not interpret them
   biologically — segmentation belongs to msi-spatial-analysis, run on the
   masked output.
2. `preview_tissue_mask(result_path=<result_npz>)` (mode `auto` uses the
   suggestion) or with explicit `mode`/`labels`. Compare the masked image with
   the TIC image: the whole section must survive.
3. `apply_tissue_mask(..., output_path=<stem>_tissue.imzML,
   **apply_arguments)`.

### 6. Export (on request)

`export_dataset(output_format=…)`: `ion_image_zarr` (default; fast ion images,
`massflow.morch` datasets), `ion_image_h5`/`ion_image_binary`, or
`spectra_zarr` (pixel-major, also for unaligned data). TIC/RMS/mean spectrum
are precomputed into the store.

### 7. Review and report

Run the **post-processing review** of `msi-dataset-evaluation` on the final
output (raw vs processed: pixels, representation, features, TIC morphology,
spectrum overlays, reference ions). Then report, in the user's language:

1. **Input** — dataset, verdict and caveats from evaluation, key acquisition
   facts (polarity + source, bin width / peak sampling).
2. **Decisions** — table: step, applied/skipped, owner (auto / verify / user /
   delegated), key evidence, final parameters; name every departure from
   MassFlow's defaults with the evidence behind it.
3. **Outputs** — each pass's path, pixel and feature counts, audit CSVs,
   warnings and automatic recoveries.
4. **Review** — what changed between raw and processed, what was verified,
   remaining limitations.
5. **Reproduce** — the exact `accepted_steps` JSON of each pass.
6. **Next step** — for regions, marker ions or compound questions, continue
   with the msi-spatial-analysis skill on the final output, carrying the goal
   and the evaluation caveats.

Stay within preprocessing: interpreting clusters, naming molecules and
literature work belong to the msi-spatial-analysis skill.

## Error handling (summary)

- Read `details.recovery` on any failure and apply it once; if the same call
  fails twice, stop and report both errors.
- `warnings` are repairs already made (tolerant loading, geometry from
  coordinates, renamed outputs, polarity filled from facts, pipeline retried
  with a smaller batch). Report them; don't redo them.
- Never fix a failure by overwriting the input, inventing polarity, or silently
  dropping a step the user accepted.

Full table: `references/error-recovery.md`.
