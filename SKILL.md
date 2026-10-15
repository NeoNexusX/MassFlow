---
name: msi-dataset-evaluation
description: "Perform a complete first-level MassFlow analysis of a real imzML MSI dataset: metadata and biometa intake, spatial/mass-axis/TIC quality control, evidence-based preprocessing with required user decisions, full-dataset execution, and a concise post-preprocessing integrity review. Use whenever a new .imzML dataset is introduced or the user asks to inspect, QC, preprocess, validate, or summarize preprocessing. Do not extend into UMAP, marker-ion, molecular-annotation, literature, or anatomical interpretation unless the user separately requests downstream analysis."
---

# First-level MSI evaluation and preprocessing

Treat mzorbit tools as measurement and execution layers. Read their numeric and
image outputs, make evidence-bound decisions, preserve user-controlled choices,
and finish with one preprocessing summary. Do not equate a successful tool call
with a passed quality check.

## Scope

Complete only the first-level workflow:

1. Load and verify the real imzML dataset.
2. Evaluate spatial layout, mass axis, TIC, and biological metadata.
3. Plan, preview, verify, and execute appropriate preprocessing.
4. Reopen the processed output and review data integrity and preprocessing effects.
5. Summarize applied/skipped steps, parameters, retries, output, before/after evidence,
   and unresolved limitations.

Do not automatically run `run_umap_clustering`, `cluster_marker_ions`,
`score_ion_structure`, targeted ion images, compound annotation, or literature
search. “Post-processing review” means validating the preprocessing output, not
starting downstream spatial or biological interpretation.

## Workflow

### 1. Intake and raw-data QC

Call in this order. Stop early only when a tool itself fails or an escalation
condition below is met.

1. `load_imzml`: record continuous/processed, profile/centroid flags,
   instrument, scan pattern/type/direction, and line direction.
2. `profile_spatial`: record image shape, acquired/grid pixels, occupancy,
   bounding box, pixel size, and physical extent.
3. `profile_mass_axis`: record axis kind, representation, m/z range, bins,
   spacing, uniformity, density, and number sampled.
4. `render_tic_image` with `structure_only=True`: inspect the PNG itself.
5. `get_biometa`: collect organism, tissue, and condition. Infer tissue from TIC
   only as a tentative guess; never infer organism or condition from morphology.
   Ask one consolidated question for missing fields, then use `set_biometa`.

Cross-check the outputs:

- Flag occupancy below about 0.3, or `None` when grid dimensions are unverified.
- Flag datasets with only a few tens of spectra as likely test/debug acquisitions.
- Verify physical extent against pixel size and grid dimensions.
- Require metadata and mass-axis representation to agree. If representation is
  missing or conflicting, report ambiguity rather than assuming profile/centroid.
- Flag a non-uniform shared axis for review.
- In the TIC, distinguish tissue-like structure from blank/noisy signal. Compare
  streaking/banding orientation with scan pattern and scan direction; describe it
  as a possible acquisition/raster artifact, not a chemical feature.
- Describe edge enhancement or isolated bright pixels without assigning a cause.
  “Possible matrix/ionization/acquisition edge effect” is acceptable; “matrix
  accumulation” requires independent evidence.
- Remember that `structure_only=True` masks unacquired coordinates, not
  histologically defined non-tissue pixels.

Report the raw numbers/image evidence first, then the judgment. End intake with
`Proceed to plan_preprocessing` or `Hold - needs user input on: ...`.

### 2. Plan preprocessing

Call `plan_preprocessing` only after intake QC. Preserve its apply order and
record the complete returned plan, including sampled SNR, execute/skip reasons,
initial parameters, attempt budget, and spectrum type.

If planning fails, report the exact error and the last successful QC evidence.
Do not silently invent metadata, parameters, or a spectrum type. Ask for user
input only after exhausting tool-supported recovery.

### 3. Decide and verify each step

#### Baseline correction

For profile data, call `preview_baseline`. Inspect low/mid/high-baseline spectra,
state an evidence-based preliminary apply/skip recommendation, and then ask the
user to confirm. The user decision wins. Never skip or apply baseline correction
without recording this confirmation.

For centroid data, skip baseline correction and state that the raw profile is no
longer available for baseline estimation.

#### Noise reduction

For profile data, keep the numeric `try_step` result as the acceptance gate and
use `preview_noise_reduction` as the visual fidelity check. Do not accept or skip
denoising from images alone.

For every candidate parameter set:

1. Run `try_step` with a fixed sample size and seed across retries. Record SNR
   before/after, score, target/fidelity result, outcome, and suggested parameters.
2. Run `preview_noise_reduction` without an explicit window to obtain matched
   10, 20, and 50 m/z views of the same noisiest/typical/cleanest spectra.
3. Inspect all three scales:
   - **10 m/z:** local peak height, position, width, and high-frequency jitter.
   - **20 m/z:** neighboring peaks, shoulders, local baseline, and peak merging.
   - **50 m/z:** broader attenuation, baseline drift, relative intensity changes,
     and loss of low-abundance peaks across a wider region.
4. Reject or reduce smoothing when real peaks are clipped, shifted, broadened,
   merged, or lost at any scale, even if the numeric result passes. If the visual
   check is acceptable but `try_step` requests a retry, follow the numeric retry.
5. Repeat within the attempt budget. On final failure, report
   `escalate_to_user`; do not quietly drop the step.

If an explicit scientific m/z region matters, call the preview again with
`mz_min`/`mz_max` after the standard 10/20/50 m/z review.

For centroid data, skip noise reduction and explain why.

#### Normalization

Normalization is off by default and user-controlled. Call
`compare_normalization`, show raw/TIC/RMS previews, explain that a TIC-normalized
TIC image becomes nearly uniform by construction, and ask the user to choose
`none`, `tic_numba`, `rms_numba`, or `ref_numba` with required parameters. Do not
describe a uniform TIC-normalized image by itself as proof of over-normalization.

#### Peak picking and alignment

Follow the plan: profile data normally needs peak picking; centroid data skips
it; peak alignment runs when planned. Record resolved parameters, not merely
“defaults”. Do not infer molecular identities from the resulting peaks.

### 4. Execute the full dataset

Build `accepted_steps` as explicit objects containing each `step` and its complete
resolved `params`. Include only confirmed/accepted steps and the fixed steps
selected by the plan. Call `finalize_preprocessing` with a deterministic output
path that does not overwrite the raw input.

Record every failure and retry. If the first call fails because of malformed
`accepted_steps` or another argument error, correct it and report both the error
and recovery in the final summary.

### 5. Post-preprocessing review

Review the processed output before declaring completion:

1. `load_imzml`: confirm the output reopens and record its representation.
2. `profile_spatial`: compare acquired pixels, shape, bounding box, pixel size,
   and occupancy against the raw input. Unexpected coordinate loss is a failure.
3. `profile_mass_axis`: compare representation, range, bins/features, spacing,
   and uniformity. Explain expected changes from peak picking/alignment.
4. `render_tic_image`: compare morphology, edge enhancement, banding, blank
   regions, and new artifacts using the same display assumptions.
5. Inspect representative raw and processed spectra at matched coordinates when
   the tools allow it. Confirm that accepted steps did not create peak shifts,
   clipping, merging, or loss inconsistent with the intended operation.

Do not run UMAP or marker-ion analysis as a substitute for this review. A
visually attractive clustering result does not prove preprocessing integrity.

### 6. Final preprocessing summary

Return one concise report with these sections:

1. **Input and raw QC:** dataset identity; core spatial/mass-axis metrics; TIC
   observations; proceed/hold verdict.
2. **Preprocessing decisions:** each step, apply/skip, evidence, user decisions,
   full parameters, numeric verification, and retries.
3. **Execution:** output path, applied order, runtime errors, and recoveries.
4. **Before/after review:** spatial integrity, representation/mass-axis changes,
   TIC/spectrum changes, resolved issues, and unresolved limitations.

Avoid unsupported anatomical labels and causal claims. Say “spatially distinct
region” rather than cortex/medulla/pelvis unless registered histology or another
independent reference supports the label. Keep molecular annotation and PFAS
interpretation out of this first-level report.

## Escalation

Hold before preprocessing when:

- occupancy is `None` or below about 0.3 and TIC does not explain it;
- metadata and mass-axis representation conflict;
- the TIC has no tissue-like signal;
- spectrum type remains ambiguous and planning cannot recover;
- an output loses coordinates or cannot be reopened.

Missing biometa is a flag, not a blocker. Severe banding can proceed only as an
explicitly limited exploratory preprocessing run; state that the current spectral
steps are not a dedicated de-striping method.

## Recommended user prompt

> Please use MassFlow-Agent to inspect and preprocess `[IMZML_PATH]`. Complete
> raw-data QC, evidence-based preprocessing, and a concise post-preprocessing
> review. Ask only when a user decision is required, and report applied/skipped
> steps, parameters, retries, output, before/after changes, and unresolved issues.
> Keep this to first-level preprocessing; do not run clustering, marker-ion,
> molecular-annotation, or literature analysis.
