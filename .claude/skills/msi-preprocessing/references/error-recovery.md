# Error recovery

Every failed tool result carries `errors`, `details.error_type`, and — when the
cause is recognised — `details.recovery` with the next step. Repairs the tools
already made appear in `warnings` on successful results. Apply a recovery once;
if the same call fails twice, stop and report both errors with what you tried.

## Repaired automatically (report, don't redo)

| Warning says | What happened | Your follow-up |
|---|---|---|
| "header is incomplete … tolerant loader" | Missing scanSettings; opened with MassFlow's `risky_imzml_loader` | Check pixel size and scan settings with the user if they matter |
| "`max_count_of_pixels_*` was … using N" | Image size missing or smaller than the coordinates | None; mention it in the evaluation |
| "Spectrum type flags were unusable … inferred" | Contradictory profile/centroid flags | Confirm the inferred type looks right in a spectrum |
| "Opened … through dot-free links" | Dots in the file name | None |
| "Output name … has dots … writing …" | Output renamed so MassFlow writes where expected | Use `artifacts.output` as the real path |
| "polarity set to … (source)" | Adduct step lacked polarity; filled from facts | Check the source is trustworthy (file/user over file name) |
| "Pipeline TimeoutError/MemoryError … retried with batch_size=…" | Stall or out of memory; retried smaller | If it recurs, lower `max_threads` too |
| "Centroid input: MassFlow skipped profile-only step(s)" | Baseline/noise/peak_pick sent for centroid data | Remove them from future step lists |
| "Pixel count changed" | Output lost spectra | Treat as a finding; check the input for empty spectra |

## Failures and what to do

| Error / recovery text | Cause | Action |
|---|---|---|
| "imzML file not found … imzML files in …" | Wrong path | Pick the intended file from the listed candidates or ask |
| "No binary data file next to …" | `.ibd` missing or renamed | Ask for the `.ibd`; both files must share a base name |
| "… is empty" | Truncated copy | Ask the user to re-copy the dataset |
| "m/z values must be sorted in ascending order" | Unsorted export | Hold; the dataset must be re-exported |
| "Ion polarity is unknown …" / polarity errors | No cvParam, no recorded fact | Ask the user → `set_acquisition` → retry |
| "needs aligned centroid data … run peak_align first" | Adduct/stats filter or ion-image export on unaligned data | Run pass 1 (`peak_pick` + `peak_align`) and use its output |
| "needs per-spectrum peaks, but this dataset is already aligned" | Isotope filter on aligned data | Drop it, or rerun pass 1 from the raw file with the isotope filter |
| "needs centroid peaks; include peak_pick" | Filter on profile data without peak picking | Add `peak_pick` (and `peak_align` for adduct/stats) |
| "Unknown preprocessing step 'X'; valid steps: …" | Misspelled step | Use one of the listed names |
| "params for X must be an object" | Params passed as a list/string | Pass a JSON object |
| "output_path would overwrite the input dataset" | Output = input | Choose `<stem>_aligned.imzML` etc. |
| "Async preprocessing pipeline made no progress" (after the retry) | Persistent stall | Lower `max_threads`, check free disk space, retry once |
| "cardinal" / "rpy2" errors | R/Cardinal not installed | Use the python backend (default) |
| `try_step` outcome `escalate_to_user` | No parameter set improved SNR without damaging peaks | Skip noise reduction (autonomous) or ask with the attempt history |
| preview returns `applicable: false` | Step not relevant for centroid data | Skip the step |
| "UMAP/KMeans needs a shared m/z axis" | Clustering unaligned data | Cluster the aligned output |
| "Clustering result not found" | Wrong `result_path` | Use `artifacts.result_npz` from `run_umap_clustering` |
| "mask does not retain any source pixels" | Labels remove everything (often `keep_clusters` with background labels) | Re-preview with corrected labels |
| "n_clusters" errors | Out of range | Use 1–20 |
| "Comparing two datasets needs an x/y coordinate" | `render_spectrum` compare by index | Pass `x` and `y` (indices differ between files) |
| "No spectrum at that pixel" | Coordinate not acquired (or masked out in the second file) | Pick a coordinate from `profile_spatial`'s bounding box inside the tissue |

## Never

- Overwrite or modify the raw input.
- Guess the polarity, organism, tissue or condition to get past an error.
- Silently drop a step the user accepted; report why it could not run.
- Loop on the same failing call with unchanged inputs.
