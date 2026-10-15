---
name: msi-spatial-analysis
description: "Analyze preprocessed (peak-picked, aligned) imzML mass spectrometry imaging data with the mzorbit tools: UMAP/KMeans segmentation and region interpretation, cluster marker ions, ion-image structure scores, and targeted or literature-driven compound searches (PubMed → PubChem monoisotopic mass → adduct m/z → the dataset's peaks → ion images), reported with ppm errors and annotation confidence. Use whenever the user asks which regions a section has, which ions mark a region, whether a compound, lipid class, drug or exposure-related molecule is present and where, what an m/z could be, or how the data relate to the literature (空间分割, 聚类, 区域, 标志离子, 特征离子, 离子图像, 化合物注释, 代谢物鉴定, 脂质注释, 文献检索, 靶向分析). Needs aligned data: preprocess first with msi-preprocessing; data QC belongs to msi-dataset-evaluation."
---

# MSI spatial analysis

Answer **what is in this section, and where?** with evidence that holds up.
The mzorbit tools compute — clusters, enrichment, structure scores, adduct m/z,
nearest peaks — and this skill decides what those numbers support and how
strongly to say it. Analysis writes no new datasets: processing belongs to
`msi-preprocessing`, data QC to `msi-dataset-evaluation`.

Read `references/annotation-evidence.md` before naming any compound or
choosing a tolerance: confidence levels, tolerances, adduct coverage,
look-alike ions, contaminants, and how to read structure scores. Matrix ions,
lipid classes by polarity and the reference-ion panel are in
`../msi-dataset-evaluation/references/msi-domain-knowledge.md`.

## Principles

- **Tools do the arithmetic.** Every m/z comes from `resolve_compound_mass` →
  `annotate_adducts`. A mass recalled from memory is a guess, and a 0.01 Da
  slip moves a match by more than 10 ppm.
- **A mass match is a candidate, not an identification.** Without MS/MS or a
  standard, this data supports at best a formula-level annotation, and any
  name is one of the isomers of that formula (PC 34:1 and PE 37:1 share one).
  Write names as candidates, with formula, adduct and ppm.
- **Exploration speaks in classes.** Describe regions by molecular class
  (phospholipids, fatty acids, small metabolites), as the tool guidance asks.
  Specific names appear only for hypotheses the user or the literature
  supplies, or when the user asks what an m/z could be.
- **"Where" needs an image you looked at.** Rankings (enrichment, effective
  pixels) nominate ions; the ion image decides.
- **Not detected is not absent.** Ion suppression, detection limits, the
  polarity, the m/z range or a preprocessing filter can hide a compound.
- **Describe regions; don't name them.** Use position and shape ("outer band",
  "central core") or the cluster number. Anatomical names (cortex, white
  matter, glomeruli) need histology, a registered annotation or the user's
  labels; otherwise offer them only as a labelled hypothesis.
- **One section is descriptive.** Differences between regions, or between two
  datasets, are observations; without replicate sections there is no
  statistical claim.
- Report `warnings` (repairs already made); apply `details.recovery` once.

## Before you start

1. `load_imzml(<dataset>)` — read `processing_history` (peak picking,
   alignment, filters, mask, normalization), polarity with its source, ion
   source and biometa. Work on the latest processed output: its sidecar is
   the handoff from preprocessing.
2. **Aligned data.** Clustering and markers need centroid data on a shared
   m/z axis (`profile_mass_axis`: `spectrum_type` centroid, `is_continuous`
   true); `annotate_adducts` refuses profile data. Otherwise say that pass 1
   of `msi-preprocessing` is needed, and stop.
3. **Polarity.** Adduct sets depend on it. Unknown or only from the file name
   → ask once, then `set_acquisition`; otherwise the tool tries both
   polarities, which invites chance matches — say so.
4. **Background.** If the raster includes bare slide (rectangular raster, an
   off-tissue margin flagged in evaluation) and no mask was applied, recommend
   the tissue-mask step of `msi-preprocessing` first: background pixels form
   their own clusters, enter every marker average and dilute structure scores.
   If the user declines, keep `suggested_background` clusters out of the
   interpretation.
5. **Question.** Exploratory (which regions, which markers), targeted (is X
   here, and where) or literature-driven (topic → compounds → this data). Ask
   in one line if unclear. Carry the evaluation caveats forward (striping,
   edge effect, re-binned spectra, normalization): they limit what any finding
   can mean.

## Exploratory: regions and their markers

1. **Cluster.** `run_umap_clustering(imzml_path, n_clusters=…, random_state=0,
   output_dir=<stem>_clusters_k<n>_s<seed>)`. Give every run its own
   `output_dir`: the default folder and file names are fixed, so a second run
   would overwrite the first `result_npz`. Start with the default 10, fewer
   for small or simple sections.
2. **Read the cluster image beside the TIC image.** For each cluster note its
   `fraction`, location, and whether it is contiguous or scattered. Flag
   technical clusters before any biology:
   - background: `suggested_background` (frame-like, low `high_mass_ratio`);
   - frame-like with normal `high_mass_ratio`: tissue at the section edge — an
     edge effect or real peripheral tissue; check the TIC edge ratio and ion
     images before treating it as biology;
   - stripes along the scan direction: raster artifact;
   - clusters that reproduce TIC bands in unnormalized data: they track total
     signal, not composition (a normalization pass on the aligned output is
     the remedy).
3. **Stability.** Re-run with another k (±2–3) or `random_state=1`. A region
   that persists is a finding; one that splits or merges with k is a
   clustering choice. Report which is which.
4. **Markers.** `cluster_marker_ions(imzml_path, result_path=<result_npz>,
   top_n=10)`. Enrichment is the cluster's mean over the mean of all cluster
   means, so an ion found in only one cluster saturates at `n_clusters`
   however weak it is, and background clusters take part in the average. The
   list nominates; it does not prove.
5. **Structure.** `score_ion_structure(imzml_path, mzs=<candidate markers>)`.
   Keep ions that are spatially coherent (Moran's I clearly above 0, roughly
   ≥ 0.3) with enough strong pixels; drop salt-and-pepper ions (I ≈ 0) and
   hot-spot ions (few effective pixels, very high `dynamic_range`). The ranking
   favours large, diffuse ions, so look at the small, specific ones as well.
6. **Look.** `render_ion_images(imzml_path, mzs=<3–6 per region>)`. Does the
   ion follow the region? Is it bright on bare slide (matrix, contamination)
   or striped (raster)? Lower `clip_percentile` when hot spots flatten the
   colour scale.
7. **Describe** each non-technical region: size, position and shape,
   stability, its structured markers (m/z, enrichment, Moran's I) and their
   likely class from mass range and polarity. If the user wants to know what a
   marker could be, continue with the targeted workflow.

## Targeted and literature-driven: compounds

1. **Candidates.**
   - From the user: names, or a class — then list representative species of
     that class that ionize in this polarity.
   - For an m/z the user asks about: propose several candidates — another
     class with the same formula, another adduct (the [M+Na]⁺ of one lipid vs
     the [M+H]⁺ of another), the isotope peak of a neighbour. One candidate in
     tolerance does not rule out the others.
   - From the literature: `search_pubmed(query, retmax=15–20)` with the topic,
     the organism/tissue from biometa and method terms, e.g. `PM2.5 AND lung
     AND (lipidomics OR "mass spectrometry imaging")`. Refine the query rather
     than raising `retmax`; every abstract costs context. Take compounds the
     abstracts explicitly report as associated with the topic, recording PMID,
     direction of change and sample type; skip those outside the acquired m/z
     range or unlikely to ionize in this polarity. Cite only what the abstracts
     say.
2. **Masses.** `resolve_compound_mass(names)`; on a miss retry with a synonym,
   the IUPAC name or a CID (`namespace="cid"`). Lipid sum-composition names
   ("PC 34:1") often fail: resolve one molecular species with that sum
   composition (e.g. POPC for PC 34:1) — they share formula and mass — and
   report the result at sum-composition level.
3. **Match.** `annotate_adducts(compounds=<resolved records>,
   imzml_path=<aligned output>, polarity="auto", ion_source=<from load_imzml>,
   ppm_tol=<per the reference>)`. Without `ion_source` only [M+H]⁺/[M−H]⁻ are
   tried. Read each compound's `hits`, not only `best_*`:
   - `within_tol` false → not detected at this tolerance (the nearest peak is
     always returned, however far away).
   - `within_tol` true → candidate. Check the ppm against the dataset's
     calibration offset (the evaluation's reference-ion check), whether other
     adducts of the same compound match with the same offset (supportive), and
     the look-alikes listed in the reference.
4. **Where.** `score_ion_structure` and `render_ion_images` for the matched
   candidates: present, structured and on tissue before "detected in region
   X". Compare with the cluster image if one exists.
5. **Absent?** Before reporting a compound as not detected, confirm its m/z
   lies in the acquired range and its expected adduct fits the polarity, and
   that no filter removed it: look it up in the adduct audit CSV, test the
   `_aligned` (pre-filter) file, or render the raw data's ion image at the
   theoretical m/z.

## Report

Write in the user's language.

1. **Question and data** — dataset and processing history (normalization,
   mask, filters), polarity with its source, inherited caveats.
2. **Regions** (exploratory) — the cluster image; a table of cluster, pixels
   and fraction, position and shape, technical flag, stability, top structured
   markers (m/z, enrichment, Moran's I), likely class.
3. **Compounds** (targeted) — a table of compound (source: user or PMID),
   formula, adduct, theoretical m/z, observed m/z, ppm, within tolerance,
   effective pixels, Moran's I, location, confidence level.
4. **Literature** (if used) — PMID, year, compound(s), reported association;
   one line each.
5. **Limitations and how to confirm** — annotation level and look-alikes,
   single section, normalization effects, artifacts; what would settle it
   (on-tissue MS/MS, standards, registered histology, replicate sections).
6. **Reproduce** — n_clusters, random_state and `result_npz` of each run,
   ppm_tol, ion_source, polarity, PubMed queries with date windows.

## Errors

- "UMAP/KMeans needs a shared m/z axis", or `needs_user_input` from
  `annotate_adducts` → the data need pass 1 of `msi-preprocessing`.
- "Clustering result not found" → pass `artifacts.result_npz` of the run you
  mean.
- `resolve_compound_mass` misses → synonym, IUPAC name, CID; report what stays
  unresolved instead of estimating a mass.
- PubMed or PubChem network errors → retry once, then report the gap; never
  fill it from memory.
