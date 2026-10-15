# Annotation evidence for MSI

How far an m/z match can be trusted, and how to say so. Masses below are
monoisotopic and were computed from elemental masses; tolerances are review
settings, not instrument specifications.

## Confidence levels

Adapted from Schymanski et al., *Environ. Sci. Technol.* 2014 (levels for
high-resolution MS). Imaging without MS/MS stays at the lower levels.

| Level | Evidence needed | With mzorbit alone | Write it as |
|---|---|---|---|
| 1 Confirmed structure | Reference standard measured the same way (m/z + MS/MS, ideally on tissue) | No | "identified as …" |
| 2 Probable structure | MS/MS library match or diagnostic fragments | No | "probably …" |
| 3 Tentative candidates | Independent evidence for a structure class (e.g. diagnostic fragment); several structures possible | Only with evidence from outside mzorbit | "tentatively a PI" |
| 4 Molecular formula | Accurate mass within tolerance, consistent with the calibration offset, no competing formula in tolerance; isotope peak supportive | Best case | "C47H83O13P as [M−H]⁻ (−1.8 ppm); candidate PI 38:4" |
| 5 Exact mass of interest | An m/z with a pattern worth reporting | Always | "m/z 885.550" |

Lipid names follow the LIPID MAPS shorthand (Liebisch et al., *J. Lipid Res.*
2020). Accurate mass reaches the species level only — the sum composition
"PC 34:1" — never the molecular species "PC 16:0_18:1" (needs MS/MS) or sn
positions "PC 16:0/18:1".

## Matching tolerance

- Set `ppm_tol` to about 2–3× the analyzer's typical single-peak error: TOF
  ~50, Orbitrap ~10, FT-ICR ~5, unknown ~25 ppm (mzorbit's typical errors are
  20 / 5 / 2 / 10 ppm).
- Judge each candidate by its distance from the dataset's systematic offset
  (the evaluation's reference-ion check), |best_ppm − offset|, not by
  |best_ppm| alone.
- **Re-binned profile data** know positions only to ± half a bin. For 0.1225 Da
  bins that is ±306 ppm at m/z 200, ±122 ppm at 500 and ±77 ppm at 800. Set
  `ppm_tol` to half a bin at the lowest m/z of interest, check each match's
  error in Da (≤ half a bin), and expect several formulas to fit, above all at
  low m/z.
- **Unaligned centroid data**: the tool pools the peaks of ~100 spectra, so
  some peak sits near almost any m/z. Match against the aligned output.
- State the tolerance with every result: below the calibration scatter it
  hides real ions, above it it admits look-alikes.

## Adducts mzorbit tries

| Polarity | Always | MALDI / ESI / DESI add | ESI / DESI add |
|---|---|---|---|
| positive | [M+H]⁺ | [M+Na]⁺, [M+K]⁺ | [M+NH₄]⁺, [M+H−H₂O]⁺ |
| negative | [M−H]⁻ | [M+Cl]⁻, [M−2H+Na]⁻, [M−2H+K]⁻ | [M+HCOO]⁻, [M+CH₃COO]⁻ |

Not tried — a compound seen only in these forms looks absent:

- intrinsic cations [M]⁺ (heme b, m/z 616.1767);
- in-source fragments and neutral losses: the phosphocholine head group
  (m/z 184.0733) from PC and SM; loss of trimethylamine (59.0735) from
  [PC+Na]⁺ / [PC+K]⁺; water losses;
- dimers and multiply charged ions.

To check one of these forms, derive its m/z in code from the tool's
`theoretical_mz` (e.g. [M+Na]⁺ − 59.0735), say that it was derived, and look
at its ion image.

## Look-alikes

| Look-alike | Example | Separation | What to do |
|---|---|---|---|
| Same formula, different class | PC 34:1 and PE 37:1 are both C42H82NO8P (759.5778) | None by mass | Decide by polarity behaviour (PE gives a strong [M−H]⁻, PC none) and biology (odd-chain PE is rare); otherwise report both |
| Na⁺ vs H⁺ adducts | [PC 32:0+Na]⁺ 756.5514 vs [PC 34:3+H]⁺ 756.5538 | 2.4 mDa (3.2 ppm); needs resolving power ≈ 315,000 | Unresolved on TOF and most imaging Orbitrap runs: report both |
| ¹³C₂ isotope of a neighbour | ¹³C₂ peak of [PC 34:2+H]⁺ 760.5761 vs [PC 34:1+H]⁺ 760.5851 | 8.9 mDa (11.8 ppm) | Inside TOF tolerance: part of an apparent PC 34:1 signal belongs to PC 34:2 |
| Wrong polarity | [M+H]⁺ lies 7.84 mDa above the ¹³C₂ isotope of [M−H]⁻ (8.8 ppm for PI 38:4) | Small | With `polarity="both"`, a near-miss in the wrong polarity is expected |
| Isomers in one sum composition | PC 34:1 covers PC 16:0_18:1, PC 18:0_16:1, PC 14:0_20:1 …, each with several double-bond positions | None | Report the sum composition |
| Xenobiotic vs endogenous | PFOS [M−H]⁻ 498.9302 (mass defect −0.07) vs lipids near m/z 500, e.g. LPE 18:1 [M−H]⁻ 478.2939, LPC 16:0 [M+H]⁺ 496.3398 (defects +0.29, +0.34) | Mass defect | A fluorinated candidate matched to a peak with a lipid-like defect is a mismatch; exogenous compounds also need a tighter tolerance and control tissue |

## Contaminants and non-biological ions

- Matrix ions and clusters (matrix table in
  `../../msi-dataset-evaluation/references/msi-domain-knowledge.md`): bright
  across the whole raster, bare slide included; clusters reach m/z 500–600.
- Polymers: polyethylene glycol (series spaced 44.0262 Da), polysiloxanes
  (74.0188 Da).
- Plasticizers: DEHP [M+H]⁺ 391.2843, [M+Na]⁺ 413.2662.
- FA 16:0 and FA 18:0 also come from handling and solvents.
- An ion as intense off tissue as on it is not tissue chemistry.

## Reading structure scores

`score_ion_structure` integrates each m/z over the dataset's tolerance:

- `effective_pixels`: pixels above `detection_quantile` × the 99th-percentile
  intensity (default: half the near-maximum) — where the ion is strong, not
  everywhere it is present.
- `tissue_pixels` and `occupancy_fraction` count all acquired pixels,
  background included unless the data were masked.
- `spatial_autocorr`: Moran's I between 4-neighbours. Near 0 behaves like
  noise; roughly ≥ 0.3 is coherent structure — but striping and edge rings
  score high too, so look at the image.
- `dynamic_range`: maximum over median of the positive pixels; very high with
  few effective pixels means hot spots (dust, salt, matrix crystals).
- Results are ranked by `effective_pixels`, which favours large, diffuse ions;
  a small, specific structure ranks low.

## Clusters and markers

- `cluster_marker_ions` enrichment = cluster mean / mean of all cluster means.
  Ions found in one cluster only saturate at `n_clusters`, background clusters
  shift every ratio, and there is no intensity floor: confirm with structure
  scores and images.
- KMeans returns k clusters whether or not the tissue has k compartments;
  stability across k and seeds is the evidence that a region is real.
- Unnormalized data can cluster by total signal (section thickness, matrix
  density) rather than composition; in TIC-normalized data a few very intense
  ions can imprint their pattern on all others (normalization section of
  `../../msi-preprocessing/references/step-guide.md`).

## Literature evidence

- Abstracts are summaries: cite what they state, with the PMID, and do not
  infer methods, tissues or effect sizes they do not report.
- A literature association says a compound matters for the topic, not that it
  is present here; a match here says an ion at that m/z is present, not that it
  is that compound.
- Prefer MS-based studies (lipidomics, metabolomics, MS imaging) of the same
  organism and tissue; serum or cell-line findings transfer poorly to sections.
- The default search window is the last ten years; widen it with `mindate` for
  long-established compounds.

## Mass ranges and typical classes (MALDI; overlapping, not diagnostic)

| m/z | Negative mode | Positive mode |
|---|---|---|
| < 200 | Small acids and amino acids (taurine, glutamate), matrix | Small amines (choline, carnitine), matrix |
| 200–400 | Free fatty acids, glutathione, first nucleotides | Acylcarnitines, matrix clusters |
| 400–600 | Nucleotides, lysophospholipids | Lysophospholipids (LPC, LPE) |
| 600–950 | Glycerophospholipids (PA, PE, PG, PS, PI), sulfatides | Heme b (616), PC, SM and PE as H⁺/Na⁺/K⁺ adducts, TG as Na⁺/K⁺ adducts toward 950 |
| > 1000 | Gangliosides, cardiolipins (~1400–1500) | Gangliosides; peptides with a suitable matrix |
