# MSI domain knowledge for dataset evaluation

Background for interpreting what the mzorbit tools measure. Values are typical
ranges for review, not specifications of any one instrument.

## imzML in brief

- Two files with one base name: `.imzML` (XML header: metadata as controlled
  vocabulary terms, plus per-spectrum offsets) and `.ibd` (binary arrays). A
  UUID in both must match; renaming one file breaks the pair.
- **Storage**: *continuous* (`IMS:1000030`) — one m/z array shared by all
  spectra; *processed* (`IMS:1000031`) — every spectrum stores its own m/z array.
- **Representation**: *profile* (`MS:1000128`) — sampled peak shapes; *centroid*
  (`MS:1000127`) — one point per peak. Aligned peak lists are written as
  continuous centroid data (one shared axis of peak positions).
- **Polarity**: `MS:1000130` positive scan, `MS:1000129` negative scan, usually in
  every spectrum's cvParams. MassFlow's writer drops it; mzorbit keeps it in the
  sidecar (`<base>.metadata.json`) so processed files inherit it.
- **Geometry**: `max count of pixels x/y`, pixel size (`IMS:1000046`/`IMS:1000047`),
  scan pattern (meandering vs flyback), scan type (horizontal vs vertical lines)
  and line scan direction. Coordinates are 1-based by convention.
- Vendor exports are often **re-binned**: a "continuous profile" file whose
  uniform bin grid (e.g. 0.1225 Da) is coarser than the true peak width. Each ion
  then occupies about one bin; mass accuracy is limited to roughly ± half a bin.

## Mass analyzers

| Class | Resolving power (FWHM) | Mass accuracy | Imaging notes |
|---|---|---|---|
| TOF (reflectron) | 10k–40k | 5–50 ppm | Fast pixels; m/z shifts with sample height/topography and temperature; calibration drifts over long runs |
| TOF (linear) | 1k–5k | 100+ ppm | Proteins/large peptides |
| Q-TOF | 20k–60k | 2–10 ppm | Common for DESI |
| Orbitrap | 60k–480k at m/z 200 (falls as √m/z) | 1–5 ppm | Slower pixels; resolution set per run |
| FT-ICR | 100k–1M+ | < 1–2 ppm | Long transients, very large files |
| Trapped ion mobility TOF | as TOF + mobility | 2–10 ppm | imzML exports often collapse the mobility dimension |

## Ion sources

| Source | Typical pixel size | Ions | Notes |
|---|---|---|---|
| MALDI | 5–100 µm (≤ 5 µm specialised) | Almost only singly charged | Matrix clusters dominate below m/z ~500–600; MALDI-2 post-ionization boosts lipids |
| DESI | 30–200 µm | ESI-like; multiply charged for peptides/proteins | No matrix; solvent and spray geometry drive sensitivity; Cl⁻ adducts with chlorinated solvents |
| IR-MALDESI | 50–150 µm | ESI-like | Ice matrix; multiply charged ions possible |
| SIMS (ToF-SIMS) | < 1 µm | Fragments, low m/z | Elemental/fragment imaging, limited intact lipids |
| LAESI / nano-DESI | 50–400 µm / 10–50 µm | ESI-like | Ambient |

## Matrices and the low-mass region (MALDI)

| Matrix | Polarity | Characteristic matrix ions (m/z) |
|---|---|---|
| DHB (2,5-dihydroxybenzoic acid) | positive | 137.02 [M+H−H₂O]⁺, 155.03 [M+H]⁺, 177.02 [M+Na]⁺, 273.04 [2M+H−2H₂O]⁺ |
| CHCA (α-cyano-4-hydroxycinnamic acid) | positive | 172.04 [M+H−H₂O]⁺, 190.05 [M+H]⁺, 212.03 [M+Na]⁺, 379.09 [2M+H]⁺ |
| 9-AA (9-aminoacridine) | negative | 193.08 [M−H]⁻ |
| DAN (1,5-diaminonaphthalene) | both | 157.08 [M−H]⁻, 159.09 [M+H]⁺ |
| Norharmane | both | 167.06 [M−H]⁻, 169.08 [M+H]⁺ |

Matrix ions are strong on and off tissue; an ion image that is bright across the
whole raster (including bare slide) is usually matrix or contamination, not
tissue chemistry. DHB forms large crystals (speckled images); sprayed CHCA/DAN
coatings are finer.

## Polarity and lipid classes

- **Positive mode**: phosphatidylcholines (PC) and sphingomyelins (SM) dominate
  as [M+H]⁺, [M+Na]⁺, [M+K]⁺; neutral lipids (DG/TG) appear as Na⁺/K⁺ adducts;
  heme b (m/z 616.18, M⁺) marks blood; m/z 184.07 is the phosphocholine head
  group (in-source fragment).
- **Negative mode**: PI, PE, PS, PG, PA, sulfatides (ST: brain white matter,
  kidney), cardiolipins, gangliosides, free fatty acids, and small metabolites
  (taurine, glutamate, nucleotides) as [M−H]⁻; Cl⁻ adducts of neutral lipids.
- **Alkali adducts vary across tissue**: [M+K]⁺ often exceeds [M+Na]⁺ because
  cells are K⁺-rich; regions with blood, damage or ischemia shift toward Na⁺.
  Adduct ratio maps can therefore look like structure without any change in the
  lipid itself.

## Reference-ion panel (calibration and polarity check)

Neutral monoisotopic masses for `annotate_adducts` (pass `ion_source="MALDI"` to
include Na/K/Cl adducts), with the main ion expected:

| Compound | Neutral mass | Main ion (m/z) | Where it is common |
|---|---|---|---|
| PC 34:1 | 759.5778 | 760.5851 [M+H]⁺; 782.5670 [M+Na]⁺; 798.5410 [M+K]⁺ | Most tissues (positive) |
| PC 32:0 | 733.5622 | 734.5694 [M+H]⁺ | Lung, brain grey matter |
| PC 36:1 | 787.6091 | 788.6164 [M+H]⁺ | Most tissues |
| SM 34:1;O2 | 702.5676 | 703.5749 [M+H]⁺ | Kidney, brain |
| PI 38:4 | 886.5571 | 885.5498 [M−H]⁻ | Most tissues (negative) |
| PE 38:4 | 767.5465 | 766.5392 [M−H]⁻ | Brain, kidney |
| PS 40:6 | 835.5363 | 834.5291 [M−H]⁻ | Brain |
| FA 16:0 | 256.2402 | 255.2330 [M−H]⁻ | Ubiquitous (also contamination) |
| FA 18:1 | 282.2559 | 281.2486 [M−H]⁻ | Ubiquitous |
| FA 20:4 | 304.2402 | 303.2330 [M−H]⁻ | Ubiquitous |
| FA 22:6 | 328.2402 | 327.2330 [M−H]⁻ | Brain, retina |
| Taurine | 125.0147 | 124.0074 [M−H]⁻ | Kidney, liver, heart |
| Glutamate | 147.0532 | 146.0459 [M−H]⁻ | Brain |

Heme b is detected as the intact cation (m/z 616.1767) rather than a protonated
adduct: check it with `render_ion_images`, not `annotate_adducts`.

Reading the check: `annotate_adducts` always returns the nearest peak, and
`within_tol` only compares it with the `ppm_tol` you passed, so read the
`best_ppm` values together. Several panel ions with one consistent ppm offset →
polarity and calibration confirmed (report the offset). A large consistent
offset → calibration drift. Values scattered at random ppm → no real match:
wrong polarity, wrong m/z range, or a sample whose chemistry differs (e.g. a
drug or plant study): say which explanation the evidence favours. On re-binned
data, offsets up to about half a bin are expected.

Two traps:

- **Dense peak pools.** On unaligned (per-spectrum) centroid data the tool pools
  the peaks of ~100 spectra, so some peak sits near almost any m/z. Trust only an
  offset shared by several panel ions, and prefer the aligned output (one
  feature per ion) when it exists.
- **Wrong polarity can still look close.** With `polarity="both"`, the best
  adduct is simply the smallest |ppm|. In negative-mode data a compound's
  [M+H]⁺ falls 7.84 mDa above the ¹³C₂ isotope of its [M−H]⁻ (8.8 ppm for
  PI 38:4), so a near-miss in the wrong polarity is expected. A polarity is
  confirmed when it gives the tighter, shared offset for most panel ions.

## Raster patterns and spatial artifacts

| Artifact | Looks like | Likely causes |
|---|---|---|
| Row/column striping | Line-to-line intensity jumps along the scan | Laser-energy drift, raster timing; in **meandering** rasters, alternate lines differ (period-2 stripes) |
| Edge effect | Bright or dark rim along the section edge | Matrix accumulation ("coffee ring") at edges, different ablation on section borders, delocalization |
| Off-tissue margin | Bright frame around the section | Bare-slide matrix ionizes without tissue suppression; its spectra are matrix ions with little high-mass signal — tell it from an edge effect by the spectra, not the TIC |
| Speckle | Pixel-scale bright/dark noise | Large matrix crystals (DHB), low ion counts |
| Delocalization | Signal smeared beyond the tissue outline | Wet matrix application, too much solvent |
| Hot pixels | Isolated extreme pixels | Dust, salt crystals, arcing, detector saturation |
| Dropouts | Near-zero pixels inside the raster | Missed shots, holes/tears, bare slide |
| Gradual gradient | Intensity ramp across the section | Section thickness, matrix thickness, long-run sensitivity decay (correlates with acquisition order) |
| m/z drift | Peak position shifts across the section (TOF) | Topography, sample charging, temperature; alignment absorbs small drifts |

## Sample-level context to ask for when missing

Organism, tissue/organ, condition or treatment group, polarity (if not in the
file), matrix and application method, and the intended analysis (exploratory
segmentation, targeted ion, quantification). The intended analysis decides
which caveats matter.
