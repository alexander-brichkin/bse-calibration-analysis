# Contribution statement

MT035A *Additive Manufacturing in Metal*, HT2026 — group project
**Automatic Analysis of Calibration Patterns from BSE Images**

<!-- FILL IN before submission: the names and the two bracketed lines below
     are the only things in this repository that are not derived from the
     data. Do not submit it with the brackets still in place. -->

| Member | Contribution |
|---|---|
| **Alekszandr Bricskin** | Registration pipeline in this repository: spot detection, the closed-form similarity fit and the ICP matching, the rotation seed and the predicted-position second pass, the variance decomposition of systematic against random deviation, the Week-3 directional focus analysis, the Week-4 robustness sweep and the synthetic-ground-truth self-test. Wrote `README.md` and built the seminar deck. |
| **Felix** | An independently written pipeline covering Weeks 1–2, including the morphological top-hat segmentation, the watershed split of merged spots and the predicted-position search. Cross-checking his results against this one is what identified the 10.25 px/mm constant as the source of the disputed 2.4 % scale discrepancy; his predicted-position idea is adopted here (see `refine_by_prediction` in `bse_register.py`). |
| **Moiz** | Independent spot-detection overlays used as a visual cross-check of the detected centroids. |
| *[name]* | *[contribution]* |

## How the cross-checks were used

The three implementations were written separately and compared rather than
merged. Two disagreements were resolved that way, and both changed a result:

1. **The scale discrepancy.** Felix's pipeline reported a 2.4 % difference in
   pixel size. Re-running his own code with the assumed `BSE_MM_PER_PX = 1/10.25`
   constant disabled removed the discrepancy entirely, which located the
   problem in the assumed constant rather than in either measurement.
2. **The reference file.** The millimetre columns supplied in the Week-1
   reference CSV assume a rounded image scale of 18.000 px/mm, while the
   pattern's own lattice pitch gives 17.9555 px/mm — a 0.25 % error. This
   repository re-derives the millimetre coordinates from the pixel columns
   and does not use the supplied millimetre columns.

Methods were adopted from the other implementations only where they measured
better on our own test set, not because they were available: the watershed
split and the predicted-position second pass are in, the morphological
top-hat is not, because it finds 13 more spots on a clean image but returns a
confident 11.5 % pixel-size error at a noise level where the median filter is
still correct. `robustness.py` measures both and the comparison is in
`README.md`.

## Open questions we could not close from the supplied data

* The provenance of the 10.25 px/mm figure in the Week-2 announcement.
* Whether the plate labelled *flat 3* is in fact flat; its focus gradient sits
  at 1.5× the flat-to-flat scatter, above both flats but below the clearly
  tilted plate.
