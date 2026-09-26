# BSE calibration-pattern analysis for PBF-EB

Automatic measurement of how far a melted calibration pattern has drifted from
the pattern the machine was programmed with, from a single backscattered-
electron (BSE) image.

The electron beam in a PBF-EB machine is both the tool and the probe: the same
deflection system melts the powder and forms the BSE image. Melting a known
lattice on a plate and imaging it therefore gives a direct check on beam
positioning — but that check is normally done by eye, which is slow and depends
on who is doing it. This code does it automatically and returns numbers you can
trend.

Written for the group project in **MT035A, Additive Manufacturing in Metal,
Mid Sweden University, HT2026**. Course-mates: clone it, drop the course files
in, run `./run_all.sh`.

---

## What it measures

Given a BSE snapshot and the programmed pattern, it fits

```
[x', y'] = s · R(θ) · [x, y] + [tx, ty]
```

and reports:

| measurand | where it comes from |
|---|---|
| ΔX, ΔY | translation of the pattern origin |
| Δθ | rotation — only defined modulo 90°, because a square lattice is 4-fold symmetric |
| scale s | and with it the **derived BSE pixel size**, which is an output, not an input |
| local distortion | the per-spot residual left after the fit is removed |
| radial / tangential | components of that residual about the centre of the field |

Four degrees of freedom is deliberately few. The model cannot describe an
anisotropic scale, a shear, perspective from a tilted plate, barrel or
pincushion distortion, or a per-spot placement error — all of which therefore
stay in the residual, which is exactly where the distortion measure comes from.
A richer model would absorb the defect into the fit and make the machine look
better than it is.

---

## Quick start

```bash
git clone <this repo>
cd bse-calibration-analysis

python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python3 selftest.py
```

`selftest.py` needs **no course data**. It builds a synthetic calibration
pattern and a synthetic BSE image with a transform it already knows, then
checks the pipeline recovers it:

```
synthetic pattern: 1517 spots, 2.0 mm pitch
ground truth: scale 10.0 px/mm, angle 0.35 deg, shift (451.0, 448.0) px

  [PASS] clean                        matched 1485/1517 (97.9%) | scale err 0.001 % | angle err 0.0000 deg | shift err 0.00 px
  [PASS] low contrast (0.18)          matched 1486/1517 (98.0%) | scale err 0.004 % | angle err 0.0004 deg | shift err 0.00 px
  [PASS] noisy (sigma 0.06)           matched 1479/1517 (97.5%) | scale err 0.004 % | angle err 0.0001 deg | shift err 0.02 px
  [PASS] blurred (sigma 1.5 px)       matched 1496/1517 (98.6%) | scale err 0.004 % | angle err 0.0007 deg | shift err 0.00 px
  [PASS] strong background gradient   matched 1484/1517 (97.8%) | scale err 0.002 % | angle err 0.0004 deg | shift err 0.01 px
  [soft] low contrast + noisy         matched   29/1517 ( 1.9%) | scale err 110.496 %
```

If that prints `ALL CHECKS PASSED`, the install is good. The last row is
deliberately degraded past the point where the method should be trusted — note
how the failure announces itself: the matched fraction collapses to 2 %. That
is the reliability criterion. **Never read the deviation numbers without
reading the matched fraction first.**

---

## Running it on real data

The course files are not in this repository — they belong to the teaching team.
Put them in the repository root (they are already in `.gitignore`):

```
pattern_image.png          the programmed reference pattern, as supplied
reference_spots.csv        the Week-1 spot table extracted from it
bse-snapshot-*.png         the flat-plate BSE snapshots
```

Then:

```bash
./run_all.sh                                  # every bse-snapshot-*.png
./run_all.sh bse-snapshot-img1.png            # or the ones you name
```

`run_all.sh` processes each snapshot, merges the results, and cross-checks the
orientation of the first two. About 5 s per image.

### One image at a time

```bash
python3 bse_register.py bse-snapshot-img1.png [reference_spots.csv]
```

prints the fit and writes, next to the image:

- `…-overlay.png` — the verification figure: programmed pattern drawn on the
  snapshot, a zoom, the residual vector field, the residual histogram, the
  numbers
- `…-deviations.csv` — one row per matched spot
- `…-summary.csv` — one row with the global fit

### The other scripts

| command | what it does |
|---|---|
| `python3 analyse_repeatability.py A-deviations.csv B-deviations.csv` | splits the deviation into a systematic and a random part by comparing two snapshots |
| `python3 check_orientation.py A.png B.png` | decides whether two snapshots are in the same orientation (the lattice alone cannot tell) |
| `python3 make_reference_figure.py` | redraws the Week-1 reference figure |
| `python3 make_workflow_figure.py bse-snapshot-img2.png` | original → processed → identified → overlay, on one region |
| `python3 make_slide_overlay.py bse-snapshot-img2.png` | a wide three-panel overlay laid out for a 16:9 slide |

---

## Output columns

`…-deviations.csv`, one row per matched spot:

| column | meaning |
|---|---|
| `ref_idx`, `obs_idx` | row in the reference table / index of the detected spot |
| `ref_x_mm`, `ref_y_mm` | programmed position, mm, relative to the reference spot |
| `pred_x_px`, `pred_y_px` | where the fit says that spot should appear |
| `obs_x_px`, `obs_y_px` | where it actually is |
| `obs_x_mm`, `obs_y_mm` | observed centre in mm **in the BSE image's own calibration**, origin on the BSE reference spot |
| `dx_um`, `dy_um` | the local deviation vector |
| `dist_um` | its magnitude |
| `radius_mm` | distance from the centre of the field |
| `dev_radial_um` | component pointing away from the centre (+ outward) |
| `dev_tangential_um` | component perpendicular to it (+ counter-clockwise) |

`…-summary.csv`, one row per image: `spots_detected`, `matched_pairs`,
`reference_points`, `lattice_pitch_px`, `scale_px_per_mm`, `pixel_size_um`,
`rotation_deg_mod90`, `dx_px`, `dy_px`, `rms_um`, `median_um`, `p95_um`,
`max_um`.

---

## How it works

1. **Isolate the plate** — threshold, fill holes, keep the largest component,
   erode by 12 px so the bright rim is excluded.
2. **Flatten the background** — BSE grey level follows atomic number, local
   tilt and topography, so the background drifts across the field and one
   global threshold fails. A 25 px median filter estimates it; subtracting
   leaves a flat map where Otsu behaves.
3. **Segment** — Otsu on the flattened map, drop blobs under 8 px, fill holes.
4. **Filter by shape** — keep blobs with an area 0.3–4× the median and
   eccentricity below 0.85. Scratches, rim fragments and merged pairs go.
5. **Centroid** — the spot centre is the centroid of the segmented area.
6. **Seed** — measure the lattice pitch from nearest-neighbour distances and
   use pitch / 2.000 mm as the starting scale.
7. **Register** — Umeyama's closed-form least-squares similarity fit, iterated
   ICP-style 30 times. Matching is mutually-nearest-neighbour inside a gate of
   0.45 × pitch, i.e. below half a lattice step, so a spot cannot lock onto its
   neighbour.
8. **Measure** — ΔX, ΔY, Δθ and s come from the fit; what is left over is the
   local distortion.

### Two things worth knowing before you trust a number

**The millimetre columns of a reference table may be wrong.** In our course
file they were written assuming a round 18.000 px/mm, while the lattice
actually sits at 35.911 px (17.9555 px/mm) — a 0.25 % stretch that lands
straight in the derived pixel size. `load_reference()` therefore ignores those
columns and re-derives millimetres from the pixel columns against the nominal
`PITCH_MM`. Two checks confirmed this: the field then measures exactly
88.000 mm (44 steps × 2.000 mm), and the derived BSE pixel size moves to within
0.03 % of a round 100 µm/px instead of 0.28 %.

**The rotation is only defined modulo 90°.** A square lattice is 4-fold
symmetric, so the fit returns −90.01° and −0.01° with identical matches and
identical residuals. Use `check_orientation.py`, or key the orientation off a
marked reference spot.

---

## What we measured with it

Two flat-plate BSE snapshots, 897 × 899 px:

| | img1 | img2 |
|---|---|---|
| spots detected / matched | 1282 / 1272 | 1441 / 1441 |
| derived pixel size | 99.969 µm | 99.977 µm |
| Δθ (mod 90°) | −0.0119° | −0.0090° |
| residual RMS | 51.0 µm | 47.8 µm |
| median / 95th / max | 36 / 93 / 353 µm | 37 / 86 / 264 µm |

Two independent images agreeing on the pixel size to 82 ppm, and both landing
0.03 % from a round 100 µm/px, is the main reason to believe the fit — nothing
in the pipeline takes a pixel size as input.

Comparing the two images spot by spot splits the deviation into a part that
repeats (systematic, 32.5 µm) and a part that does not (random, 35.5 µm). About
16 µm per axis of the random part is pixel quantisation in the supplied
reference file, not the machine. So a single image cannot resolve a drift below
roughly 35 µm; averaging N images pulls that floor down as 1/√N while the
systematic part stays put.

---

## Limitations

- Validated on two real snapshots and a synthetic set. Not a qualified
  measurement procedure.
- The derived pixel size rests on the nominal lattice pitch being true. If the
  programmed pitch is wrong, the scale is wrong by the same factor.
- Spots near the rim are measured less reliably — contrast falls off and some
  are clipped by the erosion margin. Part of the growth of the residual towards
  the edge is probably this, not the machine.
- Registering a calibration pattern successfully does **not** qualify a PBF-EB
  process. It says nothing about melt-pool behaviour, porosity, layer bonding,
  the Z axis or powder spreading. It is a single-layer, in-plane geometric
  check and needs a physically measured artefact alongside it.

---

## Method sources

- Umeyama, S. (1991). Least-squares estimation of transformation parameters
  between two point patterns. *IEEE TPAMI* 13(4), 376–380.
- Besl, P. & McKay, N. (1992). A method for registration of 3-D shapes.
  *IEEE TPAMI* 14(2), 239–256.
- Otsu, N. (1979). A threshold selection method from gray-level histograms.
  *IEEE Trans. SMC* 9(1), 62–66.
- van der Walt, S. et al. (2014). scikit-image: image processing in Python.
  *PeerJ* 2:e453.
- ISO/ASTM 52930 (installation, operation and performance qualification) and
  ISO/ASTM 52920 (requirements for industrial additive manufacturing
  processes).

## License

MIT — see [LICENSE](LICENSE). The course data is not covered by it and is not
included here.
