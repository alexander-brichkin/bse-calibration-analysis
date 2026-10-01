#!/usr/bin/env python3
"""
MT035A Week 3 - tilted plate: what the data supports, and what it does not.

Two things make this harder than it looks.

First, foreshortening is not a tilt signature in this system. In a projection
a tilted plane is compressed along the tilt axis, so one expects an anisotropic
scale. Here the beam is deflected to a commanded (x, y) and the BSE image is
formed by the same deflection, so both the writing and the reading use the same
in-plane coordinates: the tilt changes the working distance, not the geometry.
The measurement bears this out - the tilted plates show no more anisotropy than
the flat ones. Any tilt estimate therefore has to come from focus, not shape.

Second, a flat plate already shows a strong directional variation in spot area,
because the BSE detector sits to one side and the collected signal varies across
the field. That shading is of the same order as the tilt effect, so a raw
gradient cannot separate the two. The flats are used here as a baseline and
subtracted; what remains is the part that is not shading, and the scatter
between the flats is the noise floor against which it has to be judged.

Usage
-----
    python3 tilt_analysis.py --flat FLAT1.png FLAT2.png [...] \
                             --test TILT1.png [TILT2.png ...]

Requires: numpy, scipy, scikit-image, matplotlib.
"""

import os
import sys

import numpy as np
from skimage import io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bse_register as B

# Relaxed on purpose: the defaults reject smeared and merged blobs, which is
# right for positions and wrong for measuring how much the spots smear.
RELAXED = dict(area_lo=0.15, area_hi=10.0, max_ecc=0.98)


def spot_field(path, ref_mm):
    """Normalised spot area as a function of position on the plate, in mm."""
    spots = B.detect(io.imread(path), **RELAXED)
    obs, area, ecc = spots[:, :2], spots[:, 3], spots[:, 4]
    pitch, _ = B.lattice_pitch_and_angle(obs)
    scale, rot, trans, ii, jj, _, _ = B.register(ref_mm, obs, pitch)
    plate = ((rot.T @ (obs - trans).T).T / scale) * np.array([1, -1])
    return plate, area / area.mean(), ecc, scale


def gradient(plate, value):
    """Least-squares plane through value(x, y): returns d/dx, d/dy per mm."""
    design = np.column_stack([np.ones(len(plate)), plate[:, 0], plate[:, 1]])
    coef = np.linalg.lstsq(design, value, rcond=None)[0]
    return np.array([coef[1], coef[2]])


def main(flat_paths, test_paths, out_path="week3-tilt-analysis.png"):
    ref_mm, _, _ = B.load_reference("reference_spots.csv")
    data = {}
    for p in flat_paths + test_paths:
        name = os.path.splitext(os.path.basename(p))[0]
        data[name] = spot_field(p, ref_mm)

    flats = [os.path.splitext(os.path.basename(p))[0] for p in flat_paths]
    tests = [os.path.splitext(os.path.basename(p))[0] for p in test_paths]

    grads = {n: gradient(d[0], d[1]) for n, d in data.items()}
    baseline = np.mean([grads[n] for n in flats], axis=0)
    # how much the flats disagree with each other: the floor for any claim
    floor = np.mean([np.hypot(*(grads[n] - baseline)) for n in flats])

    def show(v):
        return f"{np.hypot(*v) * 100:6.3f} %/mm at {np.degrees(np.arctan2(v[1], v[0])) % 360:5.1f} deg"

    print("Directional gradient of spot area, in plate coordinates\n")
    for n in flats + tests:
        print(f"  {n:10s} raw    {show(grads[n])}")
    print(f"\n  flat baseline {show(baseline)}")
    print(f"  flat-to-flat scatter (the noise floor): {floor * 100:.3f} %/mm\n")
    print("Excess over the flat baseline - the part that is not detector shading:\n")
    verdict = {}
    for n in tests + flats:
        excess = grads[n] - baseline
        ratio = np.hypot(*excess) / max(floor, 1e-9)
        verdict[n] = ratio
        print(f"  {n:10s} excess {show(excess)}   {ratio:4.1f}x the floor")

    print("\nNo tilt ANGLE is reported. Converting a spot-area gradient into "
          "degrees\nneeds the beam's depth-of-focus characteristic - how spot "
          "area grows per mm\nof working-distance error - which is not in the "
          "supplied data. What the data\nsupports is a focus gradient across "
          "the plate, consistent with a tilt, with the\nmagnitudes above.")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    for n in flats + tests:
        plate, area, ecc, _ = data[n]
        r = np.hypot(*plate.T); r /= r.max()
        edges = np.linspace(0, 1, 7)
        mid = 0.5 * (edges[:-1] + edges[1:])
        style = "o--" if n in flats else "o-"
        ax[0].plot(mid, [area[(r >= lo) & (r < hi)].mean()
                         for lo, hi in zip(edges[:-1], edges[1:])], style, label=n)
        ax[1].plot(mid, [ecc[(r >= lo) & (r < hi)].mean()
                         for lo, hi in zip(edges[:-1], edges[1:])], style, label=n)
    ax[0].set_ylabel("spot area / image mean"); ax[0].set_title("spot area", fontsize=11)
    ax[1].set_ylabel("eccentricity"); ax[1].set_title("spot elongation", fontsize=11)
    for a in ax:
        a.set_xlabel("distance from plate centre, fraction of radius")
        a.grid(alpha=.25); a.legend(fontsize=8)
    fig.suptitle("Week 3 - spot geometry across the plate "
                 "(dashed = flat reference plates)", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight", metadata={"Software": None})
    print("\nwritten:", out_path)


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--flat" not in a or "--test" not in a:
        sys.exit(__doc__)
    f, t = a.index("--flat"), a.index("--test")
    main(a[f + 1:t], a[t + 1:])
