#!/usr/bin/env python3
"""
Spot geometry across the plate - the measurement Week 3 needs, and the control
that keeps it honest.

A tilted plate changes the working distance along one direction, so it
defocuses directionally and is foreshortened along the tilt axis. Pure
defocus, by contrast, varies with distance from the centre of the field and
leaves the geometry isotropic. The two are easy to confuse if spot size alone
is used, which is why this script profiles area, equivalent diameter and
eccentricity against radius for several images at once: a flat plate already
shows a sizeable area variation across the field from detector shading, so the
tilted case has to be read against that baseline, not against zero.

The segmentation filters are deliberately relaxed here. The defaults in
bse_register.py reject merged and smeared blobs, which is right for measuring
positions and wrong for measuring how much the spots smear.

Usage
-----
    python3 focus_profile.py IMAGE [IMAGE ...] [-o output.png]

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

RELAXED = dict(area_lo=0.15, area_hi=10.0, max_ecc=0.98)
BINS = np.array([0.0, 0.3, 0.5, 0.7, 0.85, 1.0])


def profile(path):
    spots = B.detect(io.imread(path), **RELAXED)
    x, y, diam, area, ecc = spots.T
    r = np.hypot(x - x.mean(), y - y.mean())
    r = r / r.max()
    rows = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        m = (r >= lo) & (r < hi)
        rows.append((area[m].mean(), diam[m].mean(), ecc[m].mean(), m.sum()))
    return np.array(rows), len(spots)


def main(paths, out_path="week3-focus-profile.png"):
    mid = 0.5 * (BINS[:-1] + BINS[1:])
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    print(f"{'image':24s} {'n':>6s} {'area outer/centre':>18s} "
          f"{'ecc centre -> outer':>22s}")
    for p in paths:
        rows, n = profile(p)
        name = os.path.splitext(os.path.basename(p))[0]
        ratio = rows[-1, 0] / rows[0, 0]
        print(f"{name:24s} {n:6d} {ratio:18.2f} "
              f"{rows[0, 2]:11.3f} -> {rows[-1, 2]:.3f}")
        ax[0].plot(mid, rows[:, 0], "o-", label=name)
        ax[1].plot(mid, rows[:, 1], "o-", label=name)
        ax[2].plot(mid, rows[:, 2], "o-", label=name)
    for a, t, yl in zip(ax,
                        ["spot area", "equivalent diameter", "eccentricity"],
                        ["area, px", "diameter, px", "eccentricity"]):
        a.set_xlabel("distance from field centre, fraction of plate radius")
        a.set_ylabel(yl); a.set_title(t, fontsize=11)
        a.grid(alpha=.25); a.legend(fontsize=8)
    fig.suptitle("Spot geometry across the plate - directional change points to "
                 "tilt, radial change to defocus", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight",
                metadata={"Software": None})
    print("\nwritten:", out_path)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "-o"]
    if "-o" in sys.argv:
        out = sys.argv[sys.argv.index("-o") + 1]
        main([a for a in args if a != out], out)
    else:
        main(args)
