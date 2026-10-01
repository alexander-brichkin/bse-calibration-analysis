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
    pitch, angle = B.lattice_pitch_and_angle(obs)
    # No second pass here: the recovered points have a position but no area,
    # and this function measures area against position.
    scale, rot, trans, ii, jj, _, _ = B.register(ref_mm, obs, pitch,
                                                 angle_deg=angle)
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

    # ---------------------------------------------------------------- figure
    # A radial profile is the wrong display for a tilt: the effect is
    # directional, so averaging over angle is exactly what destroys it. The
    # two panels below show the measurement the argument actually rests on -
    # where each image's gradient points, and the trend along the direction in
    # which the test plates depart from the flat baseline.
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))

    lim = max(np.hypot(*g) for g in grads.values()) * 100 * 1.25
    # One colour per image, used in BOTH panels. Reusing a colour for a
    # different image between the two halves of one figure is how a reader
    # ends up comparing the wrong curves.
    palette = {}
    for i, n in enumerate(flats):
        palette[n] = ["#9aa8ae", "#7d8f9a", "#53656e"][i % 3]
    for i, n in enumerate(tests):
        palette[n] = ["#b03030", "#8e44ad", "#c1660d", "#0a7e92"][i % 4]

    ax[0].axhline(0, color="#dddddd", lw=.8, zorder=0)
    ax[0].axvline(0, color="#dddddd", lw=.8, zorder=0)
    for n in flats:
        v = grads[n] * 100
        ax[0].annotate("", xy=v, xytext=(0, 0),
                       arrowprops=dict(arrowstyle="->", color=palette[n],
                                       lw=2))
        ax[0].annotate(n, xy=v, fontsize=8, color=palette[n])
    bv = baseline * 100
    ax[0].annotate("", xy=bv, xytext=(0, 0),
                   arrowprops=dict(arrowstyle="->", color="#16232B", lw=2.6))
    ax[0].annotate("flat baseline", xy=bv, fontsize=9, fontweight="bold")
    circle = plt.Circle(bv, floor * 100, color="#16232B", alpha=.08)
    ax[0].add_patch(circle)
    for n in tests:
        c = palette[n]
        v = grads[n] * 100
        ax[0].annotate("", xy=v, xytext=(0, 0),
                       arrowprops=dict(arrowstyle="->", color=c, lw=2.4))
        ax[0].annotate(n, xy=v, fontsize=9, color=c, fontweight="bold")
    ax[0].set_xlim(-lim, lim); ax[0].set_ylim(-lim, lim)
    ax[0].set_aspect("equal"); ax[0].grid(alpha=.25)
    ax[0].set_xlabel("gradient along plate X, % per mm")
    ax[0].set_ylabel("gradient along plate Y, % per mm")
    ax[0].set_title("where the spot-area gradient points\n"
                    "shaded circle = flat-to-flat scatter", fontsize=11)

    for n in flats + tests:
        c = palette[n]
        plate, area, _, _ = data[n]
        excess = grads[n] - baseline
        if np.hypot(*excess) < 1e-12:
            continue
        unit = excess / np.hypot(*excess)
        t = plate @ unit                      # mm along the departure direction
        edges = np.linspace(t.min(), t.max(), 9)
        mid = 0.5 * (edges[:-1] + edges[1:])
        prof = [area[(t >= lo) & (t < hi)].mean()
                if ((t >= lo) & (t < hi)).sum() > 4 else np.nan
                for lo, hi in zip(edges[:-1], edges[1:])]
        ax[1].plot(mid, prof, "o--" if n in flats else "o-", color=c, label=n)
    ax[1].set_xlabel("position along each image's own excess-gradient direction, mm")
    ax[1].set_ylabel("spot area / image mean")
    ax[1].grid(alpha=.25); ax[1].legend(fontsize=8)
    ax[1].set_title("the directional trend the radial view hides", fontsize=11)

    fig.suptitle("Week 3 - spot area across the plate: a tilt is directional, "
                 "so it is measured along a direction", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight", metadata={"Software": None})
    print("\nwritten:", out_path)


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--flat" not in a or "--test" not in a:
        sys.exit(__doc__)
    f, t = a.index("--flat"), a.index("--test")
    main(a[f + 1:t], a[t + 1:])
