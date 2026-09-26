#!/usr/bin/env python3
"""
MT035A Week 2 - is the measured deviation systematic or random?

The brief asks whether the observed deviations appear systematic or random and
whether different regions of the field behave differently. A single image
cannot answer that: any residual left after the similarity fit looks like a
pattern. Two independent snapshots of the same plate can, because a deviation
caused by the machine repeats from image to image while detection noise does
not.

For every spot detected in both images this script compares the two deviation
vectors and splits the variance:

    var(observed)  =  var(systematic)  +  var(noise)
    var(d1 - d2)   =  2 * var(noise)                  (noise independent)
    cov(d1, d2)    =  var(systematic)                 (systematic shared)

Usage
-----
    python3 analyse_repeatability.py A-deviations.csv B-deviations.csv \
            [output.png]

Requires: numpy, matplotlib.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    t = np.genfromtxt(path, delimiter=",", names=True)
    return {int(k): n for n, k in enumerate(t["ref_idx"])}, t


def main(path_a, path_b, out_path=None):
    out_path = out_path or "week2-repeatability.png"
    idx_a, a = load(path_a)
    idx_b, b = load(path_b)

    pairs = [(idx_a[int(k)], n) for n, k in enumerate(b["ref_idx"])
             if int(k) in idx_a]
    p = np.array([i for i, _ in pairs])
    q = np.array([j for _, j in pairs])
    name_a = os.path.basename(path_a).replace("-deviations.csv", "")
    name_b = os.path.basename(path_b).replace("-deviations.csv", "")
    print(f"spots detected in both images: {len(p)}")

    stats = {}
    for axis in ("dx_um", "dy_um"):
        x, y = a[axis][p], b[axis][q]
        corr = np.corrcoef(x, y)[0, 1]
        cov = corr * x.std() * y.std()                  # systematic variance
        noise_var = (x - y).var() / 2.0                 # random variance
        stats[axis] = dict(corr=corr,
                           systematic=np.sqrt(max(cov, 0.0)),
                           noise=np.sqrt(noise_var),
                           x=x, y=y)
        print(f"  {axis}: corr {corr:+.3f} | systematic sd "
              f"{stats[axis]['systematic']:.1f} um | random sd "
              f"{stats[axis]['noise']:.1f} um")

    d1 = np.hypot(a["dx_um"][p], a["dy_um"][p])
    d2 = np.hypot(b["dx_um"][q], b["dy_um"][q])
    diff = np.hypot(a["dx_um"][p] - b["dx_um"][q],
                    a["dy_um"][p] - b["dy_um"][q])
    rms1, rms2 = np.sqrt((d1 ** 2).mean()), np.sqrt((d2 ** 2).mean())
    rms_diff = np.sqrt((diff ** 2).mean())
    print(f"  RMS |deviation|: {rms1:.1f} and {rms2:.1f} um")
    print(f"  RMS |d1 - d2|:   {rms_diff:.1f} um "
          f"(would be {np.sqrt(2) * rms1:.1f} um if fully independent)")

    sys_tot = np.hypot(stats["dx_um"]["systematic"], stats["dy_um"]["systematic"])
    noi_tot = np.hypot(stats["dx_um"]["noise"], stats["dy_um"]["noise"])
    print(f"  => systematic {sys_tot:.1f} um vs random {noi_tot:.1f} um")

    fig, ax = plt.subplots(1, 3, figsize=(16, 5))

    s = stats["dx_um"]
    lim = max(abs(s["x"]).max(), abs(s["y"]).max()) * 1.05
    ax[0].plot([-lim, lim], [-lim, lim], color="#bbbbbb", lw=1, zorder=0)
    ax[0].axhline(0, color="#dddddd", lw=.8, zorder=0)
    ax[0].axvline(0, color="#dddddd", lw=.8, zorder=0)
    ax[0].scatter(s["x"], s["y"], s=6, alpha=.35, color="#1f77b4",
                  edgecolors="none")
    ax[0].set_xlabel(f"dX in {name_a}, um")
    ax[0].set_ylabel(f"dX in {name_b}, um")
    ax[0].set_aspect("equal")
    ax[0].set_title(f"same spot, two images\ncorrelation {s['corr']:+.2f} "
                    "-> a repeating component exists", fontsize=10)

    # RMS deviation against distance from the field centre
    for tab, pick, name, colour in ((a, p, name_a, "#1f77b4"),
                                    (b, q, name_b, "#d62728")):
        r = tab["radius_mm"][pick]
        dev = np.hypot(tab["dx_um"][pick], tab["dy_um"][pick])
        edges = np.linspace(0, r.max(), 9)
        mid = 0.5 * (edges[:-1] + edges[1:])
        rms = [np.sqrt((dev[(r >= lo) & (r < hi)] ** 2).mean())
               if ((r >= lo) & (r < hi)).sum() > 4 else np.nan
               for lo, hi in zip(edges[:-1], edges[1:])]
        ax[1].plot(mid, rms, "o-", color=colour, label=name)
    ax[1].set_xlabel("distance from field centre, mm")
    ax[1].set_ylabel("RMS |deviation|, um")
    ax[1].legend(fontsize=8)
    ax[1].set_title("the error grows towards the rim\n"
                    "in both images", fontsize=10)
    ax[1].grid(alpha=.25)

    # the shared (systematic) part of the field: the mean of the two images
    mx = 0.5 * (a["dx_um"][p] + b["dx_um"][q])
    my = 0.5 * (a["dy_um"][p] + b["dy_um"][q])
    mag = np.hypot(mx, my)
    qv = ax[2].quiver(a["obs_x_px"][p], a["obs_y_px"][p], mx, -my, mag,
                      cmap="viridis", angles="xy", scale_units="xy",
                      scale=1 / 0.6, width=.004)
    ax[2].set_aspect("equal"); ax[2].invert_yaxis()
    ax[2].set_xticks([]); ax[2].set_yticks([])
    ax[2].set_title("deviation averaged over both images\n"
                    "= the reproducible part", fontsize=10)
    fig.colorbar(qv, ax=ax[2], fraction=.046).set_label("|mean deviation|, um",
                                                        fontsize=9)

    fig.suptitle("MT035A Week 2 - systematic vs random deviation "
                 f"(systematic {sys_tot:.0f} um, random {noi_tot:.0f} um)",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    print("written:", out_path)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
