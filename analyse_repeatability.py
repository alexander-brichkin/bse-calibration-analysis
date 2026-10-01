#!/usr/bin/env python3
"""
MT035A - is the measured deviation systematic or random?

The brief asks whether the observed deviations appear systematic or random.
A single image cannot answer that: any residual left after the similarity fit
looks like a pattern. Several independent snapshots of the same plate can,
because a deviation caused by the machine repeats from image to image while
detection noise does not.

With N images the variance of each spot's deviation splits the usual way:

    between-image mean      ->  the reproducible (systematic) field
    scatter about that mean ->  the random part

For N = 2 this reduces to the familiar pair of identities,
var(d1 - d2) = 2*var(noise) and cov(d1, d2) = var(systematic), and the script
also prints the pairwise correlations so the two agree.

Usage
-----
    python3 analyse_repeatability.py A-deviations.csv B-deviations.csv [...] \
            [-o output.png]

Any number of deviation tables from two upwards.

Requires: numpy, matplotlib.
"""

import itertools
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    t = np.genfromtxt(path, delimiter=",", names=True)
    name = os.path.basename(path).replace("-deviations.csv", "")
    return name, {int(k): n for n, k in enumerate(t["ref_idx"])}, t


def main(paths, out_path="week2-repeatability.png"):
    sets = [load(p) for p in paths]
    names = [s[0] for s in sets]
    n_img = len(sets)
    if n_img < 2:
        sys.exit("need at least two deviation tables")

    common = sorted(set.intersection(*[set(s[1]) for s in sets]))

    # A spot recovered by the second pass is measured with a window centroid
    # instead of a blob centroid. Two estimators have two different scatters,
    # so mixing them inflates the random term and washes out the systematic
    # one. The variance split uses only spots the detector itself found in
    # every image; the recovered ones are counted but not pooled.
    n_all = len(common)
    if all("detected" in s[2].dtype.names for s in sets):
        common = [r for r in common
                  if all(t["detected"][idx[r]] > 0.5 for _, idx, t in sets)]
    n_rec = n_all - len(common)
    print(f"{n_img} images | matched in all of them: {n_all}"
          f" | measured by the detector in all of them: {len(common)}"
          f" ({n_rec} recovered by the second pass, excluded from the split)")

    # (N, M, 2) stack of deviation vectors, micrometres
    dev = np.array([[[t["dx_um"][idx[r]], t["dy_um"][idx[r]]] for r in common]
                    for _, idx, t in sets])

    mean = dev.mean(axis=0)                      # reproducible field
    resid = dev - mean                           # what does not repeat
    # Unbiased split. With d_i = s + n_i, the scatter about the per-spot mean
    # has expectation (N-1)*sigma^2, so sigma^2 is that sum divided by N-1 --
    # not multiplied by N/(N-1), which overstates the noise N-fold and drives
    # the systematic term to zero.
    noise_var = (resid ** 2).sum(axis=0).mean(axis=0) / (n_img - 1)
    total_var = (dev ** 2).mean(axis=(0, 1))
    syst_var = np.maximum(total_var - noise_var, 0.0)

    for ax, lab in enumerate("XY"):
        print(f"  d{lab}: systematic sd {np.sqrt(syst_var[ax]):5.1f} um | "
              f"random sd {np.sqrt(noise_var[ax]):5.1f} um")
    syst_tot = np.sqrt(syst_var.sum())
    noise_tot = np.sqrt(noise_var.sum())
    print(f"  => systematic {syst_tot:.1f} um vs random {noise_tot:.1f} um")

    print("\n  pairwise correlation of the same spot's deviation:")
    pairs = []
    for i, j in itertools.combinations(range(n_img), 2):
        cx = np.corrcoef(dev[i, :, 0], dev[j, :, 0])[0, 1]
        cy = np.corrcoef(dev[i, :, 1], dev[j, :, 1])[0, 1]
        pairs.append((names[i], names[j], cx, cy))
        print(f"    {names[i]} vs {names[j]}:  dX {cx:+.3f}  dY {cy:+.3f}")

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))

    i, j = 0, 1
    lim = max(np.abs(dev[[i, j], :, 0]).max(), 1) * 1.05
    ax[0].plot([-lim, lim], [-lim, lim], color="#bbbbbb", lw=1, zorder=0)
    ax[0].scatter(dev[i, :, 0], dev[j, :, 0], s=6, alpha=.35, color="#1f77b4",
                  edgecolors="none")
    ax[0].set_xlabel(f"dX in {names[i]}, um")
    ax[0].set_ylabel(f"dX in {names[j]}, um")
    ax[0].set_aspect("equal")
    ax[0].set_title(f"same spot, two images\ncorrelation {pairs[0][2]:+.2f}",
                    fontsize=10)

    first = sets[0][2]
    r_all = np.array([first["radius_mm"][sets[0][1][k]] for k in common])
    edges = np.linspace(0, r_all.max(), 9)
    mid = 0.5 * (edges[:-1] + edges[1:])
    for k in range(n_img):
        mag = np.hypot(dev[k, :, 0], dev[k, :, 1])
        rms = [np.sqrt((mag[(r_all >= lo) & (r_all < hi)] ** 2).mean())
               if ((r_all >= lo) & (r_all < hi)).sum() > 4 else np.nan
               for lo, hi in zip(edges[:-1], edges[1:])]
        ax[1].plot(mid, rms, "o-", label=names[k])
    ax[1].set_xlabel("distance from field centre, mm")
    ax[1].set_ylabel("RMS |deviation|, um")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.25)
    ax[1].set_title("the error grows towards the rim", fontsize=10)

    px = np.array([first["obs_x_px"][sets[0][1][k]] for k in common])
    py = np.array([first["obs_y_px"][sets[0][1][k]] for k in common])
    mag = np.hypot(mean[:, 0], mean[:, 1])
    q = ax[2].quiver(px, py, mean[:, 0], -mean[:, 1], mag, cmap="viridis",
                     angles="xy", scale_units="xy", scale=1 / 0.6, width=.004)
    ax[2].set_aspect("equal"); ax[2].invert_yaxis()
    ax[2].set_xticks([]); ax[2].set_yticks([])
    ax[2].set_title(f"deviation averaged over {n_img} images\n"
                    "= the reproducible part", fontsize=10)
    fig.colorbar(q, ax=ax[2], fraction=.046).set_label("|mean deviation|, um",
                                                       fontsize=9)

    fig.suptitle(f"MT035A - systematic vs random deviation over {n_img} images "
                 f"(systematic {syst_tot:.0f} um, random {noise_tot:.0f} um)",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight",
                metadata={"Software": None})
    print("\nwritten:", out_path)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "-o"]
    if "-o" in sys.argv:
        out = sys.argv[sys.argv.index("-o") + 1]
        args = [a for a in args if a != out]
        main(args, out)
    else:
        main(args)
