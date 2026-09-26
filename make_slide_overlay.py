#!/usr/bin/env python3
"""
Wide, three-panel version of the registration overlay, laid out for a 16:9
slide: the whole plate with the programmed pattern drawn on it, a zoom that
lets the audience judge the fit, and the residual field.

Usage
-----
    python3 make_slide_overlay.py SNAPSHOT.png [output.png]
"""

import os
import sys

import numpy as np
from skimage import io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bse_register as B

X0, Y0, SIZE = 330, 330, 180


def main(img_path, out_path="week2-slide-overlay.png"):
    img = io.imread(img_path)
    ref_mm, origin, _ = B.load_reference("reference_spots.csv")
    obs = B.detect(img)[:, :2]
    pitch, _ = B.lattice_pitch_and_angle(obs)
    scale, rot, trans, ii, jj, resid, ref = B.register(ref_mm, obs, pitch)
    um = 1000.0 / scale
    pred = (scale * (rot @ ref.T).T) + trans
    dev = np.hypot(*resid.T) * um

    fig, ax = plt.subplots(1, 3, figsize=(17.5, 5.6))

    ax[0].imshow(img, cmap="gray")
    ax[0].scatter(pred[:, 0], pred[:, 1], s=9, facecolors="none",
                  edgecolors="#00e5ff", linewidths=0.5)
    ax[0].add_patch(plt.Rectangle((X0, Y0), SIZE, SIZE, fill=False,
                                  edgecolor="#ff3b30", linewidth=1.8))
    ax[0].set_title(f"whole plate: {len(ii)} of {len(ref)} spots matched",
                    fontsize=13)

    ax[1].imshow(img, cmap="gray")
    ax[1].set_xlim(X0, X0 + SIZE); ax[1].set_ylim(Y0 + SIZE, Y0)
    ax[1].scatter(pred[:, 0], pred[:, 1], s=210, facecolors="none",
                  edgecolors="#00e5ff", linewidths=1.8)
    for sp in ax[1].spines.values():
        sp.set_color("#ff3b30"); sp.set_linewidth(1.8)
    ax[1].set_title("zoom: each ring sits on its own spot", fontsize=13)

    q = ax[2].quiver(obs[jj, 0], obs[jj, 1], resid[:, 0], -resid[:, 1], dev,
                     cmap="viridis", angles="xy", scale_units="xy",
                     scale=1 / 25, width=0.004)
    ax[2].set_aspect("equal"); ax[2].invert_yaxis()
    ax[2].set_title("what is left over, exaggerated x25", fontsize=13)
    cb = fig.colorbar(q, ax=ax[2], fraction=0.046, pad=0.02)
    cb.set_label("|deviation|, um", fontsize=11)
    cb.ax.tick_params(labelsize=10)

    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight", metadata={"Software": None})
    print(f"written: {out_path} | pixel {um:.3f} um | RMS "
          f"{np.sqrt((resid ** 2).sum(1).mean()) * um:.1f} um")


if __name__ == "__main__":
    main(*(sys.argv[1:3] or ["bse-snapshot-img2.png"]))
