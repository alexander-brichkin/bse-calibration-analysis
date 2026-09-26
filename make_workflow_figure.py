#!/usr/bin/env python3
"""
MT035A Week 2 - the workflow figure the brief requires:
original image -> processed image -> identified geometry -> overlay with the
reference pattern, shown on the same region of the same snapshot.

Usage
-----
    python3 make_workflow_figure.py SNAPSHOT.png [output.png]
"""

import sys

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, io, measure, segmentation
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bse_register as B

X0, Y0, SIZE = 330, 330, 180          # the region shown in every panel


def main(img_path, out_path="week2-workflow.png"):
    img = io.imread(img_path)
    g = img.astype(float)
    if g.ndim == 3:
        g = g[..., :3].mean(axis=2)
    g /= 255.0

    # reproduce the intermediate stages of detect()
    disc = ndi.binary_fill_holes(g > 0.12)
    lab = measure.label(disc)
    if lab.max() > 1:
        counts = np.bincount(lab.ravel()); counts[0] = 0
        disc = lab == counts.argmax()
    disc_in = ndi.binary_erosion(disc, structure=B._disk(12))
    filled = np.where(disc, g, np.median(g[disc]))
    resid = ndi.median_filter(filled, size=25) - g
    resid[~disc_in] = 0
    bw = (resid > filters.threshold_otsu(resid[disc_in])) & disc_in
    lab = measure.label(bw)
    sizes = np.bincount(lab.ravel())
    big = np.flatnonzero(sizes >= B.MIN_SPOT_PX)
    bw = ndi.binary_fill_holes(np.isin(lab, big[big != 0]))

    obs = B.detect(img)[:, :2]
    ref_mm, _, _ = B.load_reference("reference_spots.csv")
    pitch, _ = B.lattice_pitch_and_angle(obs)
    scale, rot, trans, ii, jj, _, ref = B.register(ref_mm, obs, pitch)
    pred = (scale * (rot @ ref.T).T) + trans

    sl = (slice(Y0, Y0 + SIZE), slice(X0, X0 + SIZE))
    inside = ((obs[:, 0] >= X0) & (obs[:, 0] < X0 + SIZE) &
              (obs[:, 1] >= Y0) & (obs[:, 1] < Y0 + SIZE))
    inside_ref = ((pred[:, 0] >= X0) & (pred[:, 0] < X0 + SIZE) &
                  (pred[:, 1] >= Y0) & (pred[:, 1] < Y0 + SIZE))

    fig, ax = plt.subplots(1, 4, figsize=(17, 4.7))

    ax[0].imshow(g[sl], cmap="gray")
    ax[0].set_title("1. original BSE image\nbackground drifts across the field",
                    fontsize=11)

    ax[1].imshow(resid[sl], cmap="magma")
    ax[1].set_title("2. background removed\n25 px median filter subtracted",
                    fontsize=11)

    bounds = segmentation.find_boundaries(measure.label(bw[sl]), mode="outer")
    shown = np.dstack([g[sl]] * 3)
    shown[bounds] = [1.0, 0.25, 0.1]
    ax[2].imshow(shown)
    ax[2].scatter(obs[inside, 0] - X0, obs[inside, 1] - Y0, s=10,
                  color="#ffd400", marker="+", linewidths=.9)
    ax[2].set_title("3. identified geometry\nspot boundaries and centroids",
                    fontsize=11)

    ax[3].imshow(g[sl], cmap="gray")
    ax[3].scatter(pred[inside_ref, 0] - X0, pred[inside_ref, 1] - Y0, s=95,
                  facecolors="none", edgecolors="#00e5ff", linewidths=1.3)
    ax[3].scatter(obs[inside, 0] - X0, obs[inside, 1] - Y0, s=10,
                  color="#ffd400", marker="+", linewidths=.9)
    ax[3].set_title("4. registered overlay\ncyan = programmed, yellow = observed",
                    fontsize=11)

    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"Workflow on a {SIZE} x {SIZE} px region "
                 f"({SIZE / scale:.1f} x {SIZE / scale:.1f} mm) of the plate",
                 fontsize=13, y=1.06)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    print("written:", out_path)


if __name__ == "__main__":
    main(*(sys.argv[1:3] or ["bse-snapshot-img2.png"]))
