#!/usr/bin/env python3
"""
MT035A Week 1 - figure for the programmed reference pattern.

Draws the extracted spot centres on the reference image, marks the chosen
origin spot, and shows how the image scale was established: the
nearest-neighbour spacing is the nominal 2.000 mm lattice pitch.

Usage
-----
    python3 make_reference_figure.py [pattern_image.png] [reference_spots.csv]
"""

import sys

import numpy as np
from scipy.spatial import cKDTree
from skimage import io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PITCH_MM = 2.0


def main(img_path="pattern_image.png", csv_path="reference_spots.csv",
         out_path="week1-reference-pattern.png"):
    img = io.imread(img_path)
    t = np.genfromtxt(csv_path, delimiter=",", names=True)
    origin = int(np.flatnonzero(t["is_reference"] == 1)[0])

    # The stored pixel coordinates are relative to the plot frame that was
    # cropped away in Week 1. Locate that frame again so the centroids land on
    # the spots when drawn over the original image.
    grey = img[..., :3].mean(axis=2) if img.ndim == 3 else img.astype(float)
    dark = grey < 128
    h, w = dark.shape
    rows = np.where(dark.sum(axis=1) > 0.5 * w)[0]
    cols = np.where(dark.sum(axis=0) > 0.5 * h)[0]
    ox, oy = cols.min() + 3, rows.min() + 3
    x, y = t["x_px"] + ox, t["y_px"] + oy

    pts = np.column_stack([t["x_px"], t["y_px"]])
    nn_px = cKDTree(pts).query(pts, k=2)[0][:, 1]
    # Marker centres are rounded to whole pixels, so nearest-neighbour
    # distances take two integer values. Fit the distinct row and column
    # positions instead to recover the true pitch.
    pitch_px = float(np.mean([
        np.linalg.lstsq(
            np.column_stack([np.round((np.unique(c) - np.unique(c)[0])
                                      / np.median(np.diff(np.unique(c)))),
                             np.ones(len(np.unique(c)))]),
            np.unique(c), rcond=None)[0][0]
        for c in (t["x_px"], t["y_px"])]))
    px_per_mm = pitch_px / PITCH_MM
    um_per_px = 1000.0 / px_per_mm

    fig = plt.figure(figsize=(15, 6.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1, 1], hspace=0.28,
                          wspace=0.2)

    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(img, cmap="gray")
    ax.scatter(x, y, s=4, color="#00bcd8", edgecolors="none")
    ax.scatter(x[origin], y[origin], s=170, facecolors="none",
               edgecolors="#d62728", linewidths=2.0)
    ax.annotate("origin spot, assigned (0, 0)",
                xy=(x[origin], y[origin]), xytext=(x[origin] + 130,
                                                   y[origin] - 230),
                color="#d62728", fontsize=11,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5),
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.4))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{len(x)} programmed spots, centroids extracted",
                 fontsize=13)

    ax = fig.add_subplot(gs[0, 1])
    cx, cy = x[origin], y[origin]
    half = 5.5 * pitch_px
    ax.imshow(img, cmap="gray")
    ax.set_xlim(cx - half, cx + half); ax.set_ylim(cy + half, cy - half)
    ax.scatter(x, y, s=26, color="#00bcd8", edgecolors="none")
    ax.plot([cx - 0.5 * pitch_px, cx + 0.5 * pitch_px],
            [cy + 4.5 * pitch_px] * 2, color="#d62728", lw=2.6,
            solid_capstyle="butt")
    ax.text(cx, cy + 4.15 * pitch_px,
            f"{PITCH_MM:.3f} mm = {pitch_px:.2f} px", color="#d62728",
            fontsize=10, ha="center", va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("zoom on the origin: the lattice step", fontsize=11)

    ax = fig.add_subplot(gs[1, 1])
    vals, counts = np.unique(nn_px, return_counts=True)
    ax.bar(vals, counts, width=0.22, color="#0a7e92")
    ax.axvline(pitch_px, color="#d62728", lw=1.8,
               label=f"fitted pitch {pitch_px:.3f} px")
    ax.set_xlim(vals.min() - 0.8, vals.max() + 0.8)
    ax.set_xlabel("nearest-neighbour distance, px")
    ax.set_ylabel("spots")
    ax.legend(fontsize=9)
    ax.set_title("marker centres are rounded to whole pixels\n"
                 "so the spacing takes two integer values", fontsize=11)

    ax = fig.add_subplot(gs[:, 2]); ax.axis("off")
    span_x = (x.max() - x.min()) / px_per_mm
    span_y = (y.max() - y.min()) / px_per_mm
    ax.text(0, 1,
            "REFERENCE PATTERN - WEEK 1\n\n"
            f"spots extracted     {len(x)}\n"
            f"origin spot id      {int(t['id'][origin])}\n"
            f"lattice pitch       {pitch_px:.4f} px\n"
            f"nominal pitch       {PITCH_MM:.3f} mm\n"
            f"=> image scale      {px_per_mm:.4f} px/mm\n"
            f"=> pixel size       {um_per_px:.3f} um\n\n"
            f"field width         {span_x:.2f} mm\n"
            f"field height        {span_y:.2f} mm\n"
            "build envelope      100 mm diameter\n\n"
            "Each spot carries an ID and a position\n"
            "in mm relative to the origin spot.\n\n"
            "NOTE 1  the spot diameter in this file\n"
            "is identical for every spot (zero\n"
            "spread), so it is a plotting marker,\n"
            "not a physical size. Shape descriptors\n"
            "must not be compared against it.\n\n"
            "NOTE 2  marker centres are rounded to\n"
            "whole pixels: +/-0.5 px = +/-28 um of\n"
            "quantisation on every reference\n"
            "coordinate, sd about 16 um per axis.\n"
            "That is a floor on any deviation\n"
            "measured against this pattern.",
            family="monospace", fontsize=13, va="top")

    fig.suptitle("MT035A Week 1 - programmed reference pattern and its scale",
                 fontsize=13, y=1.0)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    print("written:", out_path)
    print(f"  {len(x)} spots | pitch {pitch_px:.3f} px | "
          f"{px_per_mm:.4f} px/mm | pixel {um_per_px:.3f} um | "
          f"field {span_x:.2f} x {span_y:.2f} mm")


if __name__ == "__main__":
    main(*sys.argv[1:4])
