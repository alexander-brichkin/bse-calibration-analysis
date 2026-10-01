#!/usr/bin/env python3
"""
Single-panel registration overlay, drawn in the convention the group agreed:
origin (0,0) at the centre of the BSE image, axes in pixels from that centre.

Observed spot centres are drawn as open green circles, the programmed pattern
after the fitted similarity transform as filled red dots. Where the two sit on
top of each other the registration is good; a green circle with no red dot in
it is a spot the matching did not take, which is the honest signal that part of
the field is not being measured.

Usage
-----
    python3 make_registration_overlay.py IMAGE [IMAGE ...]

Writes, next to each input:
    <image>-registration.png              the overlay
    <image>-matched_spot_residuals.csv    one row per matched spot, in the same
                                          columns and the same centre-origin
                                          frame the group is using, so tables
                                          from different pipelines line up
                                          column for column

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


def main(paths, ref_path="reference_spots.csv"):
    ref_mm, _, ref_px_per_mm = B.load_reference(ref_path)
    for path in paths:
        img = io.imread(path)
        if img.ndim == 3:
            img = img[..., :3].mean(axis=2)
        h, w = img.shape
        cx, cy = w / 2.0, h / 2.0

        obs = B.detect(img)[:, :2]
        pitch, _ = B.lattice_pitch_and_angle(obs)
        s, rot, t, ii, jj, resid, ref = B.register(ref_mm, obs, pitch)
        um = 1000.0 / s
        pred = (s * (rot @ ref.T).T) + t

        # into the agreed frame: origin at the image centre, y up
        ox, oy = obs[:, 0] - cx, -(obs[:, 1] - cy)
        px, py = pred[ii, 0] - cx, -(pred[ii, 1] - cy)
        dev = np.hypot(*resid.T) * um

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(img, cmap="gray", extent=[-cx, w - cx, -(h - cy), cy])
        ax.scatter(ox, oy, s=26, facecolors="none", edgecolors="#2ca02c",
                   linewidths=0.8, label="Observed BSE centres")
        ax.scatter(px, py, s=10, color="#e8000b",
                   label="Registered programmed reference")
        ax.axhline(0, color="yellow", lw=1)
        ax.axvline(0, color="yellow", lw=1)
        ax.plot(0, 0, "+", color="yellow", ms=14, mew=2, label="Origin (0,0)")
        ax.set_xlabel("X from BSE image centre (px)")
        ax.set_ylabel("Y from BSE image centre (px)")
        name = os.path.splitext(os.path.basename(path))[0]
        ax.set_title(f"{name}: programmed-pattern registration\n"
                     f"{len(ii)} of {len(obs)} detected spots matched | "
                     f"pixel {um:.3f} um | RMS "
                     f"{np.sqrt((resid ** 2).sum(1).mean()) * um:.1f} um",
                     fontsize=12)
        ax.legend(loc="lower left", fontsize=9, framealpha=.9)
        # per-spot table in the group's agreed frame and column order,
        # with the physical units the weekly plan asks for appended
        csv_out = f"{name}-matched_spot_residuals.csv"
        # the programmed position in the REFERENCE image's own pixels,
        # origin at the pattern centre - the same quantity, in the same units,
        # that the other pipelines put in this column
        rx = ref[ii, 0] * ref_px_per_mm
        ry = -ref[ii, 1] * ref_px_per_mm
        np.savetxt(
            csv_out,
            np.column_stack([rx, ry, ox[jj], oy[jj], px, py,
                             resid[:, 0], -resid[:, 1], np.hypot(*resid.T),
                             resid[:, 0] * um, -resid[:, 1] * um, dev]),
            delimiter=",", fmt="%.4f",
            header="reference_x_px,reference_y_px,observed_x_px,observed_y_px,"
                   "registered_reference_x_px,registered_reference_y_px,"
                   "dx_residual_px,dy_residual_px,residual_magnitude_px,"
                   "dx_residual_um,dy_residual_um,residual_magnitude_um",
            comments="")

        out = f"{name}-registration.png"
        fig.savefig(out, dpi=110, bbox_inches="tight",
                    metadata={"Software": None})
        plt.close(fig)
        print(f"{name}: matched {len(ii)}/{len(obs)} | pixel {um:.3f} um | "
              f"RMS {np.sqrt((resid ** 2).sum(1).mean()) * um:.1f} um | "
              f"median {np.median(dev):.1f} | max {dev.max():.0f}\n"
              f"        -> {out}, {csv_out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
