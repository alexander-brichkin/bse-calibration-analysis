#!/usr/bin/env python3
"""
Decide whether two BSE snapshots were taken in the same orientation.

The calibration pattern is a square lattice, which is 4-fold symmetric, so the
lattice-orientation estimate in bse_register.py is only defined modulo 90 deg:
89.99 deg and -0.01 deg describe the same grid. To tell whether two snapshots
really differ by a quarter turn you have to compare the images themselves, not
their lattices - the plate rim, the illumination gradient and the surface
texture are not 4-fold symmetric.

This script normalises both images to zero mean and unit variance, then
cross-correlates image A against image B rotated by 0, 90, 180 and 270 deg and
reports the peak for each. The largest peak is the relative orientation.

Usage
-----
    python3 check_orientation.py A.png B.png

Requires: numpy, scipy, scikit-image.
"""

import sys

import numpy as np
from scipy import signal
from skimage import io


def grey(path):
    img = io.imread(path).astype(float)
    if img.ndim == 3:
        img = img[..., :3].mean(axis=2)
    return img


def normalise(a):
    a = a - a.mean()
    return a / (a.std() + 1e-9)


def main(path_a, path_b):
    a, b = grey(path_a), grey(path_b)
    n = min(a.shape + b.shape)
    a, b = a[:n, :n], b[:n, :n]
    print(f"compared over {n} x {n} px, images identical: "
          f"{np.array_equal(a, b)}\n")

    peaks = {}
    for k in range(4):
        rotated = np.rot90(b, k)
        corr = signal.fftconvolve(normalise(a),
                                  normalise(rotated)[::-1, ::-1], mode="same")
        peaks[90 * k] = corr.max() / a.size
        print(f"B rotated {90 * k:3d} deg : peak normalised "
              f"correlation {peaks[90 * k]:.4f}")

    best = max(peaks, key=peaks.get)
    print(f"\nbest match at {best} deg -> "
          f"{'same orientation' if best == 0 else f'B is rotated {best} deg relative to A'}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
