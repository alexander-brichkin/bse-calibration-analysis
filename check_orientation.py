#!/usr/bin/env python3
"""
Decide whether two BSE snapshots were taken in the same orientation.

The calibration pattern is a square lattice, which is 4-fold symmetric, so the
lattice-orientation estimate in bse_register.py is only defined modulo 90 deg:
89.99 deg and -0.01 deg describe the same grid. To tell whether two snapshots
really differ by a quarter turn you have to compare the images themselves, not
their lattices - the plate rim, the illumination gradient and the surface
texture are not 4-fold symmetric.

The two images need not share a frame size or a magnification. Each is first
reduced to its own plate: the disc is segmented, cropped to its bounding box
and resampled to a common size, so only the plate contributes and the
comparison is scale- and frame-invariant. Without that step two images of
different size produce near-identical scores for all four rotations and the
test silently returns noise.

Usage
-----
    python3 check_orientation.py A.png B.png

Requires: numpy, scipy, scikit-image.
"""

import sys

import numpy as np
from scipy import ndimage as ndi
from scipy import signal
from skimage import io, measure, transform

CANON = 512          # common size both plates are resampled to


def grey(path):
    img = io.imread(path).astype(float)
    if img.ndim == 3:
        img = img[..., :3].mean(axis=2)
    return img


def plate_only(img, size=CANON):
    """Crop to the plate and resample to a fixed size.

    Returns the canonical image and the disc diameter in the original pixels,
    which is what makes two differently scaled images comparable.
    """
    g = img / max(img.max(), 1.0)
    disc = ndi.binary_fill_holes(g > 0.12)
    lab = measure.label(disc)
    if lab.max() > 1:
        counts = np.bincount(lab.ravel())
        counts[0] = 0
        disc = lab == counts.argmax()
    rows = np.flatnonzero(disc.any(axis=1))
    cols = np.flatnonzero(disc.any(axis=0))
    crop = np.where(disc, img, np.nan)[rows[0]:rows[-1] + 1,
                                       cols[0]:cols[-1] + 1]
    crop = np.where(np.isfinite(crop), crop, np.nanmedian(crop))
    diameter = 2.0 * np.sqrt(disc.sum() / np.pi)
    return transform.resize(crop, (size, size), order=1,
                            anti_aliasing=True), diameter


def normalise(a):
    a = a - a.mean()
    return a / (a.std() + 1e-9)


def main(path_a, path_b):
    raw_a, raw_b = grey(path_a), grey(path_b)
    a, dia_a = plate_only(raw_a)
    b, dia_b = plate_only(raw_b)
    print(f"A {raw_a.shape} plate {dia_a:.1f} px | "
          f"B {raw_b.shape} plate {dia_b:.1f} px")
    print(f"both resampled to {CANON} x {CANON} before comparing "
          f"(plate scale ratio {dia_b / dia_a:.4f})\n")

    peaks = {}
    for k in range(4):
        rotated = np.rot90(b, k)
        corr = signal.fftconvolve(normalise(a), normalise(rotated)[::-1, ::-1],
                                  mode="same")
        peaks[90 * k] = corr.max() / a.size
        print(f"B rotated {90 * k:3d} deg : peak normalised "
              f"correlation {peaks[90 * k]:.4f}")

    order = sorted(peaks.values(), reverse=True)
    best = max(peaks, key=peaks.get)
    margin = (order[0] - order[1]) / max(order[0], 1e-9)
    print(f"\nbest match at {best} deg, "
          f"{margin:.1%} ahead of the next candidate")
    if margin < 0.20:
        print("INCONCLUSIVE: the four rotations score within 20 % of each "
              "other, so this test cannot separate them. Use the marked "
              "reference spot to fix the orientation instead.")
    elif best == 0:
        print("-> same orientation")
    else:
        print(f"-> B is rotated {best} deg relative to A")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
