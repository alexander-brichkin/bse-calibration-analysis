#!/usr/bin/env python3
"""
Self-test: build a synthetic calibration pattern and a synthetic BSE image
with a KNOWN similarity transform, then check that the pipeline recovers it.

Nothing from the course data set is needed, so this is the fastest way to
confirm an installation works - and it doubles as the robustness harness the
project brief asks for, because it re-runs the same recovery under degraded
image quality.

Usage
-----
    python3 selftest.py            # run every case
    python3 selftest.py --keep     # also write selftest-preview.png

Exit code is 0 if every case passes, 1 otherwise.
"""

import sys

import numpy as np
from scipy import ndimage as ndi

import bse_register as B

# ---- the ground truth we will try to recover -------------------------------
TRUE_SCALE = 10.0          # px per mm
TRUE_ANGLE = 0.35          # degrees
TRUE_TX, TRUE_TY = 451.0, 448.0
N_SIDE = 45                # 45 x 45 lattice, like the real pattern
PLATE_R_MM = 44.0
IMG = 900

# ---- how close is close enough ---------------------------------------------
TOL_SCALE_REL = 2e-3       # 0.2 % on the scale
TOL_ANGLE_DEG = 0.05
TOL_SHIFT_PX = 0.5
MIN_MATCH_FRAC = 0.90


def synthetic_reference():
    """The programmed pattern: a square lattice clipped to a round plate."""
    k = np.arange(N_SIDE) - (N_SIDE - 1) / 2.0
    gx, gy = np.meshgrid(k * B.PITCH_MM, k * B.PITCH_MM)
    xy = np.column_stack([gx.ravel(), gy.ravel()])
    xy = xy[np.hypot(xy[:, 0], xy[:, 1]) <= PLATE_R_MM]
    origin = int(np.argmin(np.hypot(xy[:, 0], xy[:, 1])))
    return xy, origin


def synthetic_image(ref_mm, scale, angle_deg, tx, ty,
                    contrast=0.45, noise=0.02, blur=0.0, gradient=0.10,
                    seed=0):
    """Render a BSE-like image of that pattern under a known transform."""
    rng = np.random.default_rng(seed)
    th = np.radians(angle_deg)
    rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    flipped = ref_mm * np.array([1, -1])              # image convention, y down
    pts = (scale * (rot @ flipped.T).T) + np.array([tx, ty])

    img = np.full((IMG, IMG), 0.62)

    # uneven illumination, as a real BSE field has
    yy, xx = np.mgrid[0:IMG, 0:IMG]
    img += gradient * ((xx - IMG / 2) / IMG + 0.6 * (yy - IMG / 2) / IMG)

    # the melted spots: darker gaussian blobs on the lattice
    sigma = 0.09 * scale * B.PITCH_MM
    for px, py in pts:
        i0, i1 = int(py - 4 * sigma), int(py + 4 * sigma) + 1
        j0, j1 = int(px - 4 * sigma), int(px + 4 * sigma) + 1
        if i0 < 0 or j0 < 0 or i1 > IMG or j1 > IMG:
            continue
        sy, sx = np.mgrid[i0:i1, j0:j1]
        img[i0:i1, j0:j1] -= contrast * np.exp(
            -(((sx - px) ** 2 + (sy - py) ** 2) / (2 * sigma ** 2)))

    # outside the plate is dark, as on the real snapshots
    plate = np.hypot(xx - IMG / 2, yy - IMG / 2) <= 0.995 * IMG / 2
    img = np.where(plate, img, 0.02)

    if blur:
        img = ndi.gaussian_filter(img, blur)
    img += rng.normal(0, noise, img.shape)
    return np.clip(img, 0, 1) * 255.0


def recover(img, ref_mm):
    obs = B.detect(img)[:, :2]
    pitch, _ = B.lattice_pitch_and_angle(obs)
    scale, rot, trans, ii, jj, resid, _ = B.register(ref_mm, obs, pitch)
    angle = np.degrees(np.arctan2(rot[1, 0], rot[0, 0]))
    return dict(n_detected=len(obs), n_matched=len(ii), scale=scale,
                angle=angle, tx=trans[0], ty=trans[1],
                rms_px=float(np.sqrt((resid ** 2).sum(1).mean())))


def check(name, res, n_ref, strict=True):
    ds = abs(res["scale"] - TRUE_SCALE) / TRUE_SCALE
    da = abs(((res["angle"] - TRUE_ANGLE + 45) % 90) - 45)
    dt = np.hypot(res["tx"] - TRUE_TX, res["ty"] - TRUE_TY)
    frac = res["n_matched"] / n_ref
    ok = (ds < TOL_SCALE_REL and da < TOL_ANGLE_DEG and dt < TOL_SHIFT_PX
          and frac > MIN_MATCH_FRAC)
    mark = "PASS" if ok else ("FAIL" if strict else "soft")
    print(f"  [{mark}] {name:<34} matched {res['n_matched']:4d}/{n_ref} "
          f"({frac:5.1%}) | scale err {ds * 1e2:6.3f} % | angle err "
          f"{da:6.4f} deg | shift err {dt:5.2f} px | RMS {res['rms_px']:.3f} px")
    return ok or not strict


def main(keep=False):
    ref_mm, _ = synthetic_reference()
    n_ref = len(ref_mm)
    print(f"synthetic pattern: {n_ref} spots, {B.PITCH_MM} mm pitch")
    print(f"ground truth: scale {TRUE_SCALE} px/mm, angle {TRUE_ANGLE} deg, "
          f"shift ({TRUE_TX}, {TRUE_TY}) px\n")

    cases = [
        ("clean",                     dict(), True),
        ("low contrast (0.18)",       dict(contrast=0.18), True),
        ("noisy (sigma 0.06)",        dict(noise=0.06), True),
        ("blurred (sigma 1.5 px)",    dict(blur=1.5), True),
        ("strong background gradient", dict(gradient=0.30), True),
        ("low contrast + noisy",      dict(contrast=0.18, noise=0.06), False),
    ]

    passed = True
    for i, (name, kw, strict) in enumerate(cases):
        img = synthetic_image(ref_mm, TRUE_SCALE, TRUE_ANGLE, TRUE_TX, TRUE_TY,
                              seed=i, **kw)
        try:
            res = recover(img, ref_mm)
        except Exception as exc:                       # noqa: BLE001
            print(f"  [FAIL] {name:<34} raised {type(exc).__name__}: {exc}")
            passed = passed and not strict
            continue
        passed = check(name, res, n_ref, strict) and passed
        if keep and i == 0:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.imshow(img, cmap="gray")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title("selftest: synthetic BSE image")
            fig.savefig("selftest-preview.png", dpi=110, bbox_inches="tight")
            print("         wrote selftest-preview.png")

    print("\n" + ("ALL CHECKS PASSED" if passed else "SOME CHECKS FAILED"))
    print("(rows marked 'soft' are deliberately degraded beyond the point "
          "where the method is expected to stay reliable)")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main(keep="--keep" in sys.argv))
