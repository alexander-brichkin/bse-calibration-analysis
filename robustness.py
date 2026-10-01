#!/usr/bin/env python3
"""
MT035A Week 4 - robustness of the method, measured on a real BSE image.

The brief asks how the measurands move when image quality or method settings
change, and where the method stops being reliable. This script takes a real
snapshot, applies controlled degradations and parameter changes one at a time,
re-runs the whole pipeline on each, and records what every measurand did.

Degradations applied to the image
    contrast     grey values compressed towards the mean
    noise        additive Gaussian noise, in grey levels
    blur         Gaussian blur, to mimic defocus
    background   an added linear intensity ramp across the field
    occlusion    a wedge of the plate blanked out
    geometry     a KNOWN rotation and shift - the fit has to recover it

Method parameters swept
    bg_size      width of the background median filter
    min_spot_px  smallest blob kept
    max_ecc      blob eccentricity limit
    gate         match radius, as a fraction of the lattice pitch

Usage
-----
    python3 robustness.py [SNAPSHOT.png] [reference_spots.csv]

Writes week4-robustness.csv and week4-robustness.png.

Requires: numpy, scipy, scikit-image, matplotlib.
"""

import csv
import sys

import numpy as np
from scipy import ndimage as ndi
from skimage import io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import bse_register as B

# A run is called unreliable when any of these is exceeded. They are our
# proposal, not values given in the data.
#
# The retained fraction is measured against the CLEAN run on the same image,
# not against the full reference pattern: part of the programmed pattern falls
# outside the imaged area, so even a perfect run on this snapshot matches only
# about 90 % of the 1597 reference points. What robustness is about is how many
# of the spots a clean run finds still survive the degradation.
MIN_RETAINED = 0.95
MAX_PIXEL_ERR_PCT = 0.10
MAX_ANGLE_ERR_DEG = 0.05


def grey(img):
    g = img.astype(float)
    return g[..., :3].mean(axis=2) if g.ndim == 3 else g


def measure(img, ref_mm, **kw):
    """Run the pipeline and return the measurands, or None if it collapsed."""
    det_kw = {k: v for k, v in kw.items() if k != "gate_fraction"}
    try:
        obs = B.detect(img, **det_kw)[:, :2]
        if len(obs) < 50:
            return None
        pitch, _ = B.lattice_pitch_and_angle(obs)
        scale, rot, trans, ii, jj, resid, ref = B.register(
            ref_mm, obs, pitch, gate_fraction=kw.get("gate_fraction",
                                                     B.GATE_FRACTION))
        if len(ii) < 50:
            return None
    except Exception:                                      # noqa: BLE001
        return None
    um = 1000.0 / scale
    return dict(n_detected=len(obs), n_matched=len(ii),
                frac=len(ii) / len(ref), pixel_um=um,
                angle=np.degrees(np.arctan2(rot[1, 0], rot[0, 0])),
                tx=trans[0], ty=trans[1],
                rms_um=float(np.sqrt((resid ** 2).sum(1).mean()) * um))


# ----------------------------------------------------------- degradations

def deg_contrast(g, c):
    return np.clip((g - g.mean()) * c + g.mean(), 0, 255)


def deg_noise(g, sigma, seed=0):
    rng = np.random.default_rng(seed)
    return np.clip(g + rng.normal(0, sigma, g.shape), 0, 255)


def deg_blur(g, sigma):
    return g if sigma == 0 else ndi.gaussian_filter(g, sigma)


def deg_background(g, amp):
    h, w = g.shape
    yy, xx = np.mgrid[0:h, 0:w]
    ramp = amp * ((xx / w - 0.5) + 0.6 * (yy / h - 0.5))
    return np.clip(g + ramp, 0, 255)


def deg_occlusion(g, frac, seed=0):
    if frac == 0:
        return g
    h, w = g.shape
    yy, xx = np.mgrid[0:h, 0:w]
    ang = np.arctan2(yy - h / 2, xx - w / 2)
    out = g.copy()
    out[(ang > -np.pi) & (ang < -np.pi + 2 * np.pi * frac)] = 0.0
    return out


def deg_geometry(g, angle_deg, shift_px):
    out = ndi.rotate(g, angle_deg, reshape=False, order=1, cval=0.0)
    return ndi.shift(out, shift_px, order=1, cval=0.0)


# ------------------------------------------------------------------- main

def plot_only(csv_path="week4-robustness.csv",
              out_path="week4-robustness.png"):
    """Redraw the figure from a finished sweep, without re-running it."""
    rows = []
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            for k in ("value", "retained", "pixel_err_pct", "angle_err_deg",
                      "rms_um"):
                r[k] = float(r[k]) if r[k] not in ("", None) else np.nan
            for k in ("complete", "accurate", "reliable"):
                r[k] = r[k] == "True"
            rows.append(r)
    draw(rows, out_path)
    print("written:", out_path)


def main(img_path="bse-snapshot-img2.png", ref_path="reference_spots.csv"):
    g0 = grey(io.imread(img_path))
    ref_mm, _, _ = B.load_reference(ref_path)

    base = measure(g0, ref_mm)
    if base is None:
        sys.exit("the baseline run failed - check the inputs")
    print(f"baseline: {base['n_matched']} matched, pixel "
          f"{base['pixel_um']:.3f} um, angle {base['angle']:+.4f} deg, "
          f"RMS {base['rms_um']:.1f} um\n")

    rows = []

    def record(family, param, value, res, angle_truth=None, note=""):
        if res is None:
            rows.append(dict(family=family, parameter=param, value=value,
                             n_detected=0, n_matched=0, frac=0.0, retained=0.0,
                             pixel_um=np.nan, pixel_err_pct=np.nan,
                             angle_err_deg=np.nan, rms_um=np.nan,
                             complete=False, accurate=False,
                             reliable=False, note="pipeline collapsed"))
            print(f"  {family:11s} {param}={value:<8} COLLAPSED")
            return
        truth = base["angle"] if angle_truth is None else angle_truth
        perr = 100.0 * (res["pixel_um"] - base["pixel_um"]) / base["pixel_um"]
        aerr = ((res["angle"] - truth + 45) % 90) - 45
        retained = res["n_matched"] / base["n_matched"]
        # Two questions, kept apart on purpose. "complete" asks whether the
        # spots are still there; "accurate" asks whether the fit built from
        # whatever survived is still right. An occlusion removes spots by
        # construction, so it fails completeness while staying accurate - that
        # is a real and useful distinction for a QA limit.
        complete = retained >= MIN_RETAINED
        accurate = (abs(perr) <= MAX_PIXEL_ERR_PCT
                    and abs(aerr) <= MAX_ANGLE_ERR_DEG)
        ok = complete and accurate
        rows.append(dict(family=family, parameter=param, value=value,
                         n_detected=res["n_detected"], n_matched=res["n_matched"],
                         frac=res["frac"], retained=retained,
                         pixel_um=res["pixel_um"],
                         pixel_err_pct=perr, angle_err_deg=aerr,
                         rms_um=res["rms_um"], complete=complete,
                         accurate=accurate, reliable=ok, note=note))
        print(f"  {family:11s} {param}={value:<8} matched {res['n_matched']:4d} "
              f"({retained:5.1%} of clean) | pixel {perr:+6.3f} % | angle {aerr:+7.4f} "
              f"deg | RMS {res['rms_um']:5.1f} um | {'ok' if ok else 'UNRELIABLE'}")

    print("image quality")
    for c in [1.0, 0.7, 0.5, 0.35, 0.25, 0.15, 0.08]:
        record("contrast", "c", c, measure(deg_contrast(g0, c), ref_mm))
    for s in [0, 5, 10, 20, 30, 45, 60]:
        record("noise", "sigma", s, measure(deg_noise(g0, s), ref_mm))
    for s in [0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0]:
        record("blur", "sigma", s, measure(deg_blur(g0, s), ref_mm))
    for a in [0, 20, 40, 60, 90, 120]:
        record("background", "amp", a, measure(deg_background(g0, a), ref_mm))
    for f in [0.0, 0.05, 0.10, 0.20, 0.30, 0.45]:
        record("occlusion", "frac", f, measure(deg_occlusion(g0, f), ref_mm))

    print("\nknown geometry - the fit must recover what we applied")
    for ang, sh in [(0.0, (0, 0)), (0.5, (0, 0)), (2.0, (0, 0)),
                    (5.0, (0, 0)), (0.0, (15, -20)), (3.0, (25, 10))]:
        res = measure(deg_geometry(g0, ang, sh), ref_mm)
        # ndi.rotate turns the image content anticlockwise in array
        # coordinates, so the fitted angle moves by -ang.
        record("geometry", f"rot{ang}_shift{sh[0]}+{sh[1]}", ang, res,
               angle_truth=base["angle"] - ang,
               note="recovery of an applied transform")

    print("\nmethod parameters")
    for v in [11, 17, 25, 35, 45, 61]:
        record("bg_size", "px", v, measure(g0, ref_mm, bg_size=v))
    for v in [2, 4, 8, 16, 32, 64]:
        record("min_spot_px", "px", v, measure(g0, ref_mm, min_spot_px=v))
    for v in [0.55, 0.70, 0.85, 0.95, 0.99]:
        record("max_ecc", "e", v, measure(g0, ref_mm, max_ecc=v))
    for v in [0.20, 0.30, 0.45, 0.60, 0.75]:
        record("gate", "frac", v, measure(g0, ref_mm, gate_fraction=v))

    # ------------------------------------------------------------- output
    cols = ["family", "parameter", "value", "n_detected", "n_matched", "frac",
            "retained", "pixel_um", "pixel_err_pct", "angle_err_deg", "rms_um",
            "complete", "accurate", "reliable", "note"]
    with open("week4-robustness.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([
                f"{r[c]:.4f}" if isinstance(r[c], float) and np.isfinite(r[c])
                else ("" if isinstance(r[c], float) else r[c])
                for c in cols])

    draw(rows)

    n_ok = sum(r["reliable"] for r in rows)
    n_acc = sum(r["accurate"] for r in rows)
    print(f"\n{n_ok} of {len(rows)} runs were both complete and accurate; "
          f"{n_acc} stayed accurate on whatever spots survived")
    print("written: week4-robustness.csv, week4-robustness.png")


def draw(rows, out_path="week4-robustness.png"):
    """Redraw the sweep.

    Four things this figure has to get right, because the first version of it
    got all four wrong:

    * The geometry family was simply missing from the panel list, and it is
      the one family that moves the measurand: a 5 deg rotation costs 0.9 %
      of pixel size, three orders of magnitude more than anything else here.
    * A run that detects nothing has no pixel size, so there is no error to
      plot. Drawing a zero there turns a total failure into a perfect score.
    * A swept parameter normally fails at one end only, so shading from the
      first failure to the edge of the axis paints most of a usable range as
      unreliable. Failing runs are marked one by one instead.
    * The retained fraction is measured against the clean run and goes ABOVE
      100 %: mild blur and reduced contrast round the blobs and let more of
      them through the shape filter. The clean run is not the best detection,
      so the axis has to show that rather than clip it away.
    """
    units = {"contrast": "c", "noise": "sigma", "blur": "sigma",
             "background": "amp", "occlusion": "frac",
             "geometry": "rotation, deg", "bg_size": "px",
             "min_spot_px": "px", "max_ecc": "e", "gate": "frac"}
    nm = lambda r: int(r["n_matched"])
    FLOOR = 1e-3

    panels = []
    for fam in units:
        sel = [r for r in rows if r["family"] == fam]
        if fam == "geometry":
            line = [r for r in sel if "shift0+0" in r["parameter"]]
            extra = [r for r in sel if "shift0+0" not in r["parameter"]]
        else:
            line, extra = sel, []
        panels.append((fam, sorted(line, key=lambda r: r["value"]), extra))

    fig, axes = plt.subplots(2, 5, figsize=(19, 7.8))
    for col, (ax, (fam, sel, extra)) in enumerate(zip(axes.ravel(), panels)):
        x = [r["value"] for r in sel]
        ax.plot(x, [100 * r["retained"] for r in sel], "o-", color="#0a7e92",
                zorder=4)
        ax.axhline(100, color="#0a7e92", lw=0.8, ls=":", alpha=0.6)
        ax.axhline(100 * MIN_RETAINED, color="#c1660d", lw=1, ls="--")
        ax.set_ylim(-7, 115)
        ax.set_xlabel(f"{fam}  ({units[fam]})", fontsize=9)
        ax.set_title(fam, fontsize=11)
        ax.grid(alpha=0.2)
        if col % 5 == 0:
            ax.set_ylabel("spots retained, % of clean run", color="#0a7e92",
                          fontsize=9)
        else:
            ax.tick_params(labelleft=False)

        ax2 = ax.twinx()
        live = [r for r in sel if nm(r) > 0 and np.isfinite(r["pixel_err_pct"])]
        ax2.plot([r["value"] for r in live],
                 [max(abs(r["pixel_err_pct"]), FLOOR) for r in live],
                 "s--", color="#b03030", ms=4, zorder=3)
        ax2.axhline(MAX_PIXEL_ERR_PCT, color="#b03030", lw=0.8, ls=":")
        ax2.set_yscale("log")
        ax2.set_ylim(FLOOR * 0.55, 300)
        if col % 5 == 4:
            ax2.set_ylabel("|pixel-size error|, %", color="#b03030",
                           fontsize=9)
        else:
            ax2.tick_params(labelright=False)

        # Runs that found nothing: no measurand exists, so nothing is drawn on
        # the red axis. The x position is marked instead.
        dead = [r for r in sel if nm(r) == 0]
        for r in dead:
            ax.axvline(r["value"], color="#777", lw=1, ls=":", zorder=1)
        if dead:
            ax.plot([r["value"] for r in dead], [0] * len(dead), "v",
                    color="#333", ms=8, zorder=6)

        # Runs that survived but fell outside the limits.
        bad = [r for r in sel if not r["reliable"] and nm(r) > 0]
        ax.plot([r["value"] for r in bad], [100 * r["retained"] for r in bad],
                "x", color="#c1660d", ms=9, mew=2, zorder=7)

        for r in extra:
            ax.plot([r["value"]], [100 * r["retained"]], "D", color="#0a7e92",
                    mfc="none", mew=1.6, ms=8, zorder=7)
        if extra:
            ax.set_title("geometry  (open: rotation + translation)",
                         fontsize=10)

    handles = [
        mlines.Line2D([], [], color="#0a7e92", marker="o",
                      label="spots retained, % of the clean run"),
        mlines.Line2D([], [], color="#b03030", marker="s", ls="--",
                      label="|pixel-size error|, % (log; floored at 1e-3)"),
        mlines.Line2D([], [], color="#c1660d", marker="x", ls="none", mew=2,
                      label=f"outside the limits (<{100*MIN_RETAINED:.0f} % "
                            f"kept or >{MAX_PIXEL_ERR_PCT} % error)"),
        mlines.Line2D([], [], color="#333", marker="v", ls="none",
                      label="detection collapsed - no measurand"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.012))
    fig.suptitle("MT035A Week 4 - robustness: where each measurand starts to "
                 "move, and where the method stops being reliable", fontsize=13)
    fig.text(0.5, 0.028, "Retention above 100 % is real: mild blur and reduced "
             "contrast round the blobs, so more of them pass the shape filter "
             "than in the clean run.", ha="center", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.savefig(out_path, dpi=120, bbox_inches="tight",
                metadata={"Software": None})
    plt.close(fig)


if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        plot_only()
    else:
        main(*[a for a in sys.argv[1:3] if not a.startswith("--")])
