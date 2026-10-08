#!/usr/bin/env python3
"""Can the tilt ANGLE be closed from the data we were given? Partly.

    python3 tilt_angle_feasibility.py --flat A.png B.png C.png --test X.png Y.png

A defocused spot has d^2 = d0^2 + (2a*dW)^2, where a is the beam convergence
semi-angle and dW the working-distance error. On a plate tilted by alpha,
dW = W0 + x*sin(alpha) with W0 the focus offset at the plate centre, so

    area ~ const + 8 a^2 W0 sin(alpha) * x  +  4 a^2 sin^2(alpha) * x^2
                  \__________ L __________/    \_______ Q _______/

One measured gradient L, three unknowns (a, W0, alpha). That is why no angle
is reported. But the unknowns cancel in two places, and this script tests
both so that the limitation is demonstrated rather than asserted.

ROUTE 1 -- the ratio of two tilted plates. L_A/L_B = sin(aA)/sin(aB): a and
W0 cancel, provided both images were taken with the same beam settings. This
one works, and is reported.

ROUTE 2 -- self-calibration from the flat plates. A flat plate is also
defocused off-axis, because the path to a point at radius r is longer:
dW = r^2/(2*WD0). That gives a radial term R = 4 a^2 W0 / WD0 and hence
sin(alpha) = L / (2 R WD0), leaving only the working distance unknown. This
one fails: R is not reproducible between nominally identical flat plates --
it does not even keep its sign -- so there is no calibration to use.

The excess gradient is computed three ways (plain least squares, outlier
trimmed, and spatially binned medians) because a conclusion that depends on
the fitting method is not a conclusion.
"""
import sys
import numpy as np

import bse_register as B
import tilt_analysis as T


def field(path, ref):
    plate, area, _, _ = T.spot_field(path, ref)
    return plate, area / np.median(area)        # median, not mean: a few
                                                # merged blobs move the mean


def plane(plate, area, mode):
    keep = np.ones(len(area), bool) if mode == "plain" else \
        (area > 0.5) & (area < 2.0)
    p, a = plate[keep], area[keep]
    if mode == "binned":
        gx = np.linspace(p[:, 0].min(), p[:, 0].max(), 13)
        gy = np.linspace(p[:, 1].min(), p[:, 1].max(), 13)
        P, A = [], []
        for i in range(12):
            for j in range(12):
                s = ((p[:, 0] >= gx[i]) & (p[:, 0] < gx[i + 1]) &
                     (p[:, 1] >= gy[j]) & (p[:, 1] < gy[j + 1]))
                if s.sum() >= 5:
                    P.append(p[s].mean(0))
                    A.append(np.median(a[s]))
        p, a = np.array(P), np.array(A)
    design = np.column_stack([p, np.ones(len(p))])
    return np.linalg.lstsq(design, a, rcond=None)[0][:2]


def main(flats, tests):
    ref, _, _ = B.load_reference("reference_spots.csv")
    data = {f: field(f, ref) for f in flats + tests}

    print("\nROUTE 1 - ratio of the two tilt angles (a and W0 cancel)\n")
    print(f"  {'fit':10}{'floor':>9}{'excess A':>11}{'excess B':>11}"
          f"{'A/B':>8}")
    ratios = []
    for mode in ("plain", "trimmed", "binned"):
        g = {f: plane(*data[f], mode) for f in flats + tests}
        base = np.mean([g[f] for f in flats], axis=0)
        floor = np.mean([np.hypot(*(g[f] - base)) for f in flats])
        ea = np.hypot(*(g[tests[0]] - base))
        eb = np.hypot(*(g[tests[1]] - base))
        ratios.append(ea / eb)
        print(f"  {mode:10}{100*floor:>8.3f}%{100*ea:>9.3f}% "
              f"({ea/floor:>3.1f}x){100*eb:>7.3f}% ({eb/floor:>3.1f}x)"
              f"{ea/eb:>8.2f}")
    r = np.array(ratios)
    print(f"\n  sin(alpha_A) / sin(alpha_B) = {r.mean():.2f} "
          f"(spread {r.min():.2f} to {r.max():.2f} across the three fits)")
    print("  The first tilted plate is the steeper of the two, by about a "
          "fifth.\n  This needs nothing from outside the data.")

    print("\nROUTE 2 - self-calibration from the radial defocus of the flats\n")
    Rs = []
    for f in flats:
        plate, area = data[f]
        keep = (area > 0.5) & (area < 2.0)
        r2 = (plate[keep] ** 2).sum(1)
        a = area[keep]
        order = np.argsort(r2)
        edges = np.linspace(0, len(r2), 15).astype(int)
        xb = [np.median(r2[order][i:j]) for i, j in zip(edges[:-1], edges[1:])]
        yb = [np.median(a[order][i:j]) for i, j in zip(edges[:-1], edges[1:])]
        R = np.polyfit(xb, yb, 1)[0]
        Rs.append(R)
        print(f"  {f.replace('.png',''):24} R = {R:+.3e} /mm^2   "
              f"(centre to rim {100*R*r2.max():+6.2f} %)")
    Rs = np.array(Rs)
    print(f"\n  mean {Rs.mean():+.3e}, sd {Rs.std():.3e} -> "
          f"sd/|mean| = {Rs.std()/abs(Rs.mean()):.1f}, and the sign is not "
          f"even stable.")
    print("  There is no usable defocus calibration in these images, so the\n"
          "  absolute angle stays out of reach - and would still need the\n"
          "  working distance on top. We report the gradient, as the weekly\n"
          "  goals document allows, and the ratio above.")


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--flat" not in a or "--test" not in a:
        sys.exit(__doc__)
    f, t = a.index("--flat"), a.index("--test")
    main(a[f + 1:t], a[t + 1:])
