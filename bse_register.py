#!/usr/bin/env python3
"""
MT035A Week 2 - BSE calibration-pattern registration.

Detects the melt spots on a BSE snapshot of the PBF-EB calibration plate,
registers them against the nominal (programmed) pattern with a similarity
transform, and reports the deviation measures required by the group project:

    dX, dY        translation
    dtheta        rotation (defined modulo 90 deg for a square lattice)
    scale error   derived pixel size vs the nominal 2.000 mm lattice pitch
    distortion    per-spot residual field left after removing the similarity

Model:  [x', y'] = s * R(theta) * [x, y] + [tx, ty]

Usage
-----
    python3 bse_register.py SNAPSHOT.png [reference_spots.csv]

If the reference file is omitted, 'reference_spots.csv' next to this script
is used.

Outputs (written next to SNAPSHOT.png)
--------------------------------------
    SNAPSHOT-overlay.png      verification figure: nominal pattern drawn on
                              top of the snapshot, zoom, residual quiver map,
                              residual histogram, numeric summary
    SNAPSHOT-deviations.csv   one row per matched spot: nominal position in mm,
                              observed position in px and in mm (own
                              calibration, origin at the BSE reference spot),
                              dX/dY/radial deviation and local deviation vector
    SNAPSHOT-summary.csv      one row with the global fit and its statistics

Requires: numpy, scipy, scikit-image, matplotlib.
"""

import os
import sys

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import filters, io, measure, morphology, segmentation
from skimage.feature import peak_local_max
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PITCH_MM = 2.0          # nominal lattice pitch of the reference pattern, mm
MIN_SPOT_PX = 8         # smallest blob kept before shape filtering
BG_FILTER_PX = 25       # median-filter width used to estimate the background
RIM_PX = 12             # how far the working area is pulled in from the rim
AREA_LO, AREA_HI = 0.3, 4.0      # kept blob area, as a fraction of the median
MAX_ECC = 0.85          # kept blob eccentricity
ICP_ITERS = 30
GATE_FRACTION = 0.45    # match radius as a fraction of the lattice pitch
BG_MODE = "median"      # background removal: "median" or "tophat".
#   The morphological black top-hat is the textbook choice and finds ~13 more
#   spots on a clean snapshot, but it is a local min-max operation and so is
#   noise-sensitive: at sigma=30 grey levels it returns a confident 11.5 %
#   pixel-size error where the median filter is still right. Robustness beats
#   13 spots, so the median stays the default and the top-hat is an option the
#   sweep measures.
SPLIT_MERGED = True     # separate touching spots instead of discarding them
REFINE = True           # second pass: measure where a spot must be
REFINE_WIN = 0.30       # search half-window, as a fraction of the pitch
REFINE_MIN_Q = 1.5      # peak must reach this multiple of the Otsu level
REFINE_MAX_FRAC = 0.20  # refuse the second pass if it recovers more than this
#   fraction of the first pass's matches -- see the guard in `solve`.


# ---------------------------------------------------------------- detection

def _disk(radius):
    """Boolean disk structuring element (no scikit-image version coupling)."""
    r = int(radius)
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def detect(img, bg_size=BG_FILTER_PX, min_spot_px=MIN_SPOT_PX,
           area_lo=AREA_LO, area_hi=AREA_HI, max_ecc=MAX_ECC, rim_px=RIM_PX,
           bg_mode=BG_MODE, split=SPLIT_MERGED, return_maps=False):
    """Find the spots on a BSE snapshot.

    Returns an (N, 6) array: centre_x_px, centre_y_px, equivalent_diameter_px,
    area_px, eccentricity, aspect_ratio. Coordinates are in image convention
    (x right, y down).

    Eccentricity and aspect ratio measure the same elongation on two scales. A
    circle is eccentricity 0 and aspect ratio 1; the aspect ratio is the more
    readable of the two for a tilted plate, where the question is how oval a
    spot has become, so both are reported.

    The keyword arguments are the method's tunable parameters, exposed so that
    `robustness.py` can sweep them and report how far each measurand moves.
    Their defaults are the module constants and are what every other script
    uses.
    """
    g = img.astype(float)
    if g.ndim == 3:
        g = g[..., :3].mean(axis=2)
    g /= 255.0

    # 1. Isolate the circular plate. Everything outside it is background.
    disc = ndi.binary_fill_holes(g > 0.12)
    lab = measure.label(disc)
    if lab.max() > 1:
        counts = np.bincount(lab.ravel())
        counts[0] = 0
        disc = lab == counts.argmax()
    # Pull the working area in from the rim: the edge is bright and would
    # otherwise be picked up as signal.
    disc_in = ndi.binary_erosion(disc, structure=_disk(rim_px))

    # 2. Flatten the background. BSE brightness depends on atomic number,
    #    local tilt and topography, so the background level drifts across the
    #    field and a single global threshold does not work.
    #
    #    Two ways to remove it. A median filter wider than a spot estimates
    #    the background and is subtracted. A morphological black top-hat
    #    instead keeps only what is darker than its surround AND smaller than
    #    the structuring element, so the element is matched to the spot rather
    #    than to an assumption that the background varies slowly. The top-hat
    #    is the default; the median is kept because the robustness sweep
    #    measures both.
    if bg_mode == "tophat":
        resid = morphology.black_tophat(g, _disk(max(3, int(bg_size) // 3)))
    else:
        filled = np.where(disc, g, np.median(g[disc]))
        resid = ndi.median_filter(filled, size=bg_size) - g
    resid = np.asarray(resid, float).copy()
    resid[~disc_in] = 0

    # 3. Threshold the flattened map and clean it up.
    bw = (resid > filters.threshold_otsu(resid[disc_in])) & disc_in
    lab = measure.label(bw)
    sizes = np.bincount(lab.ravel())
    big = np.flatnonzero(sizes >= min_spot_px)
    bw = np.isin(lab, big[big != 0])
    bw = ndi.binary_fill_holes(bw)

    # 4. Two spots that touch form one elongated blob. Discarding it loses
    #    both; splitting it along the intensity valley between them keeps
    #    both. The markers are the local maxima of the flattened map, so the
    #    cut follows the image rather than the blob's shape.
    lab = measure.label(bw)
    if split:
        props = measure.regionprops(lab)
        med_area = np.median([p.area for p in props]) if props else 0.0
        suspect = {p.label for p in props
                   if p.eccentricity >= max_ecc or p.area > area_hi * med_area}
        if suspect:
            lab = _split_merged(lab, resid, suspect,
                                min_dist=max(3, int(round(
                                    0.5 * np.sqrt(4 * med_area / np.pi)))))

    # 5. Keep blobs that look like spots: area near the median, not elongated.
    #    This drops scratches and rim fragments.
    props = measure.regionprops(lab)
    med_area = np.median([p.area for p in props])
    keep = [p for p in props
            if area_lo * med_area <= p.area <= area_hi * med_area
            and p.eccentricity < max_ecc]

    spots = np.array([[p.centroid[1], p.centroid[0],
                       2.0 * np.sqrt(p.area / np.pi), p.area, p.eccentricity,
                       _aspect(p)]
                      for p in keep])
    if return_maps:
        return spots, resid, disc_in
    return spots


def _aspect(prop):
    """Major axis over minor axis. 1.0 is a circle; a tilt makes it grow."""
    major = getattr(prop, "axis_major_length", None)
    minor = getattr(prop, "axis_minor_length", None)
    if major is None:                       # older scikit-image
        major, minor = prop.major_axis_length, prop.minor_axis_length
    return float(major / minor) if minor > 1e-9 else float("nan")


def _split_merged(lab, resid, suspect, min_dist=7):
    """Marker-based watershed on the blobs that look like merged pairs."""
    out = lab.copy()
    nxt = int(lab.max()) + 1
    for p in measure.regionprops(lab):
        if p.label not in suspect:
            continue
        r0, c0, r1, c1 = p.bbox
        blob = lab[r0:r1, c0:c1] == p.label
        sub = resid[r0:r1, c0:c1]
        pk = peak_local_max(np.where(blob, sub, 0), min_distance=min_dist,
                            labels=blob.astype(int), exclude_border=False)
        if len(pk) < 2:
            continue
        markers = np.zeros(blob.shape, int)
        for n, (y, x) in enumerate(pk, 1):
            markers[y, x] = n
        ws = segmentation.watershed(-sub, markers, mask=blob)
        region = out[r0:r1, c0:c1]
        region[blob] = 0
        for n in range(1, len(pk) + 1):
            region[ws == n] = nxt
            nxt += 1
    return measure.label(out > 0)


def lattice_pitch_and_angle(pts):
    """Estimate lattice pitch (px) and orientation (deg) from nearest neighbours.

    The orientation of a square lattice is only defined modulo 90 degrees,
    so the returned angle is folded into [0, 90).
    """
    _, idx = cKDTree(pts).query(pts, k=5)
    vec = pts[idx[:, 1:]] - pts[:, None, :]
    dist = np.hypot(vec[..., 0], vec[..., 1]).ravel()
    pitch = np.median(dist[dist < 1.6 * np.median(dist)])

    near = np.abs(dist - pitch) < 0.25 * pitch
    ang = np.arctan2(vec[..., 1], vec[..., 0]).ravel()[near] % (np.pi / 2)
    # Circular mean on the 4-fold-wrapped angle.
    angle = np.degrees(np.angle(np.exp(4j * ang).mean()) / 4) % 90
    return pitch, angle


# ------------------------------------------------------------- registration

def affine_fit(src, dst):
    """Full six-parameter affine fit, used only as a diagnostic.

    The similarity model has four parameters and cannot represent an
    anisotropic scale or a shear. Fitting the richer model and reporting how
    far it departs from a similarity answers, from the data rather than by
    assumption, whether the extra freedom is needed. A tilted plate seen from
    above is foreshortened along the tilt axis, so a real tilt shows up here as
    anisotropy: 10 deg gives 1.5 %, 15 deg gives 3.4 %.

    Returns (scale_major, scale_minor, shear_deg, rms_px).
    """
    design = np.column_stack([src, np.ones(len(src))])
    coef = np.linalg.lstsq(design, dst, rcond=None)[0]
    lin = coef[:2].T
    sv = np.linalg.svd(lin, compute_uv=False)
    _, upper = np.linalg.qr(lin)
    shear = np.degrees(np.arctan2(upper[0, 1], upper[0, 0]))
    shear = (shear + 90) % 180 - 90          # fold onto (-90, 90]
    rms = float(np.sqrt(((dst - design @ coef) ** 2).sum(1).mean()))
    return float(sv[0]), float(sv[1]), shear, rms


def umeyama(src, dst):
    """Least-squares similarity transform src -> dst (Umeyama 1991).

    Returns (scale, rotation_matrix, translation).
    """
    mu_s, mu_d = src.mean(0), dst.mean(0)
    s_c, d_c = src - mu_s, dst - mu_d
    cov = d_c.T @ s_c / len(src)
    u, d, vt = np.linalg.svd(cov)
    corr = np.eye(2)
    if np.linalg.det(u) * np.linalg.det(vt) < 0:      # forbid a reflection
        corr[1, 1] = -1
    rot = u @ corr @ vt
    var = (s_c ** 2).sum() / len(src)
    scale = np.trace(np.diag(d) @ corr) / var
    return scale, rot, mu_d - scale * rot @ mu_s


def mutual_nearest(pred, obs, gate):
    """Mutually nearest pairs within `gate` pixels.

    Requiring the match to be nearest in both directions stops one observed
    spot from claiming several reference points, and the gate (below half a
    lattice pitch) stops a point from locking onto the neighbouring node.
    """
    d_ref, i_ref = cKDTree(obs).query(pred)
    _, i_obs = cKDTree(pred).query(obs)
    ii = [k for k in range(len(pred))
          if d_ref[k] < gate and i_obs[i_ref[k]] == k]
    ii = np.asarray(ii, dtype=int)
    return ii, i_ref[ii]


def _icp(ref, obs_px, scale, rot, gate_fraction, iters):
    """Run the matching loop from one starting pose."""
    trans = obs_px.mean(0) - scale * (rot @ ref.mean(0))
    gate = gate_fraction * scale * PITCH_MM
    for _ in range(iters):
        pred = (scale * (rot @ ref.T).T) + trans
        ii, jj = mutual_nearest(pred, obs_px, gate)
        if len(ii) < 20:
            break
        scale, rot, trans = umeyama(ref[ii], obs_px[jj])
        gate = gate_fraction * scale * PITCH_MM
    pred = (scale * (rot @ ref.T).T) + trans
    ii, jj = mutual_nearest(pred, obs_px, gate)
    return scale, rot, trans, ii, jj


def register(ref_mm, obs_px, pitch_px, iters=ICP_ITERS,
             gate_fraction=GATE_FRACTION, angle_deg=None):
    """Fit the similarity transform that maps the nominal pattern onto the spots.

    Returns (scale, rotation, translation, ref_idx, obs_idx, residual, ref_xy)
    where `residual` is obs - predicted in pixels for each matched pair and
    `ref_xy` is the reference pattern in image convention (y down).

    `angle_deg` seeds the rotation from the measured lattice orientation.
    Without it the fit starts from zero rotation and converges only while the
    pattern sits within about 2 deg of the reference axes; beyond that the
    nearest-neighbour matching locks onto the wrong lattice site and the fit
    fails while still reporting a plausible-looking angle. The lattice angle
    is only known modulo 90 deg because the lattice is square, so all four
    candidates are tried and the one that matches the most spots wins.
    """
    ref = ref_mm.copy()
    ref[:, 1] *= -1                       # CSV has y up; images have y down
    scale0 = pitch_px / PITCH_MM          # seed from the measured pitch

    angles = [0.0] if angle_deg is None else [angle_deg + 90.0 * k
                                              for k in range(4)]
    best = None
    for a in angles:
        th = np.radians(a)
        rot0 = np.array([[np.cos(th), -np.sin(th)],
                         [np.sin(th), np.cos(th)]])
        out = _icp(ref, obs_px, scale0, rot0, gate_fraction, iters)
        if best is None or len(out[3]) > len(best[3]):
            best = out

    scale, rot, trans, ii, jj = best
    pred = (scale * (rot @ ref.T).T) + trans
    return scale, rot, trans, ii, jj, obs_px[jj] - pred[ii], ref


def refine_by_prediction(resid_map, ref, scale, rot, trans, matched_ref,
                         win=REFINE_WIN, min_q=REFINE_MIN_Q, disc_in=None):
    """Second pass: measure where a spot must be, not where one was found.

    The first pass has to decide from the image alone that something is a
    spot, so it loses the faint ones, the merged ones and the ones near the
    rim. After the first fit the lattice is known, so every reference point
    that went unmatched has a predicted position and only a window of about a
    third of a pitch to search. A centroid in that window recovers spots the
    detector could not commit to -- which is exactly where the detector fails.

    Returns (extra_ref_idx, extra_xy): reference indices recovered this way
    and their measured centres in image pixels.
    """
    half = max(3, int(round(win * scale * PITCH_MM)))
    level = filters.threshold_otsu(resid_map[disc_in]) if disc_in is not None \
        else filters.threshold_otsu(resid_map)
    pred = (scale * (rot @ ref.T).T) + trans
    want = np.setdiff1d(np.arange(len(ref)), matched_ref)
    h, w = resid_map.shape

    idx, xy = [], []
    for k in want:
        cx, cy = pred[k]
        c0, c1 = int(round(cx)) - half, int(round(cx)) + half + 1
        r0, r1 = int(round(cy)) - half, int(round(cy)) + half + 1
        if c0 < 0 or r0 < 0 or c1 > w or r1 > h:
            continue
        win_map = resid_map[r0:r1, c0:c1]
        peak = win_map.max()
        if peak < min_q * level:
            continue
        mask = win_map > 0.5 * peak
        if mask.sum() < 3:
            continue
        wts = np.where(mask, win_map, 0.0)
        total = wts.sum()
        yy, xx = np.mgrid[r0:r1, c0:c1]
        ox, oy = (wts * xx).sum() / total, (wts * yy).sum() / total
        # A centroid that ran to the edge of the window is a neighbouring
        # spot leaking in, not this one.
        if np.hypot(ox - cx, oy - cy) > 0.5 * half:
            continue
        idx.append(k)
        xy.append([ox, oy])

    return np.array(idx, int), (np.array(xy, float).reshape(-1, 2))


def solve(ref_mm, obs_px, pitch_px, angle_deg=None, resid_map=None,
          disc_in=None, refine=REFINE, gate_fraction=GATE_FRACTION):
    """Register, then recover the spots the detector missed, then re-fit.

    Returns the same tuple as `register`, plus the observation array actually
    used (detected spots followed by recovered ones) and how many were
    recovered.
    """
    scale, rot, trans, ii, jj, res, ref = register(
        ref_mm, obs_px, pitch_px, gate_fraction=gate_fraction,
        angle_deg=angle_deg)
    if not refine or resid_map is None:
        return scale, rot, trans, ii, jj, res, ref, obs_px, 0

    extra_i, extra_xy = refine_by_prediction(
        resid_map, ref, scale, rot, trans, ii, disc_in=disc_in)
    if len(extra_i) == 0:
        return scale, rot, trans, ii, jj, res, ref, obs_px, 0

    obs_aug = np.vstack([obs_px, extra_xy])
    ii2 = np.concatenate([ii, extra_i])
    jj2 = np.concatenate([jj, len(obs_px) + np.arange(len(extra_i))])
    order = np.argsort(ii2)
    ii2, jj2 = ii2[order], jj2[order]
    scale2, rot2, trans2 = umeyama(ref[ii2], obs_aug[jj2])
    pred2 = (scale2 * (rot2 @ ref.T).T) + trans2
    res2 = obs_aug[jj2] - pred2[ii2]

    # The second pass searches where the first fit says a spot must be, so it
    # confirms whatever lattice the first fit chose -- including a wrong one.
    # A broken fit would come back with MORE matches than a good one, which
    # would destroy the matched fraction as a warning sign. So the refined
    # solution is only accepted if it agrees with the one it started from: the
    # scale must not move by more than a tenth of a percent, and the residual
    # must not blow up. Otherwise the first-pass answer stands, with its low
    # match count visible.
    moved = abs(scale2 / scale - 1.0)
    rms1 = float(np.sqrt((res ** 2).sum(1).mean()))
    # Compare like with like: the refined RMS is recomputed over the spots the
    # detector found, so a legitimately harder set of recovered points cannot
    # by itself trip the guard.
    was_detected = jj2 < len(obs_px)
    rms2 = float(np.sqrt((res2[was_detected] ** 2).sum(1).mean()))

    # The decisive guard. When the first fit has locked onto the wrong lattice,
    # the second pass searches wrong positions -- and under noise it finds
    # something at nearly all of them, so a broken fit comes back with MORE
    # matches than a good one and the match count stops being a warning. At
    # sigma=30 the unguarded version recovered 401 points and reported a
    # confident 11.5 % pixel-size error. A pass that has to invent a fifth of
    # the pattern is not refining a fit, it is rescuing a bad one.
    too_many = len(extra_i) > REFINE_MAX_FRAC * len(ii)
    if too_many or moved > 1e-3 or rms2 > 1.5 * rms1:
        return scale, rot, trans, ii, jj, res, ref, obs_px, 0

    return (scale2, rot2, trans2, ii2, jj2, res2, ref, obs_aug, len(extra_i))


def fit(img, ref_mm, detect_kw=None, refine=REFINE,
        gate_fraction=GATE_FRACTION):
    """Detect, then register with the measured lattice angle as the seed.

    Every script goes through here, so no caller can silently skip the
    rotation seed or the second pass. Returns the tuple from `solve` followed
    by (pitch_px, angle_deg, spots), where `spots` is the full detection
    table and `obs` inside the tuple is what the fit actually used.
    """
    spots, resid_map, disc_in = detect(img, return_maps=True,
                                       **(detect_kw or {}))
    obs = spots[:, :2]
    pitch, angle = lattice_pitch_and_angle(obs)
    out = solve(ref_mm, obs, pitch, angle_deg=angle, resid_map=resid_map,
                disc_in=disc_in, refine=refine, gate_fraction=gate_fraction)
    return out + (pitch, angle, spots)


# ------------------------------------------------------------------ reading

def _reference_pitch(x_px, y_px):
    """Lattice pitch of the reference pattern, in its own pixels.

    The reference image is a rendered plot, so every marker centre is rounded
    to a whole pixel: nearest-neighbour distances come out as 35 or 36 px and
    the median alone is not the pitch. Indexing the distinct row and column
    positions and fitting a line through them recovers the true spacing.
    """
    pitches = []
    for coord in (x_px, y_px):
        u = np.unique(coord)
        step = np.median(np.diff(u))
        k = np.round((u - u[0]) / step)
        fit = np.linalg.lstsq(np.column_stack([k, np.ones_like(k)]), u,
                              rcond=None)[0]
        pitches.append(fit[0])
    return float(np.mean(pitches))


def load_reference(path):
    """Read the nominal pattern from the CSV.

    Returns (xy_mm, origin_index, px_per_mm) with xy_mm in millimetres,
    relative to the unambiguous reference spot chosen in Week 1.

    The millimetre columns stored in the CSV are NOT used. They were written
    with a rounded image scale of 18.000 px/mm, while the pattern's actual
    lattice pitch is 35.911 px; taking them at face value stretches the whole
    reference by 0.25 % and pushes the same error straight into the derived
    BSE pixel size. Instead the pattern is re-scaled here from its own lattice
    against the nominal PITCH_MM, which is the quantity the machine was
    programmed with. If the CSV has no pixel columns, the millimetre columns
    are used as a fallback and a warning is printed.
    """
    table = np.genfromtxt(path, delimiter=",", names=True)
    cols = table.dtype.names

    if "is_reference" in cols and (table["is_reference"] == 1).any():
        origin = int(np.flatnonzero(table["is_reference"] == 1)[0])
    else:
        origin = None

    if "x_px" in cols and "y_px" in cols:
        pitch_px = _reference_pitch(table["x_px"], table["y_px"])
        px_per_mm = pitch_px / PITCH_MM
        if origin is None:
            origin = int(np.argmin(np.hypot(table["x_px"] - table["x_px"].mean(),
                                            table["y_px"] - table["y_px"].mean())))
        x0, y0 = table["x_px"][origin], table["y_px"][origin]
        xy = np.column_stack([(table["x_px"] - x0) / px_per_mm,
                              (y0 - table["y_px"]) / px_per_mm])
        print(f"reference pattern   pitch {pitch_px:.4f} px -> "
              f"{px_per_mm:.4f} px/mm ({1000 / px_per_mm:.3f} um/px), "
              f"{len(xy)} spots")
        return xy, origin, px_per_mm

    print("WARNING: no pixel columns in the reference file; falling back to "
          "its stored millimetre columns and their assumed scale")
    x_key = "x_mm" if "x_mm" in cols else "x_mm_rel"
    y_key = "y_mm" if "y_mm" in cols else "y_mm_rel"
    xy = np.column_stack([table[x_key], table[y_key]])
    if origin is None:
        origin = int(np.argmin(np.hypot(xy[:, 0], xy[:, 1])))
    return xy, origin, float("nan")


def lattice_ids(xy_mm):
    """Integer (i, j) lattice coordinates for every reference spot.

    The brief asks for spots to carry an identity rather than a row number:
    the reference spot is (0, 0), its neighbour one pitch along +x is (1, 0),
    one pitch along +y is (0, 1). `load_reference` already returns millimetres
    measured from the reference spot on a 2.000 mm lattice, so the index is
    the position divided by the pitch.

    Returns (ij, worst_slip). `worst_slip` is the largest distance, in units
    of a pitch, between a spot and the lattice node it was assigned to. It is
    reported rather than ignored: anything approaching 0.5 would mean a spot
    had been given the wrong identity, and that would quietly corrupt every
    comparison built on these IDs.
    """
    ij = np.round(xy_mm / PITCH_MM).astype(int)
    slip = float(np.abs(xy_mm / PITCH_MM - ij).max())
    return ij, slip


# ------------------------------------------------------------------ figures

def overlay_figure(img, pred_all, obs, ii, jj, resid, scale, rot, trans,
                   um_per_px, title, out_path):
    dev = np.hypot(*resid.T) * um_per_px
    theta = np.degrees(np.arctan2(rot[1, 0], rot[0, 0]))
    rms = np.sqrt((resid ** 2).sum(1).mean()) * um_per_px

    fig = plt.figure(figsize=(16.5, 8.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.18, 1, 1],
                          hspace=0.22, wspace=0.16)

    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(img, cmap="gray")
    ax.scatter(pred_all[:, 0], pred_all[:, 1], s=7, facecolors="none",
               edgecolors="#00e5ff", linewidths=0.45)
    ax.add_patch(plt.Rectangle((330, 330), 180, 180, fill=False,
                               edgecolor="#ff3b30", linewidth=1.6))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{title}: {len(ii)} matched spots\n"
                 "cyan rings = nominal pattern after the similarity fit",
                 fontsize=11)

    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(img, cmap="gray")
    ax.set_xlim(330, 510); ax.set_ylim(510, 330)
    ax.scatter(pred_all[:, 0], pred_all[:, 1], s=90, facecolors="none",
               edgecolors="#00e5ff", linewidths=1.3)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("zoom, 180 x 180 px", fontsize=10)
    for spine in ax.spines.values():
        spine.set_color("#ff3b30"); spine.set_linewidth(1.6)

    ax = fig.add_subplot(gs[1, 1])
    q = ax.quiver(obs[jj, 0], obs[jj, 1], resid[:, 0], -resid[:, 1], dev,
                  cmap="viridis", angles="xy", scale_units="xy",
                  scale=1 / 25, width=0.004)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("deviation vectors, exaggerated x25", fontsize=10)
    fig.colorbar(q, ax=ax, fraction=0.046).set_label("|deviation|, um",
                                                     fontsize=9)

    ax = fig.add_subplot(gs[0, 2])
    ax.hist(dev, bins=45, color="#1f77b4", alpha=0.85)
    ax.axvline(np.median(dev), color="#d62728", lw=1.6,
               label="median %.0f um" % np.median(dev))
    ax.axvline(np.percentile(dev, 95), color="#ff7f0e", lw=1.4, ls="--",
               label="95th pct %.0f um" % np.percentile(dev, 95))
    ax.set_xlabel("|deviation|, um"); ax.set_ylabel("spots")
    ax.legend(fontsize=8)
    ax.set_title("residual distribution", fontsize=10)

    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    ax.text(0, 1,
            f"scale          {scale:.4f} px/mm\n"
            f"pixel size     {um_per_px:.2f} um\n"
            f"rotation       {theta:+.4f} deg (mod 90)\n"
            f"dX, dY         {trans[0]:.1f}, {trans[1]:.1f} px\n\n"
            f"matched pairs  {len(ii)}\n"
            f"RMS            {rms:.1f} um\n"
            f"median         {np.median(dev):.1f} um\n"
            f"95th pct       {np.percentile(dev, 95):.1f} um\n"
            f"max            {dev.max():.0f} um\n\n"
            f"nominal pitch  {PITCH_MM:.3f} mm\n"
            f"=> BSE pixel ~ 100 um (10 px/mm)",
            family="monospace", fontsize=13.5, va="top")

    fig.suptitle("MT035A Week 2 - BSE snapshot registered to the "
                 "programmed pattern", fontsize=13)
    fig.savefig(out_path, dpi=120, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


# --------------------------------------------------------------------- main

def main(img_path, ref_path=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ref_path = ref_path or os.path.join(here, "reference_spots.csv")

    img = io.imread(img_path)
    ref_mm, ref_origin, ref_px_per_mm = load_reference(ref_path)

    (scale, rot, trans, ii, jj, resid, ref, obs, n_recovered,
     pitch_px, angle, _spots) = fit(img, ref_mm)

    # Every spot gets an identity on the lattice, not a row number: the
    # reference spot is (0, 0) and the neighbours are (1, 0), (0, 1) and so
    # on. The deviation tables are keyed on this, so a spot can be followed
    # between images and between the three implementations in the group.
    ref_ij, ij_slip = lattice_ids(ref_mm)
    print(f"lattice IDs        (i, j) from ({ref_ij[:, 0].min()}, "
          f"{ref_ij[:, 1].min()}) to ({ref_ij[:, 0].max()}, "
          f"{ref_ij[:, 1].max()}), "
          f"{len(set(map(tuple, ref_ij)))} unique of {len(ref_ij)}; "
          f"worst node slip {ij_slip:.4f} pitch")

    # Shape per matched spot. A recovered spot has a position but no blob, so
    # its shape is undefined and is written as a blank rather than a zero.
    shape = np.full((len(jj), 2), np.nan)
    det = jj < len(_spots)
    shape[det] = _spots[jj[det]][:, [5, 4]]       # aspect ratio, eccentricity

    um_per_px = 1000.0 / scale
    theta = np.degrees(np.arctan2(rot[1, 0], rot[0, 0]))
    dev = np.hypot(*resid.T) * um_per_px
    base = os.path.splitext(img_path)[0]

    print(f"spots detected     {len(_spots)} by the detector, "
          f"{len(obs)} after the second pass")
    print(f"lattice pitch      {pitch_px:.3f} px   "
          f"orientation {angle:.3f} deg (mod 90)")
    print(f"matched pairs      {len(ii)} of {len(ref)}"
          f"   ({n_recovered} recovered by the second pass)")
    print(f"scale s            {scale:.5f} px/mm  ->  "
          f"pixel size {um_per_px:.3f} um")
    print(f"rotation dtheta    {theta:+.5f} deg (mod 90)")
    print(f"translation dX,dY  {trans[0]:.2f}, {trans[1]:.2f} px  "
          f"({trans[0] * um_per_px / 1000:.3f}, "
          f"{trans[1] * um_per_px / 1000:.3f} mm)")
    print(f"residual           RMS {np.sqrt((resid ** 2).sum(1).mean()) * um_per_px:.1f} um"
          f" | median {np.median(dev):.1f}"
          f" | 95th {np.percentile(dev, 95):.1f}"
          f" | max {dev.max():.0f}")

    pred_all = (scale * (rot @ ref.T).T) + trans

    # --- express the observed spots in their own physical units -------------
    # The brief requires the two datasets to be compared as physical
    # coordinates, each converted with its own calibration, not as raw pixels.
    # The fitted scale IS the BSE image's pixel-to-mm calibration: it is
    # derived from the snapshot, not borrowed from the reference image.
    # Origin: the observed spot matched to the Week-1 reference spot; if that
    # spot was not matched, the observed centroid closest to the predicted
    # position of the reference spot is used instead.
    where = np.flatnonzero(ii == ref_origin)
    if len(where):
        origin_px = obs[jj[where[0]]]
        origin_note = "matched BSE spot"
    else:
        origin_px = obs[np.argmin(np.hypot(*(obs - pred_all[ref_origin]).T))]
        origin_note = "nearest BSE spot to the predicted reference position"
    print(f"BSE reference spot  ({origin_px[0]:.2f}, {origin_px[1]:.2f}) px "
          f"-> (0, 0) mm  [{origin_note}]")

    # Rotate the pixel offsets back onto the reference axes and divide by the
    # fitted scale: observed spot centres in mm, y up, origin at the BSE
    # reference spot.
    obs_mm = ((rot.T @ (obs[jj] - origin_px).T).T / scale) * np.array([1, -1])

    # --- radial and tangential components about the plate centre -----------
    centre = obs[jj].mean(0)
    rad_px = obs[jj] - centre
    radius_mm = np.hypot(*rad_px.T) / scale
    unit = rad_px / np.maximum(np.hypot(*rad_px.T), 1e-9)[:, None]
    d_radial = (resid * unit).sum(1) * um_per_px                # outward +
    d_tangential = (resid[:, 0] * -unit[:, 1]
                    + resid[:, 1] * unit[:, 0]) * um_per_px     # CCW +

    np.savetxt(
        base + "-deviations.csv",
        np.column_stack([ii, jj, ref[ii, 0], -ref[ii, 1],
                         pred_all[ii], obs[jj], obs_mm,
                         resid * um_per_px, dev,
                         radius_mm, d_radial, d_tangential,
                         ref_ij[ii], shape,
                         (jj < len(_spots)).astype(float)]),
        delimiter=",", fmt="%.4f",
        header="ref_idx,obs_idx,ref_x_mm,ref_y_mm,pred_x_px,pred_y_px,"
               "obs_x_px,obs_y_px,obs_x_mm,obs_y_mm,dx_um,dy_um,dist_um,"
               "radius_mm,dev_radial_um,dev_tangential_um,"
               "lattice_i,lattice_j,aspect_ratio,eccentricity,detected",
        comments="")

    sx, sy, shear, aff_rms = affine_fit(ref[ii], obs[jj])
    aniso = 100.0 * (sx / sy - 1.0)
    print(f"affine diagnostic   axis scales {sx:.5f} / {sy:.5f} "
          f"(anisotropy {aniso:+.3f} %), shear {shear:+.4f} deg, "
          f"RMS {aff_rms * um_per_px:.1f} um")

    with open(base + "-summary.csv", "w") as fh:
        fh.write("image,spots_detected,matched_pairs,reference_points,"
                 "lattice_pitch_px,scale_px_per_mm,pixel_size_um,"
                 "rotation_deg_mod90,dx_px,dy_px,rms_um,median_um,p95_um,"
                 "max_um,affine_sx,affine_sy,affine_anisotropy_pct,"
                 "affine_shear_deg,affine_rms_um,"
                 "aspect_ratio_median,aspect_ratio_p95,eccentricity_median\n")
        fh.write(f"{os.path.basename(base)},{len(obs)},{len(ii)},{len(ref)},"
                 f"{pitch_px:.4f},{scale:.5f},{um_per_px:.3f},{theta:.5f},"
                 f"{trans[0]:.2f},{trans[1]:.2f},"
                 f"{np.sqrt((resid ** 2).sum(1).mean()) * um_per_px:.2f},"
                 f"{np.median(dev):.2f},{np.percentile(dev, 95):.2f},"
                 f"{dev.max():.2f},{sx:.5f},{sy:.5f},{aniso:.4f},"
                 f"{shear:.4f},{aff_rms * um_per_px:.2f},"
                 f"{np.nanmedian(shape[:, 0]):.4f},"
                 f"{np.nanpercentile(shape[:, 0], 95):.4f},"
                 f"{np.nanmedian(shape[:, 1]):.4f}\n")

    overlay_figure(img, pred_all, obs, ii, jj, resid, scale, rot, trans,
                   um_per_px, os.path.basename(base), base + "-overlay.png")
    print("written:", base + "-overlay.png,", base + "-deviations.csv,",
          base + "-summary.csv")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(*sys.argv[1:3])
