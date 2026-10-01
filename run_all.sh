#!/usr/bin/env bash
# MT035A - reproduce every result in this repository.
#
#   ./run_all.sh                 # weeks 1-3 and all figures (a few minutes)
#   ./run_all.sh --sweep         # the above plus the Week-4 robustness sweep
#   ./run_all.sh a.png b.png     # weeks 1-2 on the snapshots you name
#
# The sweep is separated only because it re-runs the whole pipeline 61 times
# and takes several minutes; nothing else here does.
#
# Produces, per snapshot: <name>-overlay.png, <name>-deviations.csv,
# <name>-summary.csv. Then merges them into week2-deviations.csv and
# week2-registration-summary.csv, runs the orientation cross-check, the
# systematic/random split, the plate-geometry profile, the Week-3 tilt
# analysis, and regenerates every figure referenced by README.md.
set -euo pipefail
cd "$(dirname "$0")"

sweep=0
args=()
for a in "$@"; do
    case "$a" in
        --sweep) sweep=1 ;;
        *) args+=("$a") ;;
    esac
done

if [ "${#args[@]}" -gt 0 ]; then
    shots=("${args[@]}")
else
    shopt -s nullglob
    shots=()
    for f in bse-snapshot-*.png; do
        case "$f" in *-overlay.png) continue ;; esac
        shots+=("$f")
    done
fi
if [ "${#shots[@]}" -eq 0 ]; then
    echo "no snapshots found (expected bse-snapshot-*.png)" >&2
    exit 1
fi

# --- 0. validation first, so a broken pipeline fails before it produces -----
#        numbers that look plausible.
echo "=== self-test against synthetic ground truth ==="
python3 selftest.py
echo

for shot in "${shots[@]}"; do
    echo "=== $shot ==="
    python3 bse_register.py "$shot"
    echo
done

# merge the per-spot deviations, prefixed with the image name
{
    echo "image,$(head -n 1 "${shots[0]%.*}-deviations.csv")"
    for shot in "${shots[@]}"; do
        base="${shot%.*}"
        tail -n +2 "${base}-deviations.csv" | sed "s|^|${base},|"
    done
} > week2-deviations.csv

# merge the one-line summaries
{
    head -n 1 "${shots[0]%.*}-summary.csv"
    for shot in "${shots[@]}"; do
        tail -n +2 "${shot%.*}-summary.csv"
    done
} > week2-registration-summary.csv

echo "=== merged ==="
if command -v column >/dev/null 2>&1; then
    column -s, -t < week2-registration-summary.csv
else
    cat week2-registration-summary.csv
fi
echo
echo "week2-deviations.csv: $(($(wc -l < week2-deviations.csv) - 1)) rows"

if [ "${#shots[@]}" -ge 2 ]; then
    echo
    echo "=== orientation cross-check ==="
    python3 check_orientation.py "${shots[0]}" "${shots[1]}"

    echo
    echo "=== systematic vs random, over every image ==="
    devs=()
    for shot in "${shots[@]}"; do devs+=("${shot%.*}-deviations.csv"); done
    python3 analyse_repeatability.py "${devs[@]}"

    echo
    echo "=== spot geometry across the plate ==="
    python3 focus_profile.py "${shots[@]}"
fi

# --- 3. Week 3: tilt, only if the tilted plates are present ----------------
if ls TILT_*.png >/dev/null 2>&1; then
    echo
    echo "=== Week 3: focus gradient on tilted plates ==="
    flats=()
    for shot in "${shots[@]}"; do flats+=("$shot"); done
    [ -f img3.png ] && flats+=(img3.png)
    python3 tilt_analysis.py --flat "${flats[@]}" --test TILT_*.png
else
    echo
    echo "(skipping Week 3: no TILT_*.png in this directory)"
fi

# --- 4. Week 4: the robustness sweep --------------------------------------
if [ "$sweep" -eq 1 ]; then
    echo
    echo "=== Week 4: robustness sweep (several minutes) ==="
    python3 robustness.py "${shots[0]}" reference_spots.csv
else
    echo
    echo "(skipping Week 4: pass --sweep to run it, or"
    echo " 'python3 robustness.py --plot-only' to redraw from week4-robustness.csv)"
fi

# --- figures --------------------------------------------------------------
echo
echo "=== figures ==="
for fig in make_reference_figure.py make_workflow_figure.py \
           make_registration_overlay.py make_slide_overlay.py; do
    [ -f "$fig" ] || continue
    echo "--- $fig"
    python3 "$fig" || echo "    (skipped: $fig needs inputs not present here)"
done

echo
echo "done. Figures referenced by README.md are in docs/."
