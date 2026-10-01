#!/usr/bin/env bash
# MT035A Week 2 - run the whole Week-2 analysis and collect the results.
#
#   ./run_all.sh                 # processes every bse-snapshot-*.png here
#   ./run_all.sh a.png b.png     # or the snapshots you name
#
# Produces, per snapshot: <name>-overlay.png, <name>-deviations.csv,
# <name>-summary.csv. Then merges them into week2-deviations.csv and
# week2-registration-summary.csv, and runs the orientation cross-check on the
# first two snapshots.
set -euo pipefail
cd "$(dirname "$0")"

if [ "$#" -gt 0 ]; then
    shots=("$@")
else
    shopt -s nullglob
    shots=()
    for f in bse-snapshot-*.png; do
        # skip the overlay figures this script produces
        case "$f" in *-overlay.png) continue ;; esac
        shots+=("$f")
    done
fi
if [ "${#shots[@]}" -eq 0 ]; then
    echo "no snapshots found (expected bse-snapshot-*.png)" >&2
    exit 1
fi

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
