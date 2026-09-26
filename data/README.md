# Data goes here (or in the repository root)

The course files are **not** in this repository. They are supplied by the
MT035A teaching team and are theirs to distribute, so they are listed in
`.gitignore`.

Put these in the repository root:

| file | what it is |
|---|---|
| `pattern_image.png` | the programmed reference pattern, as supplied |
| `reference_spots.csv` | the Week-1 spot table extracted from it |
| `bse-snapshot-*.png` | the flat-plate BSE snapshots |

Then `./run_all.sh` picks them up automatically.

You do **not** need any of these to try the code: `python3 selftest.py`
generates a synthetic pattern and a synthetic BSE image and runs the whole
pipeline on them.
