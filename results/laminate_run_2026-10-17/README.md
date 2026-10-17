# Initial laminate sweep — 17 October 2026

Copied from the supplied `results.zip` archive. The archive contains the sweep CSV, held-out predictions, metrics, LaTeX table, and fitted model artifact.

## Run settings recorded in the CSV

- 50 design cases; length, root width, and taper ratio varied.
- 100 N transverse tip load; static linear analysis.
- 20 plies and 3.50 mm total thickness.
- 40 × 6 structured SHELL181 mesh (240 shell elements per case).
- 13 held-out cases and 5 analytical calibration points from training data.

## Review notes

The observed absolute deflection-to-length ratio ranges from 0.54 to 1.77. These displacements are outside the small-deflection assumption, so interpret the values as linearized compliance results rather than valid small-deflection displacement predictions at 100 N. A lower reference load and mesh-convergence study should be run before final reporting.

The archive does not include the MAPDL release, PyMAPDL version, or solver log. Add those to the project record when available.

Raw GP held-out MAPE against the laminate analytical benchmark is 0.139%; the corrected GP MAPE is 0.150%. For this run, sparse correction did not improve the result.
