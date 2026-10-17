# Test report

## Automated checks

The pytest suite checks the 20-ply configuration, nominal 3.50 mm thickness,
and positive effective bending and shear stiffness. The suite does not launch
MAPDL and therefore does not verify solver setup or mesh convergence.

**Local status:** 7 tests passed with `python -m pytest -q`, including checks
that calibration sensitivity subsets are unique and come only from training
cases.

The non-solver fitting workflow was also exercised with a generated 120-row
smoke-test CSV. It completed the train/test split, quadratic fit, direct
least-squares stiffness-scale calibration, and fixed-section geometry-transfer
evaluation. It wrote metrics, predictions, model, and LaTeX table. The
smoke-test values were synthetic and are not FEA or project results.

## FEA and sensitivity-analysis status

A supplied archive contains a 10 N, 120-case dataset identified as PyMAPDL
SHELL181 FEA output. The dataset has 96 training cases and 24 held-out cases.
The MAPDL/PyMAPDL versions and solver log were not included, and the sweep was
not independently rerun in this environment.

The new `src/analyze_calibration_sensitivity.py` script was run on the supplied
sweep with calibration counts 1, 2, 3, and 5 and 30 unique subsets per count.
The quadratic model and 24-case geometry holdout remained fixed. The result
files are archived under `results/calibration_sensitivity_2026-10-17/`. Across
the random and centred/maximin subsets, the average held-out MAPE was worse
than the nominal model at every calibration count, though the estimated
stiffness scale became less variable as the number of calibration observations
increased. This is a sensitivity result for one train/test split, not a claim
about all possible samples or physical test data.

## Required checks after the local MAPDL run

1. Record MAPDL and PyMAPDL versions, solver logs, and the exact run command.
2. Check that all configurations remain in the small-deflection range and
   verify load proportionality on representative cases.
3. Compare representative cases at coarse, configured, and finer meshes.
4. Review the sensitivity outputs and, if needed, repeat with other fixed
   train/test seeds. Then review calibration factor distributions, held-out
   residuals, and outlier cases before interpreting model performance.
