# Test report

## Automated checks

The pytest suite checks the 20-ply configuration, nominal 3.50 mm thickness,
and positive effective bending and shear stiffness. The suite does not launch
MAPDL and therefore does not verify solver setup or mesh convergence.

**Local status:** 6 tests passed with `python -m pytest -q` after implementing
the current quadratic workflow.

The non-solver fitting workflow was also exercised with a generated 120-row
smoke-test CSV. It completed the train/test split, quadratic fit, direct
least-squares stiffness-scale calibration, and fixed-section geometry-transfer
evaluation. It wrote metrics, predictions, model, and LaTeX table. The
smoke-test values were synthetic and are not FEA or project results.

## Current FEA status

The new 10 N, 120-case PyAnsys sweep has not been run in this environment.
There are no new FEA performance metrics to report. The archived 100 N run is
historical exploratory data and is not used as evidence for the new sweep.

## Required checks after the local MAPDL run

1. Record MAPDL and PyMAPDL versions, solver logs, and the exact run command.
2. Check that all configurations remain in the small-deflection range and
   verify load proportionality on representative cases.
3. Compare representative cases at coarse, configured, and finer meshes.
4. Run the fitting command on the resulting CSV and review held-out metrics,
   calibration factor, residual plots, and any outlier cases.
