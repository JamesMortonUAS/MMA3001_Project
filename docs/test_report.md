# Test report

## Automated checks

The repository includes a pytest check in `tests/test_laminate_section.py` for
the configured 20-ply count, 3.50 mm total laminate thickness, and positive
effective laminate stiffness. The GitHub Actions workflow in
`.github/workflows/tests.yml` installs the Python requirements and executes
`python -m pytest` on pushes and pull requests.

**Local execution status:** passed: 3 tests in 0.06 s with Python 3.12 and pytest 9.1.1.
Command executed: `python -m pytest -q`. The GitHub Actions workflow runs the same test suite on pushes and pull requests.

Reproduce with:

```bash
python -m pip install -r requirements.txt -r requirements-dev.txt
python -m pytest -q
```

These checks do not launch MAPDL and do not validate the solver model, mesh
convergence, boundary conditions, or material constants.

## Preliminary model assessment

The archived 50-case run includes 37 training cases and 13 held-out cases. Five
analytical calibration cases were selected from the training partition only.
For the held-out cases, the quadratic response surface achieved 0.977% MAPE
against the FEA values. The raw Gaussian process achieved 0.139% MAPE against
the analytical beam reference; the GP with sparse analytical discrepancy
correction achieved 0.150%. The correction therefore did not improve this
particular result. See `results/laminate_run_2026-10-17/model/metrics.json` for
the full metrics and `test_predictions.csv` for case-level predictions.

The comparisons are preliminary. The data describe a linear static SHELL181
analysis at 100 N, with absolute tip deflection divided by beam length between
0.54 and 1.77. This is outside the small-deflection domain of the analytical
benchmark. The archive also lacks mesh-refinement evidence and solver/version
metadata. Do not use the current errors as a claim of physical validation or
validated engineering accuracy.

## Required follow-up before submission

1. Run the automated tests and replace the pending status above with actual
   output.
2. Run a lower-load sweep that keeps the response within the analytical model's
   assumptions; check that compliance is load-independent for the linear model.
3. Compare at least one coarser and one finer mesh and record compliance
   changes, MAPDL/PyMAPDL versions, and solver logs.
4. Re-run the fitting workflow on the accepted dataset and update the report,
   figures, and test evidence from the regenerated artifacts.
