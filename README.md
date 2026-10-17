# MMA3001 Landing Gear Optimisation

This repository develops a reproducible method for fitting a response model to finite-element results, correcting it with a small number of known reference observations, and assessing interpolation on held-out cases. The first benchmark is a straight, linearly width-tapered composite cantilever. It is a methods benchmark, not a prediction of a curved physical landing-gear leg.

## Current status

The initial composite sweep and surrogate fit are stored in [`results/laminate_run_2026-10-17/`](results/laminate_run_2026-10-17/).

- 50 design points; length, root width, and taper ratio varied.
- Fixed 20-ply `[CF-twill-0/GF-0]_{10}` layup; total cured thickness 3.50 mm.
- SHELL181 structured mesh: 40 elements along the beam and 6 across its width (240 elements per case); 3 section integration points per ply.
- Static linear analysis at 100 N. The recorded `|deflection|/length` spans 0.54–1.77, so these displacements violate the small-deflection assumption. Treat the results as a linearized compliance study only; do not interpret the 100 N deflections as physically valid small-deflection responses.
- Holdout set: 13 cases; 5 analytical calibration cases selected from training data. Raw GP MAPE against the analytical benchmark is 0.139%; corrected GP MAPE is 0.150%. The correction did not improve this run because the raw FEA surrogate already tracks the analytical benchmark closely.
- This run contains only a 40 × 6 mesh. Mesh convergence and the MAPDL release / solver log still need to be recorded.

The above scores demonstrate the pipeline on this idealized benchmark. They are not physical validation of the material card or the actual landing gear. The material constants are representative placeholders, and no physical test data are included.

## Repository layout

```text
config.example.json           Example SI-unit geometry and sweep configuration
layup.json                    20-ply material sequence and representative properties
src/                          FEA, laminate, analytical, fitting, and plotting code
tests/                        Lightweight checks of laminate setup and stiffness
latex/                        Report-ready method, layup, and result tables
results/laminate_run_.../     User-provided first sweep and fitted-model outputs
results/synthetic_demo/       Pipeline demonstration data; explicitly not FEA
docs/                         Requirement traceability, tests, generated API docs
latex/project_report.tex      Report draft with methods, results, and limitations
```

## Setup

Use Python 3.10–3.13 and install a licensed local Ansys MAPDL release. Install the Python dependencies in a virtual environment. PyMAPDL is installed separately from the project requirements.

### Windows PowerShell

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install ansys-mapdl-core
Copy-Item config.example.json config.json
```

If PyMAPDL cannot locate MAPDL, set `PYMAPDL_MAPDL_EXEC` to your installed MAPDL executable before launching the sweep. See the [PyMAPDL launch documentation](https://mapdl.docs.pyansys.com/version/stable/getting_started/launcher.html).

## Run the composite sweep

Edit `config.json` to set the project’s baseline length, root width, load, and design ranges. Values use SI units. The 20-ply layup and material constants are in `layup.json`.

`baseline.thickness_m`, `youngs_modulus_Pa`, and `poissons_ratio` belong to the legacy isotropic BEAM188 reference. They do not set the SHELL181 laminate properties. Laminate thickness is calculated from `layup.json` and is currently fixed at 3.50 mm.

Run the 50-case sweep:

```powershell
python src/run_laminate_shell_sweep.py `
  --config config.json `
  --layup layup.json `
  --output results/laminate_sweep.csv `
  --nx 40 `
  --ny 6
```

`--nx` and `--ny` are mesh divisions along the beam and across its width. Check mesh convergence by repeating the sweep with finer settings, for example `--nx 60 --ny 8`, and comparing compliance. Record the Ansys/MAPDL release, PyMAPDL version, and solver output alongside the results.

The starter configuration uses a 100 N load, which produces large linearized deflections for the current example geometry. For a small-deflection compliance benchmark, reduce the load (for example to 0.5 N) and confirm that normalized compliance remains consistent. Update the sweep output path so each run is preserved.

## Fit and assess the model

```powershell
python src/fit_corrected_surrogate.py `
  --input results/laminate_sweep.csv `
  --config config.json `
  --layup layup.json `
  --output-dir results/model
```

This fits a quadratic baseline and Gaussian-process surrogate to the FEA sweep, samples up to five analytical calibration points from training cases, and evaluates raw and corrected models on held-out cases. This tests interpolation. Keep the test set out of both fitting steps. Extrapolation requires a separate spatially held-out design region.

## Plot the results

```powershell
python src/plot_results.py `
  --sweep results/laminate_run_2026-10-17/laminate_sweep.csv `
  --predictions results/laminate_run_2026-10-17/model/test_predictions.csv `
  --output-dir results/laminate_run_2026-10-17/figures
```

This writes a PNG summary with compliance/design trends, the mass–compliance trade-off, held-out surrogate agreement, and deflection-to-length ratios.

## Test the non-solver code

```powershell
pip install -r requirements-dev.txt
python -m pytest
```

The test suite checks the configured ply count, total thickness, and positive effective stiffness. It does not launch MAPDL.

## Project documentation

- [Brief traceability](docs/requirements_traceability.md) lists the required evidence and remaining submission work.
- [Test report](docs/test_report.md) separates automated checks from preliminary model comparisons and FEA limitations.
- [LaTeX report draft](latex/project_report.tex) documents the problem, inputs, method, validation, limitations, and AI disclosure. Compile from `latex/` with LaTeX installed; resolve the marked approval and student-reflection items before submission.
- Generate HTML API docs with `python docs/generate_api_docs.py`; the output is in `docs/api/`.
- [MIT license](LICENSE) covers original project code and documentation. Third-party datasets, Ansys software, and solver outputs remain under their own terms.

## Assumptions and limitations

- Straight cantilever, clamped root, free tip, transverse tip load; linear elastic, small-deflection formulation.
- No geometric nonlinearity, contact, damage, failure, impact, or dynamic response.
- Representative orthotropic constants are used to demonstrate the workflow, not to represent a measured composite system.
- The alternating stack is unsymmetric and can exhibit extension–bending coupling.
- Analytical correction uses a laminate CLT A/B/D reduction with zero membrane resultants and an approximate transverse shear stiffness. It is benchmark truth for this simplified study, not a physical test result.
- The current sweep varies geometry only. Layup, ply count, thickness, material card, loading, and boundary conditions are fixed.

## Synthetic pipeline example

`results/synthetic_demo/` and `results/synthetic_proxy.csv` exercise the fitting scripts without Ansys. They use a deliberately coarse synthetic proxy and must not be presented as FEA results.
