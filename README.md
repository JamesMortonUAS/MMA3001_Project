# MMA3001 Landing Gear Optimisation

This repository demonstrates a reproducible workflow for fitting a response model to PyAnsys finite-element results, correcting that model with a small number of known reference observations, and checking interpolation on held-out cases. The benchmark is a straight, linearly width-tapered, 20-ply composite cantilever. It is intended to validate the modelling method before applying it to a curved leaf-spring landing gear.

## Current workflow

The current input is [`config.quadratic.example.json`](config.quadratic.example.json). It sets a 10 N transverse tip load and sweeps five factors: beam length, root width, taper ratio, total laminate thickness, and a shared effective laminate stiffness scale. For every design, the runner also records the derived laminate bending stiffness $D_{\mathrm{eff}}$ and transverse shear stiffness $S_{\mathrm{eff}}$.

The default 120-case Latin-hypercube design supports a full degree-two response surface in five inputs (21 coefficients) with cases remaining for a held-out evaluation. Thickness scaling changes every cured ply thickness and therefore total thickness and mass. The stiffness scale multiplies ply elastic and shear moduli while keeping Poisson ratios and density fixed; it provides a controlled range of effective section properties for this methods benchmark.

### Run the PyAnsys sweep

Install Python 3.10–3.13, the project requirements, `ansys-mapdl-core`, and a licensed local Ansys MAPDL release. For Windows PowerShell:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install ansys-mapdl-core
```

Run the 120-case layered-shell sweep:

```powershell
python src/run_laminate_shell_sweep.py `
  --config config.quadratic.example.json `
  --layup layup.json `
  --output results/quadratic_sweep.csv `
  --nx 40 `
  --ny 6
```

The model uses SHELL181, 20 explicit plies, three through-thickness integration points per ply, a clamped root, and an evenly distributed 10 N transverse force on the tip edge. `--nx` and `--ny` are mesh divisions along the beam and across its width. Before interpreting results, repeat representative cases with a finer mesh and record the Ansys/MAPDL release and PyMAPDL version.

### Fit and calibrate the quadratic model

```powershell
python src/fit_corrected_surrogate.py `
  --input results/quadratic_sweep.csv `
  --config config.quadratic.example.json `
  --output-dir results/quadratic_model
```

The model fits a full quadratic response surface to FEA compliance $C=|\delta|/P$. Calibration directly estimates one effective stiffness scale by bounded least squares against analytical observations. The default uses one calibration geometry, the mathematical minimum for one scalar parameter, and holds thickness at its nominal scale. The 20% test set is excluded from fitting and calibration. It tests predictions at unseen geometries while keeping thickness and effective stiffness fixed, which specifically tests geometry transfer. Set `model.calibration_points` higher to fit the same scalar using multiple analytical observations.

The fitted scale multiplies all ply elastic and shear moduli, so the corresponding laminate bending stiffness $D_{\mathrm{eff}}$ and transverse shear stiffness $S_{\mathrm{eff}}$ change together. This single-parameter model does not identify separate bending and shear stiffness corrections. The config sets the initial scale and the scale used to generate the analytical reference; both are 1.0 for the nominal benchmark. For future physical calibration, the reference scale is replaced by measured responses, not treated as a known physical value.

### Outputs and checks

The sweep CSV includes design variables, actual thickness, $D_{\mathrm{eff}}$, $S_{\mathrm{eff}}$, tip deflection, compliance, and mass. Fitting writes `metrics.json`, `test_predictions.csv`, `quadratic_model.joblib`, and a LaTeX results table. Plot held-out predictions with:

```powershell
python src/plot_results.py `
  --sweep results/quadratic_sweep.csv `
  --predictions results/quadratic_model/test_predictions.csv `
  --output-dir results/quadratic_model/figures
```

Run the solver-independent checks with `python -m pytest`. These checks do not launch MAPDL.

## Archived preliminary 10 N run

A 120-case 10 N SHELL181 sweep and one-point calibration output are preserved in [`results/quadratic_run_2026-10-17_one_point_calibration/`](results/quadratic_run_2026-10-17_one_point_calibration/). This was produced with the earlier compliance-multiplier fitter, before the revised direct effective-stiffness calibration. On 24 held-out cases, the raw model had 5.95% MAPE against the analytical reference; the one-point correction increased it to 7.27%. Treat this as a baseline only and rerun with the current scripts. The run's deflection-to-length ratio reached 0.330, so the highest-displacement cases may violate the benchmark's small-deflection assumption.

## Assumptions and limits

- Static linear elasticity and small-deflection beam theory are used. Confirm the 10 N response is within the linear small-deflection range for every design; the applied load alone does not guarantee that.
- The beam is straight, linearly width-tapered, clamped at the root, and loaded transversely at the tip. The model does not include geometric nonlinearity, contact, damage, failure, or impact dynamics.
- Material cards are representative numerical inputs for demonstrating the workflow, not measured properties of the project’s composite.
- The stack is unsymmetric, and the CLT effective bending calculation retains extension–bending coupling. The shear term uses the laminate through-thickness $G_{xz}t$ sum and the configured correction factor.
- The analytical equation is the ground truth for this benchmark only. It is not a physical test result and is not an exact solution for a curved landing-gear leg.
- The default DOE has five inputs. Quadratic interpolation is assessed within the sampled domain; extrapolation requires a separate, deliberately held-out region.

## Repository contents

`layup.json` defines the 20-ply material sequence and input properties. `src/` contains the PyMAPDL sweep, laminate section setup, analytical reference, quadratic fit, plotting, and mass optimization tools. `latex/` contains report material; `docs/requirements_traceability.md` tracks the brief. The archived `results/laminate_run_2026-10-17/` folder contains the earlier 100 N exploratory run and is superseded by the 10 N sweep configuration for the current workflow. The archived data are not the output of the new configuration.
