# MMA3001 Landing Gear Optimisation

This repository demonstrates a limited-observation model calibration workflow using PyAnsys MAPDL finite-element data. The current extended transfer experiment has completed 80 Dataset A simulations and 80 target simulations. It fits mechanics-informed compliance regression, quadratic response-surface and Gaussian-process models to Dataset A, calibrates two carbon/glass stiffness multipliers from two observations per target set, then evaluates held-out predictions.

The experiment is a straight-beam, FEA-to-FEA method study. It uses no physical test data and does not yet model the curved landing-gear leg. The results therefore demonstrate transfer within this assumed FEA setup; they do not validate physical landing-gear response.

## Latest results

- Full analysis: [extended transfer study results and interpretation](docs/extended_transfer_results_report_2026-10-22.md).
- Run design and commands: [extended transfer experiment guide](docs/extended_transfer_experiment.md).
- FEA inputs: 0.25 m cantilever; 10 N total tip-edge load; fully fixed root edge; layered SHELL181 shell mesh, 40 × 6, with 3 integration points per ply; linear-static analysis.
- Dataset A: 80 cases across 5 widths, 4 ply counts and 4 carbon/glass stiffness combinations.
- Target data: 80 cases across material interpolation, material extrapolation, alternate layup and geometry extrapolation sets.
- Validation: three calibration repeats per target set; two observations calibrate and 18 are held out each time.

Result CSVs and the summary JSON are in [the dated result folder](results/property_transfer_extended_run_2026-10-22/). The fitted joblib objects are not included; per-case predictions, metrics, case designs and FEA outputs are.

## Main findings and limits

The quadratic response surface had the lowest Dataset A five-fold CV MAPE (0.0053%), followed by the Gaussian process (0.0195%) and mechanics-informed regression (0.1083%). On the same-layup material-transfer sets, informative two-point calibration greatly reduced held-out error. One random material-extrapolation pair could not distinguish carbon from glass stiffness effects: both calibration cases had 32 plies and differed only in width, producing a rank-deficient mechanics calibration. This motivates selecting calibration cases based on sensitivity and checking Jacobian conditioning.

The linear-static runs include cases with large deflection relative to length: the maximum target `|delta|/L` is 1.70, and 21/80 targets exceed 0.1. Those cases reproduce the linear FEA result but do not support a small-deflection physical interpretation. The material cards are illustrative assumptions; the mesh and solver setup require verification. Dynamic landing, failure, impact-energy absorption, physical testing and curved-gear optimization remain outside the demonstrated scope.

## Run the extended study

Requires Python dependencies, PyMAPDL and a licensed local MAPDL installation. In PowerShell at the repository root:

```powershell
python src/run_extended_transfer_study.py --stage design
python src/run_extended_transfer_study.py --stage dataset-a
python src/run_extended_transfer_study.py --stage targets
python src/run_extended_transfer_study.py --stage fit
```

Use `--stage dataset-a --max-cases 1` or `--stage targets --max-cases 1` to smoke-test MAPDL startup first. See the experiment guide for resume behavior and output details.

## Earlier project work

Earlier analytical-benchmark, mesh-convergence and hidden-property studies remain in the repository history and prior results folders. Their configurations and metrics differ from this extended experiment and should be cited separately.
