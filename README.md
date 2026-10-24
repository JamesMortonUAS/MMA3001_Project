# MMA3001 Landing Gear Optimisation

This repository demonstrates a limited-observation calibration workflow using PyAnsys MAPDL finite-element data. The current runner is configured for a large-deflection static formulation, with 80 Dataset A training cases and 80 target cases. It compares mechanics-informed compliance regression, a quadratic response surface, and a Gaussian process; each target set uses three two-observation calibration repeats followed by held-out validation.

The study uses a straight composite cantilever as a method benchmark. It contains no physical test data and does not model the curved landing-gear leg, so its results demonstrate transfer within the selected FEA setup rather than physical gear accuracy.

## Current methodology and run instructions

- [Current experiment guide](docs/extended_transfer_experiment.md): design, FEA settings, accepted input domain, calibration/validation method, runtime fields, and commands.
- [Current project report PDF](latex/project_report.pdf) and [LaTeX source](latex/project_report.tex).
- [Test/build report](docs/test_report.md).
- [API documentation](docs/api/index.html).

The model uses a 0.25 m beam, fixed 10 N transverse tip-edge load, fixed root, SHELL181 layered shell elements, 40 x 6 mesh, and three through-thickness integration points per ply. The new runner enables MAPDL `NLGEOM,ON`, automatic time stepping, full Newton-Raphson, and configurable substep limits. In geometric nonlinear analysis, compliance depends on load; the example is valid only for the configured 10 N without regenerating and refitting.

## Results status

The checked-in numerical results in [the dated report](docs/extended_transfer_results_report_2026-10-22.md) were produced by the **previous linear-static formulation**. They are historical baseline results, not validation of the newly configured large-deflection solver. The nonlinear 160-case sweep has not yet been run. Do not use the old CSVs as nonlinear output; the new default output directory is `results/property_transfer_large_deflection/`, and resume validation checks that case definitions match.

The updated source and focused tests are on branch `feat/large-deflection-transfer-study`. After updating your checkout, run from the repository root in PowerShell:

```powershell
python -m unittest discover -s tests -v
python src/run_extended_transfer_study.py --stage design
python src/run_extended_transfer_study.py --stage dataset-a --max-cases 1
python src/run_extended_transfer_study.py --stage targets --max-cases 1
```

The one-case commands check MAPDL startup. Then run the complete FEA and model workflow:

```powershell
python src/run_extended_transfer_study.py --stage all
```

A licensed local MAPDL installation and the Python dependencies in `requirements.txt` are required. Full-run timing, nonlinear convergence, reaction balance, and mesh verification must be assessed from that installation before interpreting new prediction metrics.

## Earlier project work

Earlier analytical-beam, mesh-convergence, and hidden-property studies remain available in the repository history and result folders. Their assumptions and metrics differ and should be cited separately.
