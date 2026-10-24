# MMA3001 Landing Gear Optimisation

This repository demonstrates limited-observation calibration of FEA-based composite-beam compliance models. The current experiment uses large-deflection static MAPDL analyses: 80 Dataset A training cases and 100 Dataset B target cases across five 20-case sets. It compares a mechanics-informed regression, quadratic response surface and Gaussian process, calibrating each model with two observations per set and validating on the remaining 18 cases. Three seeded calibration repeats are run per set.

The benchmark is a straight cantilever, not a curved landing gear. It uses no physical test data and does not establish real gear accuracy, impact performance or energy absorption.

## Current experiment files

- [Run guide](docs/extended_transfer_experiment.md): solver definition, set designs, calibration/validation, assumptions and commands.
- [Project report PDF](latex/project_report.pdf) and [LaTeX source](latex/project_report.tex).
- [Test/build record](docs/test_report.md).
- [API documentation](docs/api/index.html).

The FEA is a 0.25 m beam, fixed root, 10 N total downward tip-edge load, SHELL181 layered shell model, 40 x 6 mesh, three through-thickness integration points per ply, and `NLGEOM,ON` with automatic substepping and full Newton-Raphson. Compliance is load-dependent for this geometrically nonlinear analysis; the example predictions apply to the configured 10 N load.

## Set 5 and resin-ratio assumption

The attached Tab 5 change log adds `material_extrapolation_v2`: a twill/glass target with carbon/glass stiffness factors (2.25, 0.75) and a 2x resin-to-fibre volume-ratio multiplier. The runner implements this as a documented stiffness-only proxy: it starts from 50 vol% resin, doubles resin:fibre ratio to obtain 66.7 vol% resin, and remixes effective ply moduli with an illustrative epoxy matrix. Density and ply thickness remain fixed. The model is not told the changed ratio and calibrates only carbon/glass multipliers. This is a compliance-transfer stress test, not a physical resin-content or mass model.

## Results status

The published numeric results in [the dated results report](docs/extended_transfer_results_report_2026-10-22.md) are from the earlier linear-static study with four target sets. They do not include Set 5 and are not validation of the large-deflection five-set workflow. The current sweep has not yet been run. New output is written to `results/property_transfer_large_deflection_5sets/`; it cannot resume from earlier linear or four-set outputs.

## Run the latest sweep

Use the current `main` branch and run these PowerShell commands at the repository root:

```powershell
python -m unittest discover -s tests -v
python src/run_extended_transfer_study.py --stage design
python src/run_extended_transfer_study.py --stage dataset-a --max-cases 1
python src/run_extended_transfer_study.py --stage targets --target-set material_extrapolation_v2 --max-cases 1
```

The last two commands smoke-test MAPDL startup, including the new resin-ratio target. Then run the full 180-case experiment and model fitting; the runner will resume the smoke-tested case from the same output folder:

```powershell
python src/run_extended_transfer_study.py --stage all
```

A licensed local MAPDL installation and the Python dependencies are required. Check nonlinear convergence, reaction balance and mesh sensitivity on the local solver before interpreting results. Full instructions are in the run guide.
