# MMA3001 Landing Gear Optimisation

This repository develops a reproducible FEA-to-model workflow for estimating the response of a composite beam when observations are limited. The current study trains a mechanics-informed response model on one PyAnsys MAPDL dataset, calibrates one shared effective-stiffness scale from one response in a second dataset, then evaluates predictions on held-out interpolation and extrapolation cases.

The benchmark is a straight, prismatic, 20-ply composite cantilever. It is a method demonstration, not a model of the curved landing-gear leg. No physical test data are used, and the new hidden-property sweep has not yet been run. The report therefore documents the current method and planned evaluation; it does not claim a calibration accuracy result.

## Finite-element model

- Flat layered-shell mesh using SHELL181; this is a 2D shell discretization with explicit layered section data, not a 3D solid model.
- Beam length 0.25 m, constant width per case, 40 elements along the beam and 6 across the width.
- Twenty alternating zero-degree plies: ten 200 gsm carbon 2x2 twill and ten 100 gsm UD E-glass/epoxy. Nominal cured ply thicknesses are 0.20 mm and 0.15 mm, for 3.50 mm total.
- Three through-thickness integration points per ply. Representative orthotropic material cards are in `layup.json`; values are assumptions rather than measured project material properties.
- Linear static analysis; root edge fixed in all degrees of freedom. Total 10 N global negative-Z force distributed evenly over the free tip edge. Output is mean tip-edge vertical displacement and compliance `abs(delta)/P`.

The default mesh has not been shown to be converged for the hidden-property study. Record MAPDL/PyMAPDL versions and solver logs, and check mesh sensitivity, shell orientation, load/support implementation, and `abs(delta)/L` before interpreting results.

## Dataset design

**Dataset A** is a 4D, three-level equally spaced tensor grid over width (0.03, 0.05, 0.07 m), thickness scale (0.7, 1.0, 1.3), in-plane stiffness scale (0.7, 0.9, 1.1), and transverse-shear stiffness scale (0.7, 1.15, 1.6). It would contain 81 cases. The target property pair `(0.9, 1.15)` is removed for all nine width/thickness combinations, leaving **72 FEA cases for training**.

**Dataset B** uses that withheld property pair and a 5x5 grid: widths 0.025, 0.0375, 0.05, 0.0625, 0.075 m and thickness scales 0.6, 0.8, 1.0, 1.2, 1.4. This gives **25 FEA cases**. The central case (ID 12; width 0.05 m, thickness scale 1.0) is used for the single calibration observation. Of the remaining 24, 8 lie inside both Dataset A geometry bounds (interpolation) and 16 lie beyond at least one bound (extrapolation). Here “withheld laminate” means a withheld pair of effective property scales applied to the same fixed stack; it is not a new stacking sequence.

The configured complete run requires 97 MAPDL analyses (72 in A plus 25 in B). A regular tensor grid gives transparent, deterministic factor coverage and makes withholding the whole property pair explicit. Three levels per axis limit cost but provide only modest resolution of curvature and interactions.

## Model, calibration, and evaluation

The code computes laminate effective bending rigidity per unit width `D_eff` from the full Classical Lamination Theory ABD matrix and transverse-shear rigidity `S_eff = sum(G_xz,k * t_k)`. Width and the shear correction factor give beam rigidities `EI = b * D_eff` and `kGA = kappa * b * S_eff`, where `kappa = 5/6`. For each Dataset A case the regression inputs are `L^3/(3 EI)` and `L/(kGA)`. The fitted model is a non-negative, no-intercept least-squares combination of those bending and shear compliance terms. Five-fold cross-validation on A checks within-grid predictive performance.

One Dataset B compliance observation estimates one bounded common stiffness multiplier (0.65 to 1.35) applied to both in-plane and transverse-shear moduli. The code selects the central case first. The correction minimizes squared relative residual at that observation. A single scalar response cannot identify independent bending and shear corrections, so the fitted value is a response-equivalent correction, not recovery of both hidden scales.

The calibration case is excluded from validation. Nominal and calibrated predictions are compared against the remaining Dataset B FEA values, separately for interpolation and extrapolation, using RMSE, MAE, MAPE, maximum absolute percentage error, and R-squared where defined. Dataset B hidden scales are retained for post-fit audit, but not used by calibration. Since both datasets use the same FEA setup, this validates the workflow under controlled simulation conditions; it does not establish physical validity or transfer to curved gear geometry.

## Run the study

Requires Python project dependencies, PyMAPDL, and a licensed MAPDL installation. From the repository root in PowerShell:

```powershell
python src/run_property_transfer_study.py --stage all `
  --config config.hidden_property.example.json `
  --layup layup.json
```

Stages can be run separately with `--stage dataset-a`, `--stage dataset-b`, or `--stage calibrate`. Calibration uses saved CSVs and does not launch MAPDL. Outputs are written to `results/property_transfer_grid/`; the model, predictions, and metrics go under its `model/` directory. Do not report new calibration or validation metrics until this grid run is completed and its solver provenance recorded.

See [`latex/project_report.pdf`](latex/project_report.pdf) for the full method and current status, [`docs/property_transfer_study.md`](docs/property_transfer_study.md) for the run protocol, and [`docs/requirements_traceability.md`](docs/requirements_traceability.md) for assignment requirements and remaining work.

## Historical and complementary work

Earlier quadratic-surrogate, analytical-reference, mesh-convergence, and calibration studies are preserved under `results/`. They used different configurations and are not results from the current 97-case hidden-property grid. The prior analytical and mesh work can inform separate method checks, but its performance metrics must not be substituted for the pending Dataset A/Dataset B evaluation.
