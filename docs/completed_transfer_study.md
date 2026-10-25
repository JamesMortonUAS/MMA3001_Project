# Completed transfer experiment

## Scope and evidence

The submission is based on `results/property_transfer_extended/` from the completed extended sweep. It contains 80 Dataset A cases and 80 Dataset B cases. The report and detailed dated results document use only this sweep and its model evaluation.

## Finite element definition

- Straight rectangular composite cantilever, length 0.25 m; global X along its length.
- Fully fixed root edge; 10 N total force in global -Z divided equally across free-tip-edge nodes.
- Linear-static orthotropic elastic analysis.
- Layered SHELL181 discretisation: 40 divisions along length, 6 across width, 3 integration points per ply.
- Response: absolute mean tip-edge vertical displacement divided by 10 N, in m/N.

## Training design: 80 cases

Full factorial grid: widths 0.02, 0.04, 0.06, 0.08, 0.10 m; ply counts 8, 16, 24, 32; independent carbon/glass modulus multipliers (1,1), (1,2), (2,1), (2,2). The symmetric half-stack repeats twill/glass and is mirrored. Carbon twill is 0.20 mm per ply and glass is 0.10 mm. Total thicknesses are 1.2, 2.4, 3.6, 4.8 mm.

## Target design: four sets of 20 cases

1. Material interpolation: base layup and training geometry levels, factors (1.5,1.5).
2. Material extrapolation: base layup and training geometry levels, factors (2.25,0.75).
3. Layup transfer: repeating UD/twill/glass half-stack, truncated and mirrored, training widths and ply counts, factors (1,1). UD and twill plies are 0.20 mm; glass is 0.10 mm. Thicknesses are 1.4, 2.8, 4.0, 5.4 mm.
4. Geometry extrapolation: base layup, widths 0.015, 0.04, 0.06, 0.08, 0.15 m; ply counts 8, 16, 32, 100; factors (1,1).

## Model fitting and evaluation

Dataset A alone trains nonnegative mechanics-informed compliance regression, a quadratic response surface in log predictors/output, and a Gaussian process. Five-fold cross-validation estimates within-training-domain performance.

For each target set, three seeded pairs of two cases supply compliance observations. Carbon and glass stiffness multipliers are estimated by least squares in log compliance; calculated section properties are recomputed while fitted surrogate coefficients remain fixed. The other 18 cases per repeat are withheld. Calibration pairs are shared across models. Nominal and calibrated predictions are compared using MAPE, RMSE, maximum error, and R-squared. Jacobian conditioning diagnoses ambiguous parameter estimates.

Nominal predictions require configured material cards. The alternate-layup set therefore tests transfer with known nominal cards, rather than blind identification of all material properties.

## Archived outputs and reproduction boundary

The archived CSV/JSON files record designs, responses, predictions, calibration splits, fitted factors, and metrics. Use these files to audit the report tables. Later source development is retained in the repository; the latest runner configuration is not the definition of this completed experiment. Do not mix outputs across configurations or describe a new execution as a reproduction without matching the original solver definition and case design.

## Limitations and future work

Large displacement-to-length ratios limit physical small-deflection interpretation of the computed linear responses. Mesh convergence is not established for every transfer case. Future studies should reduce load or restrict geometry, choose informative calibration pairs, and assess measurement noise. Large-deflection modelling and physical gear validation are possible subsequent extensions.
