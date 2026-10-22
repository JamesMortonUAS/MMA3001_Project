# Hidden-property FEA transfer study

## Purpose

This workflow tests whether a response model trained on a designed FEA grid can
be calibrated from one response observation and then predict unseen responses
for an excluded laminate property pair. It provides separate validation sets
for geometry interpolation and extrapolation. Dataset A and Dataset B are both
FEA, so this is a methodological FEA-to-FEA exercise, not physical validation.

## Modelled quantity and equivalent section properties

The benchmark is a straight, prismatic, 20-ply composite cantilever at fixed
length. A transverse tip-load resultant is distributed evenly along the free
tip edge; the output is mean vertical displacement on that edge. Compliance
`C = abs(delta) / P` is the model target.

For each FEA case, the script records total laminate thickness, mass, and
laminate bending and transverse-shear stiffnesses. The CLT effective bending
stiffness per unit width is `D_eff`; transverse-shear stiffness per unit width
is `S_eff = sum(G_xz,k * t_k)`. The full ABD matrix retains extension-bending
coupling. Beam rigidities are `EI = b * D_eff` and
`kGA = kappa * b * S_eff`.

The fitted response is
`C_hat = c_b * L^3/(3*EI) + c_s * L/(kGA)`, with `c_b,c_s >= 0`.
The coefficients are learned from Dataset A FEA.

## Grid and data partition

Dataset A is a four-dimensional tensor-product grid. Each input axis has three
evenly spaced levels, yielding `3^4 = 81` design points before exclusions:
width, ply-thickness scale, in-plane stiffness scale, and transverse-shear
stiffness scale. The configured target pair `(0.9, 1.15)` is excluded from
training, leaving **72 Dataset A FEA cases**. Here, "laminate" means an effective
property-scale pair applied to the supplied 20-ply sequence. The stacking
sequence itself does not change.

Dataset B uses the excluded target pair and a `5 x 5` evenly spaced grid in
width and thickness: **25 FEA cases**. The central case is selected as the
single calibration observation. The other 24 cases are validation candidates.
A case is labelled interpolation when both width and thickness lie inside the
Dataset A bounds. Otherwise it is labelled extrapolation. After using the
central case for calibration, the example configuration has eight interpolation
and sixteen extrapolation validation cases.

The one observation fits one **common stiffness correction scale**, applied to
both in-plane and transverse-shear stiffnesses. A single compliance value cannot
identify independent in-plane and shear scales. The fitted scale is therefore
a response-equivalent correction, not a separate recovery of both material
properties. Hidden Dataset B scales are saved only for post-fit audit; they are
not read by the calibration objective.

Five-fold cross-validation on Dataset A reports within-training-grid model
performance. Dataset B scores are reported separately for interpolation and
extrapolation cases, before and after one-point calibration.

## Run

From the repository root, with project requirements, PyMAPDL, and a licensed
local MAPDL installation:

```powershell
python src/run_property_transfer_study.py `
  --stage all `
  --config config.hidden_property.example.json `
  --layup layup.json
```

The default full run performs 97 FEA cases: 72 training cases plus 25 Dataset B
cases. Separate stages are `dataset-a`, `dataset-b`, and `calibrate`.
The calibration stage reads existing CSV outputs and does not launch MAPDL.
New outputs are written under `results/property_transfer_grid/`, preserving the
earlier Latin-hypercube results under `results/property_transfer/`.

## Scope and limitations

- The laminate stack is fixed. The excluded "laminates" are stiffness-scale
  combinations, not new ply angles, material proportions, or stacking orders.
- The common one-point correction does not identify independent in-plane and
  shear properties. Report prediction errors and the fitted common scale.
- Geometry extrapolation is defined using width and thickness bounds. It does
  not establish extrapolation to a curved landing-gear leg or changed support,
  length, load direction, or dynamic response.
- The 10 N static linear FEA outputs remain subject to mesh-convergence and
  small-deflection checks.
- Record MAPDL/PyMAPDL versions and solver logs when available.
