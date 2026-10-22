# Methodology review and recommendations

## Purpose and scope

The project is a proof of concept for a future workflow in which a response model trained on inexpensive simulations is corrected using a small number of observations of a new material or structure. The present evidence is entirely FEA-to-FEA. There are no physical test observations in the study, and the finite-element model has not been validated against a physical landing gear.

The benchmark is a straight composite cantilever under a transverse tip load. It is useful for controlled testing of model fitting and calibration, but it does not establish landing-impact performance, energy absorption, failure, damage, fatigue, or transfer to the curved undercarriage geometry.

## Current completed methodology

### Finite-element setup

The archived sweep uses a flat layered-shell model in PyMAPDL/MAPDL:

- SHELL181 layered-shell elements on a structured 2D shell mesh.
- Beam length 0.25 m, 40 divisions along the length and 6 across the width.
- A 20-ply, alternating, all-zero-degree carbon twill and E-glass stack; nominal total thickness 3.5 mm.
- Three through-thickness integration points per ply.
- Root edge fixed in all degrees of freedom; a total 10 N vertical force is distributed across the free tip edge.
- Static linear-elastic analysis. The output is tip-edge mean vertical displacement and compliance, C = |delta|/P.

Material cards are representative modelling inputs, not measured properties of the physical landing gear. This run used one mesh; mesh convergence and independent model verification remain separate requirements.

### Datasets and partition

Dataset A comprised 72 simulations from a four-factor, three-level tensor grid over width, ply-thickness scale, in-plane stiffness scale and transverse-shear stiffness scale. The property pair (0.9, 1.15) was withheld across all width/thickness combinations, leaving 72 training cases from the possible 81.

Dataset B comprised 25 simulations at that withheld property pair on a 5 by 5 width/thickness grid. The central point (case 12) was used as the sole calibration observation. The remaining 24 cases were validation: eight lay within both Dataset A geometry bounds and 16 exceeded at least one bound.

### Current model and calibration

The training model is a two-term mechanics-informed, nonnegative, no-intercept least-squares regression:

\[
\widehat C = c_b\frac{L^3}{3EI}+c_s\frac{L}{kGA},
\]

where EI is the beam-equivalent laminate bending rigidity and kGA is the shear rigidity including the shear correction factor. The training fit uses all Dataset A FEA results; five-fold cross-validation checks its within-grid performance.

A single bounded common stiffness multiplier is fitted to the central Dataset B response, then the corrected model predicts the other 24 Dataset B cases. One observation can identify one scalar response correction; it cannot separately recover in-plane and transverse-shear material scales.

### Current result and interpretation

The completed 97-case run gave Dataset A cross-validation MAPE of 0.0571%. On the 24 held-out Dataset B cases, MAPE fell from 9.9447% before calibration to 0.05393% after calibration. The calibrated MAPE was 0.03369% on the eight geometry-interpolation points and 0.06405% on the 16 tested, bounded geometry-extrapolation points.

This is strong transfer within a controlled, shared FEA setup. It is not physical validation: training and target data came from the same FE formulation and material-card framework. Nor does the fitted common factor identify the actual hidden material properties. The softest old cases also had large normalized deflections, so linear small-deflection assumptions deserve attention.

## Recommended next experimental design

The next experiment should separate three questions: interpolation to a new stiffness combination, extrapolation to new stiffness values or a new layup, and extrapolation in geometry. Use a fixed length and load, and hold mesh, support, load application, solver settings and output extraction constant across all runs.

### Dataset A: training grid

Run a complete 5 x 4 x 4 factorial design, for 80 cases:

| Factor | Levels |
|---|---|
| Width | 0.02, 0.04, 0.06, 0.08, 0.10 m |
| Total ply count | 8, 16, 24, 32 |
| Carbon/glass stiffness scale pair | (1,1), (1,2), (2,1), (2,2) |

The half-stack repeats carbon twill (T), glass (G), then mirrors about the mid-plane. Use a cured thickness of 0.20 mm for each T ply and 0.10 mm for each G ply. Full sequences and nominal thicknesses are:

| Total plies | Full sequence, bottom to top | T/G counts | Thickness |
|---:|---|---:|---:|
| 8 | T-G-T-G-G-T-G-T | 4/4 | 1.2 mm |
| 16 | (T-G repeated 4 times) + its reverse | 8/8 | 2.4 mm |
| 24 | (T-G repeated 6 times) + its reverse | 12/12 | 3.6 mm |
| 32 | (T-G repeated 8 times) + its reverse | 16/16 | 4.8 mm |

For a target total ply count N, repeat T-G for N/2 plies in the half-stack and mirror the half-stack. This gives a symmetric stack.

Scale every elastic modulus and shear modulus of carbon plies with the carbon factor and every elastic/shear modulus of glass plies with the glass factor. Keep ply thickness, density and Poisson ratios fixed. Record both factors and calculate the laminate bending and transverse-shear section properties for each configuration. Flexural stiffness alone is insufficient if shear deformation is included.

### Target sets and calibration

Create 20 geometry cases per target set using the same 5 widths and 4 ply counts as Dataset A:

1. **Material interpolation:** carbon/glass factors (1.5, 1.5), inside the training factors' range.
2. **Material extrapolation:** factors (2.25, 0.75), outside at least one training factor bound.
3. **Layup transfer:** a symmetric half-stack made by repeating CF unidirectional (U), CF twill (T), and E-glass (G), truncated to half the specified ply count and mirrored. Use 0.20 mm per carbon ply and 0.10 mm per glass ply. This changes layup architecture, not only stiffness. Record the complete sequence, ply counts, total thickness, mass, bending rigidity and shear rigidity.
4. **Geometry extrapolation:** use widths [0.015, 0.04, 0.06, 0.08, 0.15] m and ply counts [8, 16, 32, 100]. The 100-ply alternating T/G symmetric stack is 15 mm thick. Label each case by whether width and ply count are within or outside Dataset A bounds; the 100-ply and 0.015/0.15 m extremes are extrapolations.

For each of the first three target sets, select two distinct response observations to calibrate two unknown material multipliers (carbon and glass) and validate against the other 18 cases. Repeat this with three seeded random calibration pairs. Use identical selected pairs for every model being compared. The two-observation approach is only meaningful if the two measurements provide independent sensitivity to the two fitted factors; record the calibration Jacobian rank/condition number and warn about poorly conditioned pairs.

The geometry extrapolation set should likewise use three repeated, seeded two-point calibration partitions, while keeping its target properties nominal. This makes calibration variability visible without additional FEA runs.

### Models and evaluation

Compare three approaches trained only on Dataset A:

1. The two-term mechanics-informed compliance regression using bending and shear predictors.
2. A quadratic response surface, with inputs scaled before fitting.
3. A Gaussian-process regressor with standardized inputs and a documented kernel/noise treatment.

Use the same FEA compliance target, training data, and target calibration observations for all methods. Select/tune each model using Dataset A only. A practical choice for the quadratic and Gaussian-process models is to fit log compliance, which enforces positive predictions after exponentiation. For target calibration, estimate the carbon and glass multipliers by nonlinear least squares against the two observations; validation responses must not enter that fit.

Report nominal errors before calibration and held-out errors after calibration for each set, model and random repeat. Include RMSE, MAE, MAPE, maximum absolute percentage error, R-squared where informative, calibration IDs, fitted factors, Jacobian rank/condition number, and prediction plots. Separate interpolation from extrapolation by each input axis rather than using one broad label.

### Important clarification: fit versus calibration

The mechanics-informed response model is estimated from the full training dataset. The two target observations do not train the response surface; they estimate the two unknown target material multipliers used to calculate the target's effective section properties. This keeps the exercise focused on transfer and calibration while retaining a physically interpretable response model.

## Run size and execution recommendations

The design requests 80 Dataset A cases and 80 unique target FEA cases (four target sets of 20 each), for 160 unique simulations. Three calibration repeats reuse those target outputs; they do not require extra FEA. Run one smoke-test case first, then all cases. Store per-case configuration and solver version, preserve MAPDL logs, and write incremental CSV results so a long sweep can resume or be audited.

Before interpreting model errors, run mesh checks at a soft and stiff configuration, verify shell local axes and ply order, inspect reaction-force balance, and check displacement divided by length. If the sweep uses a linear static model, describe it as a small-deflection compliance study. The legacy 10 N and 0.25 m settings may be retained initially for comparability, but adjust if they violate the intended small-deflection range.

## Recommendations and limits

- Preserve the current 97-case run as a reproducible baseline; do not replace its scripts or metrics with the new design.
- Treat the new results as a simulation-method validation, not evidence of physical landing-gear accuracy.
- Report full stacking sequences, actual thickness, mass, section properties and the exact definition of every material-property scale.
- Keep calibration and validation partitions strictly separate. Repeated random pairs expose dependence on calibration-point choice.
- Analyze identifiability before interpreting fitted carbon/glass scales. Two points are the algebraic minimum for two unknowns, not a guarantee of a stable estimate.
- Compare models on identical held-out cases and report uncertainty or variation across calibration repeats.
- Consider physical testing, curved-gear geometry, large-deflection contact, energy absorption, failure and manufacturability as later stages beyond this benchmark.

## Archived result locations

The original implementation is `src/run_property_transfer_study.py`; the associated data and metrics are under `results/property_transfer_grid/`. The project report is `latex/project_report.tex` with the rendered report under `latex/project_report.pdf`. The new proposed experiment will live on its own branch so the completed baseline remains unchanged.
