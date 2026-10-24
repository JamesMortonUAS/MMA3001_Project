# Extended transfer study: results and interpretation

> **Historical baseline notice (24 October 2026):** All numerical results in this report come from the earlier linear-static analysis. The runner has since been changed to large-deflection static analysis (`NLGEOM,ON`), but that 180-case sweep with five target sets has not yet been run. These scores must not be presented as validation of the updated formulation. See the [current experiment guide](extended_transfer_experiment.md) and [test/build record](test_report.md) for the new setup and status.

**Run analysed:** `results/property_transfer_extended/` from the supplied `results(5).zip` archive. 
The archive also contains an earlier `property_transfer_grid` run; it is not combined with the extended-study results below.

## Executive summary

The 80-case Dataset A response models fit the training regime extremely closely. Five-fold cross-validation favoured the quadratic response surface, followed by the Gaussian process and then the mechanics-informed regression. On the two same-layup material-transfer sets, two-point calibration usually recovered the hidden stiffness changes and reduced held-out error from roughly 50–99% to well below 1%. One material-extrapolation calibration pair failed to identify carbon and glass effects separately; this is a key design result, not a result to average away.

The new-layup and geometry-extrapolation sets were already predicted accurately before calibration because their effective section properties were close to the Dataset A property basis. Calibration offered little benefit and sometimes worsened validation error. Overall this is a successful **FEA-to-FEA workflow demonstration**, not evidence that the FEA or the material cards represent a physical landing gear.

A critical validity issue is that many cases have large tip displacement relative to beam length despite the linear-static setup: 22/80 training cases and 21/80 target cases have `|delta|/L > 0.1`; the largest target ratio is 1.70. Predictions still reproduce the computed linear FEA responses, but these cases are outside a small-deflection interpretation and should be rerun at a lower load or excluded from a small-deflection study.

## What was run

- Dataset A: 80 FEA cases on a 5-width × 4-ply-count × 4-carbon/glass-stiffness-pair grid.
- Targets: 80 FEA cases, 20 each in material interpolation, material extrapolation, layup transfer, and geometry extrapolation sets.
- FEA: 0.25 m straight cantilever; 10 N total distributed over the free tip edge; root edge fully fixed; layered SHELL181 mesh, 40 × 6, 3 through-thickness integration points per ply; linear-static analysis.
- Response: tip-edge mean vertical deflection and compliance `C = |delta|/P`.
- Models trained only on Dataset A: nonnegative mechanics-informed bending/shear regression, quadratic response surface on log inputs/output, and Gaussian process on the same log inputs/output. Five-fold CV was also run on Dataset A.
- Calibration/evaluation: three random repeats per target set, each using two observations to estimate carbon and glass stiffness multipliers; the remaining 18 cases per repeat are validation. Calibration samples are reused across models within a repeat.

## Dataset A cross-validation

| Model | 5-fold CV MAPE | CV RMSE (m/N) | CV R² |
|---|---:|---:|---:|
| Mechanics | 0.1083% | 4.364e-06 | 0.999999442 |
| Quadratic | 0.0053% | 2.797e-07 | 0.999999998 |
| Gaussian Process | 0.0195% | 2.171e-06 | 0.999999862 |

The in-sample MAPE values are also very small (mechanics 0.1054%, quadratic 0.00435%, GP 0.0162%), but the CV figures are more informative. The quadratic model was best on this sampled grid; the mechanics model is interpretable but has the highest CV error. CV does not test transfer to the withheld target materials or to physical measurements.

## Held-out target validation

The table reports mean MAPE over the 18 held-out cases, averaged across the three calibration repeats. “Nominal” uses unit carbon/glass scale factors; “calibrated” uses the two factors estimated from that repeat’s two observations. Errors are percentages. The parenthetical range is the minimum to maximum MAPE across the three repeats.

| Target set | Model | Nominal MAPE, mean ± SD | Calibrated MAPE, mean ± SD (range) |
|---|---|---:|---:|
| material interpolation | Mechanics | 50.150 ± 0.008% | 0.125 ± 0.061% (0.083–0.195%) |
| material interpolation | Quadratic | 50.000 ± 0.001% | 0.006 ± 0.001% (0.005–0.008%) |
| material interpolation | Gaussian Process | 49.974 ± 0.003% | 0.011 ± 0.002% (0.009–0.013%) |
| material extrapolation | Mechanics | 98.805 ± 0.409% | 1.007 ± 1.553% (0.102–2.800%) |
| material extrapolation | Quadratic | 98.619 ± 0.422% | 0.967 ± 1.656% (0.010–2.879%) |
| material extrapolation | Gaussian Process | 98.586 ± 0.422% | 0.947 ± 1.622% (0.008–2.820%) |
| layup transfer | Mechanics | 0.067 ± 0.005% | 0.122 ± 0.078% (0.067–0.211%) |
| layup transfer | Quadratic | 0.204 ± 0.002% | 0.066 ± 0.036% (0.040–0.107%) |
| layup transfer | Gaussian Process | 0.081 ± 0.002% | 0.060 ± 0.010% (0.048–0.066%) |
| geometry extrapolation | Mechanics | 0.254 ± 0.009% | 0.305 ± 0.155% (0.182–0.479%) |
| geometry extrapolation | Quadratic | 0.283 ± 0.032% | 0.359 ± 0.102% (0.300–0.476%) |
| geometry extrapolation | Gaussian Process | 0.159 ± 0.021% | 0.246 ± 0.090% (0.184–0.349%) |

### Material interpolation: calibration transfers well

The nominal model assumes the base (1,1) material scales, so it misses the target’s (1.5,1.5) stiffness and gives about 50% held-out MAPE for all three model families. After two-point calibration, mean MAPE falls to 0.1246% for mechanics, 0.0063% for quadratic, and 0.0112% for GP. The estimated scales are close to (1.5,1.5) and stable across the three pairs. This is the clearest positive result for the proposed limited-observation transfer method.

### Material extrapolation: good when calibration points are informative

For target factors (2.25,0.75), nominal MAPE is about 98.6–98.8%. Calibrated results are generally excellent, but one of the three pairs produces an ambiguous estimate and raises the across-repeat mean. For repeat 2 the selected cases are IDs 37 and 39: both have 32 plies, at widths 0.06 m and 0.10 m. The mechanics calibration Jacobian is rank 1 with condition number about 2.2×10¹⁶; quadratic and GP also give unstable scales (about 1.93,1.94 instead of 2.25,0.75). That pair changes width but not ply count. Since width scales bending and shear compliance similarly, it supplies little information to separate carbon and glass stiffness effects. On the other two pairs, fitted factors are close to the configured target and validation MAPE is small.

This shows why two randomly selected points are not enough as a design rule. Select calibration geometries to make the sensitivities to the two unknown factors distinct—for example, vary ply count as well as width—and check Jacobian rank/conditioning before trusting the fitted factors. The 3-pair experiment exposes this failure mode but is too few repeats to quantify its probability.

### Layup transfer: effective section properties largely carry over

The target stack changes to UD/twill/glass while material scales remain nominal. Before calibration, MAPE is already 0.067% (mechanics), 0.204% (quadratic) and 0.081% (GP). Calibration modestly improves quadratic and GP on average, but worsens the mechanics model from 0.067% to 0.122%. This suggests that the computed effective bending/shear properties represent the alternate layup well for this response, and that a two-point calibration can add noise or absorb small surrogate residuals when no meaningful material correction is needed.

### Geometry extrapolation: small errors, little gain from calibration

With base material factors, nominal MAPE is 0.159–0.283% across models. After calibration it is 0.246–0.359%; the adjustment generally worsens the mean. Region-specific MAPE averaged over models and repeats remains below about 1.1% for the width- and ply-count-extrapolation partitions, but those partitions contain only about 5–7 validation points per repeat and overlap. These results are promising for this particular linear-static FEA design, but should not be described as general extrapolation capability without broader geometries and verification.

R² is near 1 for most target results because the compliance range is broad; it can look excellent while percentage errors or individual errors still matter. Use MAPE/maximum percentage error and direct plots of predicted versus FEA compliance alongside R².

## FEA validity and limits

The sweep is a linear-static shell analysis. It does not include impact dynamics, progressive damage, failure, contact, or the curved landing-gear geometry. The material cards are computational assumptions. The fitted models are being evaluated against FEA generated using the same modeling framework, so agreement demonstrates surrogate fitting and property-transfer behavior relative to that solver setup only.

The reported displacement ratio reaches 1.27 in Dataset A and 1.70 in targets; the most extreme target is the 8-ply, 0.015 m-wide case (tip deflection ≈0.425 m for a 0.25 m beam). 21 target cases and 22 training cases exceed `|delta|/L = 0.1`. The linear analysis still returns a compliance, but a small-deflection physical interpretation is not justified for these points. Because the solve is linear, compliance scales independently of applied force, so reducing the load scales down `|delta|/L` without changing the ideal linear compliance; alternatively, exclude/redefine the thinnest and narrowest combinations.

A 40 × 6 mesh was used for all cases. The archive does not include a mesh-convergence result for these exact property-transfer cases, and the dataset alone cannot certify mesh convergence, boundary-condition fidelity, or absence of locking/singularities. Check representative baseline and extreme cases before relying on sub-percent surrogate errors.

## Recommendations

1. **Fix the small-deflection validity range.** Reduce the applied load or restrict the combinations so `|delta|/L` remains suitably small for the intended linear benchmark; re-run the FEA and refit if the target domain changes.
2. **Use informative calibration pairs.** Choose two points differing in ply count/thickness and width, chosen to separate carbon- and glass-scale sensitivities. Add Jacobian singular values/rank as a precondition for reporting parameter estimates.
3. **Keep calibration and validation separate.** Continue excluding calibration cases from metrics. Add more random/space-filling calibration pairs if evaluating robustness to observation selection.
4. **Do mesh checks on extremes.** Compare compliance for refined meshes at the narrow/thin and wide/thick extremes, and verify total reaction equals 10 N and sign/local axis conventions are correct.
5. **Report model comparison with held-out metrics.** Quadratic is strongest on Dataset A CV and material interpolation; GP is competitive; mechanics-informed regression is easiest to interpret. No single model dominates every target type.
6. **Frame scope carefully.** Present this as a synthetic FEA-to-FEA methodology demonstration. Applying the same workflow to physical tests will require physical observations and a defensible FEA model; these results do not establish physical landing-gear accuracy.

## Files included with this report

- `dataset_A.csv`, `targets.csv`: case-level FEA responses and calculated section properties.
- `model/calibration_splits.csv`: selected calibration case IDs.
- `model/metrics_by_repeat.csv`: all held-out metrics, fitted factors and identifiability diagnostics.
- `model/predictions.csv`: case-level nominal and calibrated predictions.
- `model/summary.json`: dataset counts and Dataset A fit/CV metrics.
- `validation_summary.csv`: compact summary table used in this report.

The fitted `.joblib` objects from the ZIP are not included in GitHub results; the CSV outputs and JSON summary preserve the reported predictions and evaluation evidence.
