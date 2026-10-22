# Hidden-property grid sweep: results and interpretation

## Summary

The archived run completed all **97 planned MAPDL cases**: 72 Dataset A training cases and 25 Dataset B cases at the withheld effective-property pair. The two-term mechanics-informed model fits Dataset A closely. One Dataset B observation reduced the mean held-out absolute percentage error from **9.945% to 0.0539%** across the other 24 cases. It performed similarly on the eight geometry-interpolation points and sixteen geometry-extrapolation points.

This is a strong demonstration of response transfer within this controlled setup. It is **FEA-to-FEA only**: both datasets used the same shell model, solver approach, materials framework, and mesh. It does not validate material properties against physical gears or establish accuracy for a curved landing gear.

## Study configuration

- Dataset A: full 3-level grid in width (0.03, 0.05, 0.07 m), thickness scale (0.7, 1.0, 1.3), in-plane stiffness scale (0.7, 0.9, 1.1), and transverse-shear scale (0.7, 1.15, 1.6). The target property pair (0.9, 1.15) was withheld, leaving 72 of 81 combinations.
- Dataset B: 5×5 width–thickness grid at the withheld pair (0.9, 1.15). One central point (case 12; width 0.05 m, thickness scale 1.0) calibrated a **single shared stiffness multiplier**. It did not independently identify in-plane and shear scales.
- Training response: compliance, $C=|u_{z,tip}|/10\,\mathrm{N}$, from a static, linear 10 N SHELL181 model, length 0.25 m, mesh 40×6, layered 20-ply shell.
- Model: non-negative, no-intercept least squares in bending and transverse-shear compliance basis terms. The fitted Dataset A coefficients were 0.99854 (bending) and 1.28480 (shear). The calibrated shared stiffness scale was 0.90054.

The nominal (uncalibrated) model underpredicted Dataset B compliance by about 9.94% throughout the held-out set. The one-point fit removed this common bias. The fitted scale is near the withheld in-plane scale of 0.9, but the true shear scale was 1.15. That does **not** mean the calibration recovered either material property: one scalar observation cannot identify two independent property scales. A common correction is response-equivalent only. The good transfer suggests that the error over this range was largely a shared response bias and that the fitted mechanics basis captured the width/thickness trend.

## Results

Dataset A five-fold cross-validation and held-out Dataset B metrics are summarized below. Validation excludes the calibration observation.

| Data / region | n | RMSE (m/N) | MAE (m/N) | MAPE | Max absolute percentage error | R² |
|---|---:|---:|---:|---:|---:|---:|
| Dataset A, 5-fold CV | 72 | 6.82×10⁻⁷ | 5.22×10⁻⁷ | 0.0571% | 0.1156% | 0.9999995 |
| Dataset B, nominal, all validation | 24 | 1.864×10⁻⁴ | 1.272×10⁻⁴ | 9.9447% | 10.0361% | 0.981214 |
| Dataset B, calibrated, all validation | 24 | 1.442×10⁻⁶ | 7.584×10⁻⁷ | 0.05393% | 0.10646% | 0.9999989 |
| Interpolation, nominal | 8 | 9.069×10⁻⁵ | 7.915×10⁻⁵ | 9.9454% | 9.9916% | 0.958086 |
| Interpolation, calibrated | 8 | 3.591×10⁻⁷ | 2.709×10⁻⁷ | 0.03369% | 0.05253% | 0.9999993 |
| Extrapolation, nominal | 16 | 2.191×10⁻⁴ | 1.513×10⁻⁴ | 9.9443% | 10.0361% | 0.980821 |
| Extrapolation, calibrated | 16 | 1.748×10⁻⁶ | 1.002×10⁻⁶ | 0.06405% | 0.10646% | 0.9999988 |

The calibrated RMSE fell by approximately 99.23% overall. Absolute errors are reported in m/N; percentage metrics use pointwise relative errors. Case 12 was matched by calibration and is excluded from all validation rows. Its nominal error was −9.9457%; its calibrated residual was effectively zero by construction.

The “extrapolation” result is limited to the tested rectangular domain: width 0.025–0.075 m and thickness scale 0.6–1.4, versus Dataset A bounds 0.03–0.07 m and 0.7–1.3. Sixteen points lie outside at least one bound. It does not support arbitrary extrapolation to other geometry, length, layup, boundary conditions, load regime, or curvature.

![Held-out parity and residuals](grid_sweep_validation.svg)

## Interpretation and limitations

1. **Calibration worked for this controlled test.** There was a nearly uniform approximately 10% nominal underprediction. A single common correction removed it, with sub-0.11% maximum relative error on the held-out grid.
2. **The result tests the methodology under favorable shared-model conditions.** Both target and training values come from the same finite-element formulation and parameterization. Real physical data would add material, manufacturing, boundary-condition, and model-form discrepancies not tested here.
3. **Do not claim recovered material properties.** The fitted 0.90054 multiplier is not a separate estimate of the hidden in-plane and shear scales. More independent observations and/or a parameterized correction with identifiability analysis would be required to infer both.
4. **Large normalized deflections limit physical interpretation.** Dataset B $|u|/L$ spans about 0.006–0.234; Dataset A spans about 0.0067–0.158. The softest cases are not within a convincing small-deflection range. MAPDL was run as linear static; these results are numerical transfer evidence, not a large-deflection gear prediction.
5. **Mesh sensitivity and solver provenance remain necessary.** This sweep used one 40×6 mesh. Mesh convergence, MAPDL/PyMAPDL versions, solver logs, element settings and load/support details should accompany any final submission.

## Conclusion and next steps

The experiment demonstrates the intended sequence—fit from a systematic FEA grid, calibrate one response-equivalent stiffness factor from one target observation, then score withheld interpolation and bounded extrapolation cases. Within this synthetic FEA-to-FEA setting, calibration greatly improves predictions. The next useful checks are mesh refinement at soft and stiff extremes, repeat with a second hidden property pair and different calibration locations, examine a response/load level that keeps $|u|/L$ small, and preserve an independent physical-data application as future work. No physical validation is claimed.
