# 10 N quadratic sweep — direct stiffness calibration

Results from the revised direct effective-stiffness calibration, run on 17 October 2026. This uses the same 120-case PyMAPDL SHELL181 sweep as the preceding one-point compliance-multiplier analysis, but reruns the fitting/calibration step with the current direct-stiffness workflow.

## Run and split

- Static linear FEA, 10 N transverse tip load.
- 120 cases; 96 training and 24 held-out cases.
- Five inputs: beam length, root width, taper ratio, thickness scale, and effective stiffness scale.
- Full quadratic response surface with 21 coefficients.
- One analytical calibration observation, case 109; the test cases were not used for fitting or calibration.
- 20 plies; structured SHELL181 mesh with 40 by 6 elements.

## Calibration and held-out results

Least squares estimated a shared effective stiffness scale of 0.97609, relative to the analytical reference value 1.0. The geometry test holds thickness and reference stiffness fixed and compares predictions on 24 unseen geometries.

| Held-out comparison | RMSE (m/N) | MAE (m/N) | R² | MAPE | Maximum absolute percentage error |
|---|---:|---:|---:|---:|---:|
| Quadratic vs. FEA at DOE inputs | 3.624e-4 | 2.712e-4 | 0.9876 | 5.929% | 16.401% |
| Nominal stiffness vs. analytical | 2.668e-4 | 2.246e-4 | 0.9777 | 5.323% | 15.359% |
| Calibrated stiffness vs. analytical | 3.040e-4 | 2.518e-4 | 0.9710 | 5.128% | 12.168% |

Calibration slightly improves percentage-based errors (MAPE and maximum absolute percentage error), but worsens RMSE, MAE, and R². Thus, it improves some held-out error measures, not performance across all measures. One calibration point is the minimum for one scalar and does not establish robustness; repeat with several calibration points and alternative held-out splits before drawing a general conclusion.

## Limits

The maximum absolute deflection-to-length ratio is approximately 0.330, so the most flexible cases may be outside the small-deflection theory assumption. Only one mesh density is recorded. The archive has no MAPDL/PyMAPDL version record or solver log. Treat the numbers as preliminary method results, not physical validation or landing-gear performance predictions.
