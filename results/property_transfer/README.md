# Hidden-property transfer study results

**Run date:** 21 October 2026  
**Branch:** `feat/hidden-property-calibration`

## Summary

The calibrated model predicts the eight held-out Dataset B compliances very accurately. Held-out MAPE falls from **9.97%** at nominal properties to **0.0408%** after calibration. However, this result does **not** show that both hidden stiffness properties were identified: the fitted transverse-shear scale is 0.741, while the hidden target is 1.15. The chosen beam cases have very little shear contribution to total compliance, so the shear scale is weakly identifiable from tip deflection.

The result supports transfer of the predicted tip compliance within this FEA-generated case set. It does not establish recovery of all material properties, physical accuracy, or transfer to curved landing gear.

## Study setup

- Straight cantilever, length **0.25 m**, with a **10 N** transverse tip load.
- Static linear SHELL181 FEA, **40 × 6** mesh, three through-thickness integration points per ply, and the project's 20-ply layup.
- Dataset A: 120 Latin-hypercube cases. Width: 0.03–0.07 m; thickness, in-plane stiffness, and transverse-shear stiffness scales: 0.7–1.3.
- Dataset B: 12 cases from a grid of three widths (0.035, 0.050, 0.065 m) and four thickness scales (0.7, 0.9, 1.1, 1.3). Its hidden target scales were 0.90 in-plane and 1.15 transverse shear.
- Four Dataset B cases (IDs 4, 11, 9, 2) were selected for calibration. The other eight (IDs 0, 1, 3, 5, 6, 7, 8, 10) were held out.
- The fitted model is a positive two-term mechanics-informed regression using bending and shear compliance bases. Its coefficients are learned from Dataset A FEA.

## Results

| Measure | Nominal properties | Calibrated properties |
|---|---:|---:|
| Dataset B holdout RMSE (m/N) | 1.2841 × 10⁻⁴ | 8.6445 × 10⁻⁷ |
| Dataset B holdout MAE (m/N) | 1.0218 × 10⁻⁴ | 5.5216 × 10⁻⁷ |
| Dataset B holdout MAPE | 9.973% | 0.0408% |
| Maximum absolute percentage error | 10.016% | 0.0814% |
| Dataset B holdout R² | 0.97265 | 0.999999 |

Dataset A five-fold cross-validation gave RMSE **3.45 × 10⁻⁷ m/N**, MAPE **0.0357%**, and R² **0.9999995**. This measures interpolation within the same FEA data-generating setup.

| Property scale | Hidden target | Fitted from 4 B cases | Relative recovery error |
|---|---:|---:|---:|
| In-plane stiffness | 0.900 | 0.9007 | +0.07% |
| Transverse-shear stiffness | 1.150 | 0.7409 | −35.6% |

The response model attributes only about **0.036–0.125%** of predicted compliance to its shear term across Dataset B. Changing the shear scale therefore has little effect on the predicted tip response. The excellent holdout scores are mainly evidence that the in-plane scale and total compliance transferred well; they are not evidence that the shear scale was recovered. The calibration objective is poorly conditioned for separating these two scales under this design.

## Interpretation and limitations

1. **Calibration improved the requested response prediction.** It corrected the roughly 10% nominal bias on the eight untouched Dataset B cases.
2. **Property recovery was incomplete.** The in-plane scale was recovered closely, but the shear scale was not. Report the predicted response and fitted parameters separately.
3. **The validation is FEA-to-FEA only.** Both datasets use the same idealized solver setup and material formulation. No physical testing or physical material validation is represented.
4. **The parameter domain is narrow.** Dataset B tests only three widths and four thicknesses at one beam length and one load. The results do not show extrapolation to other lengths, curved gears, or different boundary conditions.
5. **Small-deflection validity needs checking.** Maximum deflection/span is about 12.6% in A and 10.5% in B. These linear static results should not be described as verified small-deflection behavior.
6. **Mesh and solver metadata are incomplete.** These runs use 40 × 6 elements. The earlier mesh-convergence file was for a 0.5 m beam, so it does not establish mesh convergence for this 0.25 m sweep. The archive does not record MAPDL/PyMAPDL versions or solver logs.

## Recommended next study

- Design cases to make the shear term measurable. For a Timoshenko beam, shortening the span or increasing thickness raises the relative shear contribution. If length becomes a variable, include that range in Dataset A as well as Dataset B.
- Select calibration cases using sensitivity to both property scales, not geometry-only maximin spacing. Check the sensitivity/Jacobian matrix before running the full sweep.
- Run a mesh-convergence check for the 0.25 m, 40 × 6 configuration and confirm the linear-response assumption at the largest deflections.
- Keep some Dataset B cases fully hidden until calibration choices and model settings are frozen. Continue to report both holdout response errors and parameter-recovery errors.

## Files

- `dataset_A.csv`: FEA design and responses for Dataset A.
- `dataset_B.csv`: geometry and FEA responses for Dataset B.
- `dataset_B_hidden_truth.json`: target scales for post-calibration audit only.
- `model/metrics.json`: calibration split and summary metrics.
- `model/dataset_b_predictions.csv`: nominal and calibrated predictions for all B cases.
- `model/model_a_calibrated.joblib`: serialized fitted model and calibration scales.
