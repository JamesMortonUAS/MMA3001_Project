# Preliminary 10 N quadratic sweep — 17 October 2026

This run is archived before the revised direct effective-stiffness calibration was used. The fitting metadata identifies the calibration as a one-point least-squares fit of a compliance multiplier (case 79); it is not the output of the current direct-stiffness calibration workflow. Keep it for comparison when rerunning the sweep.

## Contents and run settings

- `quadratic_sweep.csv`: 120 PyMAPDL SHELL181 cases; static, linear, 10 N tip load.
- Five inputs varied: length, root width, taper ratio, thickness scale, and effective stiffness scale.
- 20 plies; 40 by 6 shell mesh.
- `model/metrics.json`: 96 training cases, 24 held out, one analytical calibration case (case 79).
- `model/test_predictions.csv`, `model/results_table.tex`, and `model/quadratic_model.joblib`: saved fit outputs.

## Main finding

The uncorrected quadratic surface had 5.946% MAPE against the analytical reference on the held-out set. Applying the single fitted compliance multiplier increased MAPE to 7.265%. The one-point correction therefore did not transfer successfully in this run. It estimated a compliance multiplier of 1.06267, equivalent in the model to a shared effective stiffness factor of 0.94102.

The held-out deflection-to-length ratio ranged from approximately 0.034 to 0.330. The largest cases may not satisfy small-deflection beam assumptions. The archive does not include the MAPDL/PyMAPDL version record or solver log. Treat these outputs as preliminary and rerun with the current repository scripts before using them as final results.
