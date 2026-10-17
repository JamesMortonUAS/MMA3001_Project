# Calibration observation sensitivity — 17 October 2026

Post-processing of the 120-case 10 N SHELL181 sweep in `results/quadratic_run_2026-10-17_direct_stiffness_calibration/`. No new MAPDL analyses were run.

## Method

The script fits the five-input quadratic response surface once using the configured 96-case training split. It then estimates one shared effective stiffness scale for 1, 2, 3, and 5 analytical observations. Each count has 30 unique training-only subsets: one centred/maximin selection and 29 random selections. Every subset is scored against the same 24 held-out geometries. Held-out analytical references are used only for scoring.

## Findings

The nominal held-out MAPE is 5.323%. Across 30 subsets at each calibration count, mean and median MAPE increased for all four counts. Only 30.0%, 33.3%, 30.0%, and 36.7% of subsets improved MAPE for 1, 2, 3, and 5 observations respectively. Mean RMSE rose by 53.0%, 22.8%, 20.5%, and 8.4%, respectively. The estimated scale became more stable as calibration observations increased (scale standard deviation 0.0462 for one observation and 0.0164 for five), but that did not translate to better prediction accuracy on average.

This analysis is conditional on one FEA-trained surface and one train/test split. It demonstrates sensitivity to calibration-case choice, not universal behavior or physical validation.

## Files

- `calibration_sensitivity_runs.csv`: one row per calibration subset.
- `calibration_sensitivity_summary.csv`: grouped summaries by calibration count.
- `calibration_sensitivity_metrics.json`: fixed split, case IDs, settings, and baseline metrics.
- `calibration_sensitivity.png`: error-change and fitted-scale distributions.
