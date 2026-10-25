# Archived partial large-deflection sweep

These CSV files are the run outputs supplied on 25 October 2026. They are archived unchanged. This is an incomplete run, not the final model-validation dataset.

- `dataset_A.csv`: 80 successful training FEA responses.
- `dataset_A_design.csv`: the 80 planned Dataset A cases.
- `targets.csv`: 81 successful target responses (case IDs 0–80).
- `target_design.csv`: all 100 planned target cases; IDs 81–99 have no result.

The run stopped during target case 81 after MAPDL reported excessive distortion for elements 7 and 121. No model fit or calibration results were produced. See [the run report](../../../docs/large_deflection_transfer_partial_run_report_2026-10-25.md) for analysis, response summaries, and recommended follow-up.