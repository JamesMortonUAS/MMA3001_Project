# Project brief traceability

This table maps the MMA3001 project brief to the current hidden-property FEA transfer study. “In progress” means the repository contains a method or draft, but further analysis or student review is still required.

| Brief requirement | Repository evidence | Status / remaining action |
|---|---|---|
| Define the engineering problem and scope | `README.md`; `latex/project_report.tex` | Benchmark and limits are defined. Confirm coordinator approval of the chosen scope and record the decision. |
| Define inputs, outputs, units, and valid domain | `config.hidden_property.example.json`; `layup.json`; report Sections 2–4 | SI parameters and 97-case study design are specified. Check the engineering relevance of bounds and assess small-deflection assumptions after running FEA. |
| Implement and explain a computational method | `src/run_property_transfer_study.py`; `src/run_laminate_shell_sweep.py`; `src/laminate_section.py`; `docs/property_transfer_study.md` | Scripts implement the regular-grid FEA generation, section-property regression, one-point common-scale calibration, and regional validation workflow. Verify locally with the intended MAPDL installation. |
| Compare against a baseline or reference | `src/run_property_transfer_study.py` | Dataset A five-fold cross-validation and nominal-vs-calibrated predictions are specified. Generate Dataset A/B and report actual metrics; no current grid results are committed. |
| Validate and state limitations | `latex/project_report.tex`; `docs/test_report.md` | Validation protocol separates interpolation and extrapolation; this is FEA-to-FEA only. Run the sweep, examine mesh sensitivity and `abs(delta)/L`, and do not claim physical validation. |
| Analyse performance or optimisation | `src/optimize_mass.py`; report scope | Optimisation code is supplementary. No mass-optimal design is claimed by the hidden-property transfer experiment; define requirements and verify the predictor first. |
| Provide functional source and documentation | `src/`; `docs/`; `README.md` | Included. Student should inspect scripts, API documentation, and local results. |
| Provide tests and test evidence | `tests/`; `docs/test_report.md`; `.github/workflows/` | Solver-independent tests do not launch MAPDL. Run them and retain the current test report. |
| Provide licensing | `LICENSE` | MIT license covers original project material as stated in the file; check third-party material and solver-output provenance separately. |
| Provide written project documentation | `latex/project_report.tex`; compiled `latex/project_report.pdf` | Updated method draft is six pages before adding current FEA results. Recheck length and pagination after results are added. |
| Disclose AI use and reflect | Report Section 8 | Use disclosure is included; student should add a personal reflection on assistance, verification, and decisions. |
| Prepare and present the project | Project brief; report and repository | Prepare a short presentation with problem, model, data design, calibration, evaluation, limitations, and evidence from the completed run. |

## Evidence boundary

The configured current workflow contains 72 Dataset A training analyses and 25 Dataset B analyses. The current hidden-property grid has not yet been run in the repository, so no fitted correction or validation metrics are available. Earlier studies in `results/` used different setups and should remain identified as historical or complementary. No physical test data are included.
