# Project brief traceability

This table maps the MMA3001 2026 project brief to the current repository state.
“In progress” means the repository contains a draft or initial evidence that
still needs student review or additional engineering work before submission.

| Brief requirement | Repository evidence | Status / remaining action |
|---|---|---|
| Define a self-selected engineering problem and its scope | `README.md`; `latex/project_report.tex` | Drafted. Confirm coordinator approval for the chosen project scope and record the decision. |
| Define inputs, outputs, units, valid domains, and invalid cases | `config.example.json`; `layup.json`; report Section 3 | Initial definitions included. Finalize defensible design bounds and input validation rules. |
| Implement and explain a computational method | `src/`; generated HTML in `docs/api/` | Scripts and API documentation included. Verify locally against the intended MAPDL installation. |
| Compare with a credible alternative or baseline | Quadratic response surface and GP results in `results/laminate_run_2026-10-17/model/metrics.json` | Completed for the current 50-case preliminary run. |
| Validate the method and report metrics and limits | `latex/project_report.tex`; `docs/test_report.md`; archived metrics and predictions | Initial analytical benchmark comparison is reported. The 100 N cases violate small-deflection assumptions, so re-run at a justified lower load before making engineering claims. |
| Analyse optimisation / performance and profile computation where relevant | `src/optimize_mass.py`; report Section 7 | Optimisation code is available but no final optimal design is claimed. Runtime profiling and a design requirement are still needed. |
| Provide functional source code, comments, and documentation | `src/`; `docs/api/`; docstrings | Included; student should inspect generated documentation and scripts. |
| Provide test functions and a test report | `tests/test_laminate_section.py`; `docs/test_report.md`; GitHub Actions workflow | 3 local pytest checks pass; CI also runs pytest. Tests do not launch MAPDL. |
| Provide a licensing agreement | `LICENSE` | MIT license added for original project code and docs. Third-party materials and solver outputs are excluded from that grant. |
| Submit written project documentation within the brief's length limit | `latex/project_report.tex` | Draft targets 5–10 pages after compilation. Review pagination, references, and all bracketed student actions. |
| Disclose AI use and reflect on it | Report Section 9 | Disclosure scaffold is included. Student must add a personal reflection based on their actual use and verify all technical claims. |
| Prepare a short presentation and answer questions | Project brief; report and repository evidence | Not yet prepared. Build a 5-minute presentation with validation, decisions, limitations, AI use, and GitHub demonstration. |

## Evidence boundary

The archived sweep is recorded as “PyMAPDL FEA (assumed)” in its metrics because
the supplied archive did not include the solver log, MAPDL release, or a run
manifest. It is useful preliminary evidence for exercising the workflow, but it
does not establish mesh convergence, physical accuracy, or a valid small-
deflection response at 100 N. No physical test data are part of this project.
