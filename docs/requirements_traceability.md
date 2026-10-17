# Project brief traceability

This table maps the MMA3001 2026 project brief to the current repository state.
“In progress” means the repository contains a draft or initial evidence that
still needs student review or additional engineering work before submission.

| Brief requirement | Repository evidence | Status / remaining action |
|---|---|---|
| Define a self-selected engineering problem and its scope | `README.md`; `latex/project_report.tex` | Drafted. Confirm coordinator approval for the chosen project scope and record the decision. |
| Define inputs, outputs, units, valid domains, and invalid cases | `config.quadratic.example.json`; `layup.json`; report Section 3 | Five-dimensional SI design domain specified. Review engineering bounds and confirm small-deflection validity after FEA. |
| Implement and explain a computational method | `src/`; generated HTML in `docs/api/` | Scripts and API documentation included. Verify locally against the intended MAPDL installation. |
| Compare with a credible alternative or baseline | `src/fit_corrected_surrogate.py`; generated held-out metrics | The quadratic response surface is the selected model. New 10 N FEA metrics are pending. |
| Validate the method and report metrics and limits | `latex/project_report.tex`; `docs/test_report.md`; generated metrics and predictions | Direct least-squares estimation of one effective stiffness scale and fixed-section geometry-transfer evaluation are implemented. Run the new sweep and assess small-deflection validity. |
| Analyse optimisation / performance and profile computation where relevant | `src/optimize_mass.py`; report Section 7 | Optimisation code is available but no final optimal design is claimed. Runtime profiling and a design requirement are still needed. |
| Provide functional source code, comments, and documentation | `src/`; `docs/api/`; docstrings | Included; student should inspect generated documentation and scripts. |
| Provide test functions and a test report | `tests/test_laminate_section.py`; `docs/test_report.md`; GitHub Actions workflow | 3 local pytest checks pass; CI also runs pytest. Tests do not launch MAPDL. |
| Provide a licensing agreement | `LICENSE` | MIT license added for original project code and docs. Third-party materials and solver outputs are excluded from that grant. |
| Submit written project documentation within the brief's length limit | `latex/project_report.tex` | Draft targets 5–10 pages after compilation. Review pagination, references, and all bracketed student actions. |
| Disclose AI use and reflect on it | Report Section 9 | Disclosure scaffold is included. Student must add a personal reflection based on their actual use and verify all technical claims. |
| Prepare a short presentation and answer questions | Project brief; report and repository evidence | Not yet prepared. Build a 5-minute presentation with validation, decisions, limitations, AI use, and GitHub demonstration. |

## Evidence boundary

The archived 100 N sweep is historical exploratory evidence and is not used as
validation for the current 10 N configuration. The new PyAnsys sweep has not
yet been run, so no new FEA metrics or solver provenance are available. No
physical test data are part of this project.
