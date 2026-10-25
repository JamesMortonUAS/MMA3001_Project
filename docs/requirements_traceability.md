# Project brief traceability

This table maps the MMA3001 project brief to the submission snapshot and distinguishes completed evidence from proposed work.

| Brief requirement | Repository evidence | Status |
|---|---|---|
| Define an engineering problem and project scope | `README.md`; report Introduction | Complete. The report frames a straight composite cantilever as a methodology benchmark and identifies curved gear application as future work. |
| Define inputs, outputs, assumptions, and domain | `docs/extended_transfer_experiment.md`; report Sections 2–4; example configuration | Specified. The nonlinear design uses length, load, width, ply count, material scales, layup, and section properties; output is tip compliance/deflection. Assumed material cards are numerical inputs, not measured properties. |
| Implement a computational method | `src/run_extended_transfer_study.py`; `src/laminate_section.py`; `src/fit_corrected_surrogate.py` | Implemented and documented. The MAPDL model uses layered SHELL181 elements with large-deflection static analysis. |
| Create a planned design of experiments | `docs/extended_transfer_experiment.md`; design-generation code and tests | Specified: 80 training cases and five 20-case target sets (180 planned analyses total). |
| Fit and compare predictive models | `src/fit_corrected_surrogate.py`; report Methods and historical Results | Method implemented for mechanics-informed regression, quadratic response surface, and Gaussian process. Comparative quantitative metrics in the report are from the historical linear-static baseline only. |
| Calibrate with limited observations and assess unseen cases | Report Section 5; historical result tables | Demonstrated on the historical FEA-to-FEA baseline using two observations per repeat and 18 held-out cases. Pair identifiability problems are discussed. |
| Validate and state limits | Report Discussion; `docs/test_report.md` | Unit-level checks and historical baseline validation are documented. The full large-deflection experiment has no complete dataset or validation results in this submission; no physical validation is claimed. |
| Analyse engineering performance or optimise | `src/optimize_mass.py`; report scope and next steps | Optimisation is supplementary code only. No mass-optimal design, landing stroke, impact-energy absorption, or curved-gear performance is claimed. |
| Provide source, documentation, license, and tests | `src/`; `docs/`; `tests/`; `docs/api/`; `LICENSE` | Included. Tests are solver-independent and do not establish MAPDL solution accuracy. |
| Provide a written report of suitable length | `latex/project_report.pdf`; `latex/project_report.tex` | Included. The compiled report is 7 pages total, within the brief's 5–10 page range (excluding references/appendices). |
| Reflect on AI use | Report Section 8 | Included. It describes assistance used, human decisions, verification, and limitations; the student should confirm the reflection matches their own experience. |
| Prepare presentation and submit repository | Repository root and project brief | Repository content is organized for handover. The student should include the public repository URL and a repository ZIP in the Moodle submission as requested by the brief. |

## Evidence boundary

The report includes only the completed historical linear-static, four-set FEA-to-FEA results as quantitative model evidence. The new large-deflection experiment is documented as a reproducible design and implementation, but its complete dataset and model evaluation are unavailable in this submission. No incomplete nonlinear-run outputs are used in the report. No physical test data are included.
