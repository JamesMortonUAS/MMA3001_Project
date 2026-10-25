# Project brief traceability

The submitted project is the completed 160-case linear-static composite-beam transfer experiment.

| Requirement | Evidence and scope |
|---|---|
| Engineering problem | Report Introduction: predict composite compliance using an FEA-trained model calibrated with limited observations. |
| Inputs, outputs and data | Report finite-element and dataset sections; `docs/completed_transfer_study.md`; archived Dataset A and target CSVs. |
| Computational approach | Three model families; section-property predictors; two-factor least-squares calibration; documented source and HTML API pages. |
| Validation and comparison | Five-fold training cross-validation; four target sets; three two-observation repeats and 18 held-out cases per repeat; nominal/calibrated metrics. |
| Interpretation and limitations | Report discusses identifiability, known nominal material cards, displacement ratios, mesh evidence, and FEA-to-FEA scope. |
| Optimisation/application | Mass/stroke optimisation and physical curved-gear application remain future work; no optimal gear is claimed. |
| Report | `latex/project_report.pdf` and matching source, 7 pages including references. |
| Code, documentation and license | `src/`, `docs/api/`, README and MIT license. Later exploratory source is distinguished from the submitted completed experiment. |
| Tests | `tests/` and `docs/test_report.md`; 11 unittest checks pass, without launching MAPDL. |
| AI reflection | Report describes assistance, student responsibility, verification and limitations. Student should confirm the reflection before submission. |
| Submission | Provide the repository URL and a repository ZIP in Moodle as required by the brief. |

The results establish model agreement with the completed linear FEA dataset. They do not establish physical material properties or landing-gear performance.
