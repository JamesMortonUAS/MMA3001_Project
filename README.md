# MMA3001 Landing Gear Optimisation

This repository demonstrates a limited-observation calibration workflow for FEA-based composite-beam compliance models. The benchmark is a straight cantilever; the intended future application is a curved composite landing-gear leg.

## Submission documents

- [Final project report (PDF)](latex/project_report.pdf)
- [Report source (LaTeX)](latex/project_report.tex)
- [Current experiment definition and run guide](docs/extended_transfer_experiment.md)
- [Project brief traceability](docs/requirements_traceability.md)
- [Test and documentation build record](docs/test_report.md)
- [API documentation](docs/api/index.html)
- [MIT license](LICENSE)

## Scope and evidence status

The completed quantitative evidence in the report is from an earlier linear-static FEA-to-FEA study with four target sets. The report labels these results as historical. They are not results for the updated large-deflection experiment and do not validate physical material properties, landing impact, a curved gear, or mass optimisation.

The updated PyAnsys MAPDL workflow specifies 80 training cases and 100 target cases (five 20-case sets), using a 0.25 m cantilever, a 10 N tip-edge load, layered SHELL181 elements, a 40 by 6 mesh, and large-deflection static analysis. A complete dataset and validation results for that updated experiment are not included in the submission report. No output from an incomplete nonlinear execution is presented as evidence.

The study compares a mechanics-informed compliance regression, a quadratic response surface, and a Gaussian process. Two target observations calibrate carbon and glass stiffness multipliers; the remaining 18 target cases are held out in each repeat. See the run guide for the full design, assumptions, and limits.

## Reproduce the source-level checks

From the repository root:

```powershell
python -m unittest discover -s tests -v
python scripts/build_api_docs.py
```

Running the MAPDL sweep requires a licensed local MAPDL installation and the dependencies in `requirements.txt`. The unit tests do not launch MAPDL, so they do not verify solver convergence or FEA accuracy. Follow the run guide and verify nonlinear convergence, reaction balance, and mesh sensitivity before interpreting any newly generated results.
