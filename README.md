# MMA3001 Landing Gear Optimisation

This project evaluates limited-observation calibration of composite cantilever compliance models using a completed linear-static FEA sweep.

## Submission

- [Project report PDF](latex/project_report.pdf) and [LaTeX source](latex/project_report.tex)
- [Completed experiment definition](docs/completed_transfer_study.md)
- [Detailed results and interpretation](docs/extended_transfer_results_report_2026-10-22.md)
- [Archived sweep and evaluation outputs](results/property_transfer_extended/)
- [Brief traceability](docs/requirements_traceability.md)
- [Test/build record](docs/test_report.md)
- [HTML API documentation](docs/api/index.html)
- [License](LICENSE)

## Completed study

The experiment contains 80 training analyses and 80 target analyses: material interpolation, material extrapolation, layup transfer, and geometry extrapolation. Each target set has 20 cases. A mechanics-informed regression, quadratic response surface, and Gaussian process are trained on Dataset A. Each model is calibrated using two target observations and evaluated on the remaining 18 cases, with three calibration repeats per target set.

The FEA uses a 0.25 m cantilever, a fully fixed root, a 10 N load distributed across the tip edge, and layered SHELL181 elements on a 40 by 6 mesh. Material cards are numerical assumptions. Results demonstrate agreement with linear FEA, not physical landing-gear accuracy. Extreme displacement ratios and calibration identifiability are discussed in the report.

The report and archived evaluation outputs define the submitted study. The source tree also retains later exploratory development; those changes must not be used to regenerate the completed sweep under a different analysis definition.

## Checks

```powershell
python -m unittest discover -s tests -v
python scripts/build_api_docs.py
```

The recorded unit check passes 11 tests. Tests do not launch MAPDL or verify physical accuracy. Python dependencies are in `requirements.txt`; additional test dependencies are in `requirements-dev.txt`.

## Future work

Potential extensions include informative calibration-pair selection, noise robustness, mesh verification, large-deflection analysis, physical testing, curved gear geometry, and mass/stroke optimisation.
