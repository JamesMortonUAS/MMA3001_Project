# Test and documentation build record

**Recorded:** 24 October 2026  
**Scope:** source-level and model-workflow unit checks for the large-deflection transfer runner.

## Results

| Check | Command | Result |
|---|---|---|
| Python syntax | `python -m py_compile src/run_extended_transfer_study.py scripts/build_api_docs.py` | Passed |
| Unit suite | `python -m unittest discover -s tests -v` | 10 passed, 0 failed |
| API HTML generation | `python scripts/build_api_docs.py` | Passed; runner page and index regenerated |
| LaTeX report build | `pdflatex -interaction=nonstopmode -halt-on-error -output-directory output/pdf project_report_current.tex` | Passed; 7-page PDF generated |

The unit tests cover the 80-case training design and five 20-case target sets, symmetric stacks, effective section values, the resin-ratio stiffness proxy and baseline feature assumption, configuration rejection, invalid predictor handling, large-deflection MAPDL command generation, FEA case-signature changes, and model fit/cross-validation timing instrumentation.

## Limits of these checks

No MAPDL solve was run in this environment. The solver-control unit test checks that the runner issues `NLGEOM,ON`, `AUTOTS,ON`, `NSUBST`, `NEQIT`, and `NROPT,FULL`; it does not establish that the full nonlinear FEA converges or that its displacement is accurate. The 180-case nonlinear sweep, mesh verification, reaction-balance check, and timing comparison must be run on the licensed MAPDL installation. The historical linear-static CSVs are not accepted as nonlinear results.
