# Test and report build record

Recorded 25 October 2026.

| Check | Result |
|---|---|
| `python -m unittest discover -s tests -v` | 11 passed, 0 failed |
| Python syntax check of the transfer runner and API builder | Passed |
| `python scripts/build_api_docs.py` | Passed |
| Two-pass LaTeX build of `latex/project_report.tex` | Passed; 7-page PDF |

The unit checks exercise source-level design and model logic in the maintained source tree. They do not launch MAPDL, verify solver accuracy, or independently reproduce the archived sweep.

The submission report describes the completed 80-training/80-target linear-static experiment and its model evaluation, archived under `results/property_transfer_extended/`. Source-level checks are distinct from validation against the held-out target responses. The report documents limits on mesh evidence and physical small-deflection interpretation.
