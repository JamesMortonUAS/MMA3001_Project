"""Generate static HTML API documentation using Python's standard pydoc."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path(__file__).resolve().parent / "api"
MODULES = (
    "analytic", "fit_corrected_surrogate", "generate_reference",
    "laminate_section", "optimize_mass", "plot_results",
    "run_laminate_shell_sweep", "run_pymapdl_sweep",
)


def main() -> None:
    """Write one pydoc HTML page for each documented project module."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    for module in MODULES:
        subprocess.run([sys.executable, "-m", "pydoc", "-w", module],
                       cwd=OUTPUT, env=env, check=True)
        page = OUTPUT / f"{module}.html"
        html = page.read_text(encoding="utf-8")
        source = (ROOT / "src" / f"{module}.py").resolve()
        html = html.replace(f'file:{source}', f'../../src/{module}.py')
        html = html.replace(str(source), f"../../src/{module}.py")
        page.write_text(html, encoding="utf-8")
    links = "\n".join(f'<li><a href="{name}.html">{name}</a></li>' for name in MODULES)
    (OUTPUT / "index.html").write_text(
        "<!doctype html>\n<html><head><meta charset=\"utf-8\"><title>Project API documentation</title></head>"
        f"<body><h1>Project API documentation</h1><ul>{links}</ul></body></html>\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
