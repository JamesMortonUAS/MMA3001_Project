"""Generate lightweight, dependency-free HTML API docs from Python docstrings.

Run from the repository root with ``python scripts/build_api_docs.py``. The
builder parses source with the standard-library AST, so documentation can be
regenerated without importing PyMAPDL or launching Ansys.
"""
from __future__ import annotations

import ast
from html import escape
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "run_extended_transfer_study.py"
OUT_DIR = ROOT / "docs" / "api"
PAGE = OUT_DIR / "run_extended_transfer_study.html"
INDEX = OUT_DIR / "index.html"


def render_api_page(source: str) -> str:
    """Render module and top-level function docstrings as a standalone page."""
    tree = ast.parse(source)
    module_doc = ast.get_docstring(tree) or ""
    entries = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            signature = f"{node.name}({ast.unparse(node.args)})"
            if node.returns:
                signature += f" -> {ast.unparse(node.returns)}"
            doc = ast.get_docstring(node) or "No function docstring is currently provided."
            entries.append(
                "<section><h2><code>" + escape(signature) + "</code></h2><pre>"
                + escape(doc) + "</pre></section>"
            )
    body = "\n".join(entries)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>run_extended_transfer_study API</title>
<style>body{{font:16px/1.55 system-ui,sans-serif;max-width:960px;margin:2rem auto;padding:0 1rem;color:#18212b}}
h1,h2{{line-height:1.25}}section{{border-top:1px solid #d6dde5;padding:.6rem 0}}code,pre{{font-family:ui-monospace,monospace}}
pre{{white-space:pre-wrap;background:#f4f6f8;padding:.8rem;border-radius:5px}}a{{color:#075a9c}}</style></head>
<body><p><a href="index.html">API documentation index</a></p>
<h1>run_extended_transfer_study</h1><pre>{escape(module_doc)}</pre>
<p>Generated from <code>src/run_extended_transfer_study.py</code> by
<code>scripts/build_api_docs.py</code>. Regenerate after source or docstring changes.</p>
{body}</body></html>"""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PAGE.write_text(render_api_page(SOURCE.read_text(encoding="utf-8")), encoding="utf-8")
    names = {p.stem for p in OUT_DIR.glob("*.html") if p.name != "index.html"}
    if INDEX.exists():
        names.update(Path(href).stem for href in re.findall(r'href="([^"]+\.html)"', INDEX.read_text(encoding="utf-8")))
    else:
        names.update({"analytic", "fit_corrected_surrogate", "laminate_section", "optimize_mass", "plot_results", "run_laminate_shell_sweep"})
    names.add(PAGE.stem)
    names = sorted(names)
    links = "\n".join(f'<li><a href="{escape(name)}.html">{escape(name)}</a></li>' for name in names)
    INDEX.write_text(
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><title>Project API documentation</title></head>"
        f"<body><h1>Project API documentation</h1><ul>{links}</ul></body></html>\n",
        encoding="utf-8",
    )
    print(f"Generated {PAGE.relative_to(ROOT)} and updated docs/api/index.html ({len(names)} pages).")


if __name__ == "__main__":
    main()
