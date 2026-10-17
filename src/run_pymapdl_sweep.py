"""Run a Latin-hypercube PyMAPDL/BEAM188 sweep for the beam benchmark."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc


def make_design(config: dict) -> pd.DataFrame:
    base, sweep = config["baseline"], config["sweep"]
    n = int(sweep["n_design_points"])
    if n < 12:
        raise ValueError("Use at least 12 design points for a 3-input initial fit")
    unit = qmc.LatinHypercube(d=3, seed=int(sweep["seed"])).random(n)
    lscale = sweep["length_scale"][0] + unit[:, 0] * (sweep["length_scale"][1] - sweep["length_scale"][0])
    bscale = sweep["root_width_scale"][0] + unit[:, 1] * (sweep["root_width_scale"][1] - sweep["root_width_scale"][0])
    rb = sweep["taper_ratio"][0] + unit[:, 2] * (sweep["taper_ratio"][1] - sweep["taper_ratio"][0])
    return pd.DataFrame({
        "case_id": np.arange(n),
        "length_m": base["length_m"] * lscale,
        "root_width_m": base["root_width_m"] * bscale,
        "thickness_m": base["thickness_m"],
        "taper_ratio": rb,
        "tip_load_N": base["tip_load_N"],
        "youngs_modulus_Pa": base["youngs_modulus_Pa"],
        "poissons_ratio": base["poissons_ratio"],
    })


def solve_case(mapdl, row, n_elements: int) -> float:
    """Solve one BEAM188 model; return signed tip UZ in metres."""
    mapdl.clear()
    mapdl.prep7()
    mapdl.et(1, "BEAM188")
    mapdl.keyopt(1, 3, 3)  # cubic interpolation along beam axis
    mapdl.mp("EX", 1, float(row.youngs_modulus_Pa))
    mapdl.mp("PRXY", 1, float(row.poissons_ratio))

    length = float(row.length_m)
    b0 = float(row.root_width_m)
    t = float(row.thickness_m)
    rb = float(row.taper_ratio)
    p = float(row.tip_load_N)

    # Nodes along global X; section width is local Y (global Y), thickness
    # is local Z (global Z), so FZ produces bending through thickness t.
    for i in range(n_elements + 1):
        mapdl.n(i + 1, length * i / n_elements, 0.0, 0.0)
    orient = n_elements + 2
    mapdl.n(orient, 0.0, 1.0, 0.0)

    # Piecewise-constant cross-section per element, width evaluated at the
    # element midpoint. Mesh-convergence checks are required before use.
    for i in range(n_elements):
        xmid = length * (i + 0.5) / n_elements
        width = b0 * (1.0 - (1.0 - rb) * xmid / length)
        sec = i + 1
        mapdl.sectype(sec, "BEAM", "RECT")
        mapdl.secdata(width, t, 2, 2)
        mapdl.type(1)
        mapdl.mat(1)
        mapdl.secnum(sec)
        mapdl.e(i + 1, i + 2, orient)

    mapdl.finish()
    mapdl.slashsolu()
    mapdl.antype("STATIC")
    mapdl.d(1, "ALL", 0.0)
    mapdl.f(n_elements + 1, "FZ", -p)
    mapdl.solve()
    mapdl.finish()
    mapdl.post1()
    mapdl.set(1, 1)
    uz = float(mapdl.get_value("NODE", n_elements + 1, "U", "Z"))
    mapdl.finish()
    return uz


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--output", default="results/fea_sweep.csv")
    ap.add_argument("--elements", type=int, default=None)
    args = ap.parse_args()
    config = json.loads(Path(args.config).read_text())
    n_elements = args.elements or int(config["fea"]["elements"])
    if n_elements < 4:
        raise ValueError("Use at least four beam elements")

    from ansys.mapdl.core import launch_mapdl

    design = make_design(config)
    mapdl = launch_mapdl(nproc=int(config["fea"].get("nproc", 2)))
    results = []
    try:
        for row in design.itertuples(index=False):
            uz = solve_case(mapdl, row, n_elements)
            results.append({**row._asdict(), "n_elements": n_elements,
                            "tip_deflection_m": uz,
                            "compliance_m_per_N": abs(uz) / float(row.tip_load_N)})
            print(f"case {row.case_id:03d}: |UZ|={abs(uz):.6e} m")
    finally:
        mapdl.exit()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(out, index=False)
    print(f"Saved {len(results)} cases to {out}")


if __name__ == "__main__":
    main()
