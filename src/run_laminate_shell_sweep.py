"""PyMAPDL sweep using the explicit 20-ply SHELL181 laminate section."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

from laminate_section import (define_twenty_ply_shell_section, layup_plies,
                              scaled_layup, effective_bending_shear_stiffness)


def make_design(config: dict) -> pd.DataFrame:
    """Create a five-factor Latin-hypercube design for the laminate model."""
    base, sweep = config["baseline"], config["sweep"]
    n = int(sweep["n_design_points"])
    if n < 42:
        raise ValueError("Use at least 42 cases for a five-input quadratic surface")
    unit = qmc.LatinHypercube(d=5, seed=int(sweep["seed"])).random(n)
    ranges = [sweep[k] for k in ("length_scale", "root_width_scale", "taper_ratio",
                                  "thickness_scale", "effective_stiffness_scale")]
    values = [lo + unit[:, i] * (hi - lo) for i, (lo, hi) in enumerate(ranges)]
    return pd.DataFrame({
        "case_id": np.arange(n),
        "length_m": base["length_m"] * values[0],
        "root_width_m": base["root_width_m"] * values[1],
        "taper_ratio": values[2],
        "thickness_scale": values[3],
        "effective_stiffness_scale": values[4],
        "tip_load_N": float(base.get("tip_load_N", 10.0)),
    })


def _node_id(i: int, j: int, n_width: int) -> int:
    return i * (n_width + 1) + j + 1


def solve_case(mapdl, row, layup: dict, n_length: int, n_width: int,
               integration_points: int) -> tuple[float, float, float]:
    """Solve one structured-planform shell model; return tip UZ, mass, h_total."""
    mapdl.clear()
    mapdl.prep7()
    case_layup = scaled_layup(
        layup, float(row.thickness_scale), float(row.effective_stiffness_scale))
    total_thickness = define_twenty_ply_shell_section(
        mapdl, case_layup["material_cards"], section_id=1,
        integration_points_per_ply=integration_points)

    length = float(row.length_m)
    b0 = float(row.root_width_m)
    rb = float(row.taper_ratio)
    load = float(row.tip_load_N)

    # Structured shell grid on a straight, linearly width-tapered planform.
    # Local shell axis 1 follows global +X (beam axis); axis 2 follows +Y.
    for i in range(n_length + 1):
        x = length * i / n_length
        width = b0 * (1.0 - (1.0 - rb) * x / length)
        for j in range(n_width + 1):
            y = -0.5 * width + width * j / n_width
            mapdl.n(_node_id(i, j, n_width), x, y, 0.0)

    mapdl.type(1)
    mapdl.secnum(1)
    for i in range(n_length):
        for j in range(n_width):
            # Counter-clockwise connectivity gives shell normal +Z.
            mapdl.e(_node_id(i, j, n_width),
                    _node_id(i + 1, j, n_width),
                    _node_id(i + 1, j + 1, n_width),
                    _node_id(i, j + 1, n_width))

    root_nodes = [_node_id(0, j, n_width) for j in range(n_width + 1)]
    tip_nodes = [_node_id(n_length, j, n_width) for j in range(n_width + 1)]
    mapdl.finish()
    mapdl.slashsolu()
    mapdl.antype("STATIC")
    for node in root_nodes:
        mapdl.d(node, "ALL", 0.0)
    for node in tip_nodes:
        mapdl.f(node, "FZ", -load / len(tip_nodes))
    mapdl.solve()
    mapdl.finish()
    mapdl.post1()
    mapdl.set(1, 1)
    tip_uz = float(np.mean([
        float(mapdl.get_value("NODE", node, "U", "Z")) for node in tip_nodes
    ]))
    mapdl.finish()

    # Mass per panel area is summed from the cured-composite layer densities.
    plies = layup_plies(case_layup)
    areal_mass = sum(
        float(case_layup["material_cards"][mat]["density_kg_m3"])
        * float(case_layup["material_cards"][mat]["ply_thickness_m"])
        for mat, _angle in plies
    )
    planform_area = length * b0 * (1.0 + rb) / 2.0
    return tip_uz, areal_mass * planform_area, total_thickness


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--layup", default="layup.json")
    ap.add_argument("--output", default="results/laminate_shell_sweep.csv")
    ap.add_argument("--nx", type=int, default=40)
    ap.add_argument("--ny", type=int, default=6)
    ap.add_argument("--integration-points", type=int, default=3)
    args = ap.parse_args()
    if args.nx < 4 or args.ny < 2:
        raise ValueError("Use nx >= 4 and ny >= 2 for the initial shell mesh")
    cfg = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text())
    design = make_design(cfg)

    from ansys.mapdl.core import launch_mapdl
    mapdl = launch_mapdl(nproc=int(cfg["fea"].get("nproc", 2)))
    results = []
    try:
        for row in design.itertuples(index=False):
            uz, mass, thickness = solve_case(
                mapdl, row, layup, args.nx, args.ny, args.integration_points)
            case_layup = scaled_layup(
                layup, float(row.thickness_scale), float(row.effective_stiffness_scale))
            d_eff, s_eff = effective_bending_shear_stiffness(case_layup)
            results.append({**row._asdict(), "thickness_m": thickness,
                            "D_eff_Nm": d_eff, "S_eff_N_per_m": s_eff,
                            "tip_deflection_m": uz,
                            "compliance_m_per_N": abs(uz) / float(row.tip_load_N),
                            "mass_kg": mass, "laminate_thickness_m": thickness,
                            "nx": args.nx, "ny": args.ny,
                            "ply_count": len(layup["plies"])})
            print(f"case {row.case_id:03d}: |UZ|={abs(uz):.6e} m, mass={mass:.6g} kg")
    finally:
        mapdl.exit()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(out, index=False)
    print(f"Saved {len(results)} laminate cases to {out}")


if __name__ == "__main__":
    main()
