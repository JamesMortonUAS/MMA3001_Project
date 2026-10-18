"""Run a structured SHELL181 mesh-convergence study on representative beams.

Run this on a machine with Ansys MAPDL and PyMAPDL available. The script uses
the same material, support, and tip-load setup as the main laminate sweep.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analytic import tip_compliance_from_stiffness
from laminate_section import effective_bending_shear_stiffness, scaled_layup
from run_laminate_shell_sweep import solve_case


DEFAULT_MESHES = ((20, 3), (40, 6), (80, 12), (160, 24))


def parse_meshes(values: list[str]) -> list[tuple[int, int]]:
    """Parse mesh tokens such as ``40x6`` into (length, width) divisions."""
    meshes = []
    for value in values:
        try:
            nx_text, ny_text = value.lower().split("x", maxsplit=1)
            nx, ny = int(nx_text), int(ny_text)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Invalid mesh {value!r}; expected NXxNY, e.g. 40x6") from exc
        if nx < 4 or ny < 2:
            raise ValueError(f"Mesh {value!r} is too coarse; require NX >= 4 and NY >= 2")
        meshes.append((nx, ny))
    if len(meshes) < 2:
        raise ValueError("At least two mesh levels are required")
    if len(set(meshes)) != len(meshes):
        raise ValueError("Mesh levels must be unique")
    if meshes != sorted(meshes, key=lambda mesh: mesh[0] * mesh[1]):
        raise ValueError("List mesh levels from coarsest to finest")
    return meshes


def representative_cases(config: dict) -> list[dict]:
    """Return baseline and two geometry-extreme cases at the shorter length."""
    base = config["baseline"]
    length = float(base["length_m"])
    width = float(base["root_width_m"])
    load = float(base.get("tip_load_N", 10.0))
    # Same laminate and load in every case; vary geometry to check transfer.
    return [
        {"case": "baseline", "length_m": length, "root_width_m": width,
         "taper_ratio": float(base["taper_ratio"])},
        {"case": "narrow_high_taper", "length_m": length,
         "root_width_m": 0.8 * width, "taper_ratio": 0.4},
        {"case": "wide_low_taper", "length_m": length,
         "root_width_m": 1.2 * width, "taper_ratio": 0.9},
    ]


def add_convergence_columns(results: pd.DataFrame) -> pd.DataFrame:
    """Add successive and finest-mesh relative compliance changes by case."""
    results = results.sort_values(["case", "nx", "ny"]).copy()
    results["change_from_previous_percent"] = results.groupby("case")[
        "compliance_m_per_N"].pct_change() * 100.0
    finest = results.groupby("case")["compliance_m_per_N"].transform("last")
    results["difference_from_finest_percent"] = (
        (results["compliance_m_per_N"] - finest).abs() / finest * 100.0
    )
    return results


def write_plot(results: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    for case, group in results.groupby("case", sort=False):
        ax.plot(group["element_count"], group["compliance_m_per_N"],
                marker="o", label=case.replace("_", " "))
    ax.set_xscale("log")
    ax.set_xlabel("SHELL181 element count")
    ax.set_ylabel("Tip compliance |UZ| / P (m/N)")
    ax.set_title("Mesh convergence at L = 0.25 m")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.quadratic.example.json")
    parser.add_argument("--layup", default="layup.json")
    parser.add_argument("--output", default="results/mesh_convergence.csv")
    parser.add_argument("--meshes", nargs="+", default=["20x3", "40x6", "80x12", "160x24"])
    parser.add_argument("--integration-points", type=int, default=3)
    args = parser.parse_args()

    meshes = parse_meshes(args.meshes)
    config = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text())
    cases = representative_cases(config)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    from ansys.mapdl.core import launch_mapdl

    mapdl = launch_mapdl(nproc=int(config["fea"].get("nproc", 2)))
    rows = []
    try:
        for case in cases:
            design_row = SimpleNamespace(
                length_m=case["length_m"],
                root_width_m=case["root_width_m"],
                taper_ratio=case["taper_ratio"],
                tip_load_N=float(config["baseline"].get("tip_load_N", 10.0)),
                thickness_scale=1.0,
                effective_stiffness_scale=1.0,
            )
            for nx, ny in meshes:
                uz, mass, thickness = solve_case(
                    mapdl, design_row, layup, nx, ny, args.integration_points)
                d_eff, s_eff = effective_bending_shear_stiffness(
                    scaled_layup(layup, 1.0, 1.0))
                compliance = abs(uz) / design_row.tip_load_N
                analytic = tip_compliance_from_stiffness(
                    design_row.length_m, design_row.root_width_m,
                    design_row.taper_ratio, d_eff, s_eff,
                    float(config["fea"].get("shear_correction_factor_analytic", 5 / 6)))
                rows.append({
                    **case,
                    "tip_load_N": design_row.tip_load_N,
                    "thickness_scale": 1.0,
                    "effective_stiffness_scale": 1.0,
                    "nx": nx,
                    "ny": ny,
                    "element_count": nx * ny,
                    "node_count": (nx + 1) * (ny + 1),
                    "integration_points_per_ply": args.integration_points,
                    "ply_count": len(layup["plies"]),
                    "laminate_thickness_m": thickness,
                    "mass_kg": mass,
                    "tip_deflection_m": uz,
                    "compliance_m_per_N": compliance,
                    "analytical_compliance_m_per_N": analytic,
                    "fea_vs_analytical_percent": 100 * (compliance - analytic) / analytic,
                    "max_deflection_over_length": abs(uz) / design_row.length_m,
                })
                print(f"{case['case']}: {nx}x{ny}, |UZ|={abs(uz):.6e} m, "
                      f"C={compliance:.6e} m/N")
    finally:
        mapdl.exit()

    results = add_convergence_columns(pd.DataFrame(rows))
    results.to_csv(output, index=False)
    plot_path = output.with_name(f"{output.stem}_compliance.png")
    write_plot(results, plot_path)
    print(f"Saved results: {output}")
    print(f"Saved plot: {plot_path}")
    print("Finest-step relative changes (%):")
    print(results.groupby("case").tail(2)[["case", "nx", "ny", "change_from_previous_percent"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
