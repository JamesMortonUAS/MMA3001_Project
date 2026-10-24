"""Run the factorial FEA-to-FEA material-transfer experiment.

The script creates an 80-case training grid and five 20-case target sets,
including a resin-ratio perturbation that is hidden from the prediction model.
It compares mechanics-informed, quadratic and Gaussian-process compliance models.
Two target observations estimate carbon/glass stiffness multipliers. Three
random calibration pairs are reused for every model. Run MAPDL stages with a
licensed local MAPDL installation.

Examples (PowerShell, from the repository root):
    python src/run_extended_transfer_study.py --stage design
    python src/run_extended_transfer_study.py --stage dataset-a
    python src/run_extended_transfer_study.py --stage targets
    python src/run_extended_transfer_study.py --stage fit
    python src/run_extended_transfer_study.py --stage all
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path
from time import perf_counter

import joblib
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.model_selection import KFold
from sklearn.base import clone

BASE_DIR = Path("results/property_transfer_large_deflection_5sets")
ORTHO_FIELDS = {
    "EX_Pa": "EX", "EY_Pa": "EY", "EZ_Pa": "EZ",
    "PRXY": "PRXY", "PRYZ": "PRYZ", "PRXZ": "PRXZ",
    "GXY_Pa": "GXY", "GYZ_Pa": "GYZ", "GXZ_Pa": "GXZ",
    "density_kg_m3": "DENS",
}
ELASTIC_KEYS = ("EX_Pa", "EY_Pa", "EZ_Pa", "GXY_Pa", "GYZ_Pa", "GXZ_Pa")
FEATURES = ["width_m", "total_thickness_m", "D_eff_Nm_per_m", "S_eff_N_per_m"]


def effective_bending_shear_stiffness(layup: dict) -> tuple[float, float]:
    """Return CLT effective D11 [N m] and x-z shear rigidity [N/m].

    This small helper is kept in the runner so the new experiment is runnable
    even when the legacy laminate_section module is absent from a local clone.
    The experiment uses zero-degree plies; non-zero orientations are rejected.
    """
    cards = layup["material_cards"]
    plies = layup["plies"]
    height = sum(float(cards[p["material"]]["ply_thickness_m"]) for p in plies)
    z = -0.5 * height
    abd = np.zeros((6, 6), dtype=float)
    shear_xz = 0.0
    for ply in plies:
        name = str(ply["material"])
        if abs(float(ply.get("orientation_deg", 0.0))) > 1e-12:
            raise NotImplementedError("The extended transfer study currently uses 0-degree plies only")
        card = cards[name]
        t = float(card["ply_thickness_m"])
        e1, e2, nu12, g12 = (float(card[k]) for k in ("EX_Pa", "EY_Pa", "PRXY", "GXY_Pa"))
        nu21 = nu12 * e2 / e1
        den = 1.0 - nu12 * nu21
        if den <= 0.0:
            raise ValueError(f"Invalid in-plane elastic constants for {name}")
        q = np.zeros((3, 3), dtype=float)
        q[0, 0] = e1 / den
        q[1, 1] = e2 / den
        q[0, 1] = q[1, 0] = nu12 * e2 / den
        q[2, 2] = g12
        z_next = z + t
        abd[:3, :3] += q * t
        abd[:3, 3:] += q * (z_next**2 - z**2) / 2.0
        abd[3:, :3] += q * (z_next**2 - z**2) / 2.0
        abd[3:, 3:] += q * (z_next**3 - z**3) / 3.0
        shear_xz += float(card["GXZ_Pa"]) * t
        z = z_next
    compliance = float(np.linalg.inv(abd)[3, 3])
    if compliance <= 0.0 or shear_xz <= 0.0:
        raise ValueError("Laminate stiffness matrix is not positive definite")
    return 1.0 / compliance, shear_xz


def model_inputs(frame: pd.DataFrame) -> np.ndarray:
    """Log-transform positive geometry and section-stiffness predictors."""
    x = frame[FEATURES].to_numpy(float)
    if np.any(~np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("Model geometry/section predictors must be positive and finite")
    return np.log(x)


def validate_config(cfg: dict) -> None:
    """Validate the supported input domain before generating designs or solving.

    Supported runs use positive finite SI dimensions and loads, even ply counts,
    zero-degree plies, the two built-in stack families, and large-deflection
    static MAPDL analysis. The example design is capped at 100 plies.
    """
    required_sections = {"beam", "mesh", "dataset_a", "targets", "calibration", "models", "materials", "solver"}
    missing_sections = required_sections - set(cfg)
    if missing_sections:
        raise ValueError(f"Configuration missing sections: {sorted(missing_sections)}")

    beam = cfg["beam"]
    for key in ("length_m", "tip_load_N", "shear_correction_factor"):
        value = float(beam[key])
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"beam.{key} must be positive and finite")
    mesh = cfg["mesh"]
    for key in ("nx", "ny", "integration_points_per_ply", "nproc"):
        if int(mesh[key]) < 1:
            raise ValueError(f"mesh.{key} must be a positive integer")

    train = cfg["dataset_a"]
    widths = [float(v) for v in train["widths_m"]]
    if len(widths) != 5 or len(set(widths)) != 5 or any(not math.isfinite(v) or v <= 0 for v in widths):
        raise ValueError("dataset_a.widths_m must contain five unique positive finite widths in metres")
    if train["ply_counts"] != [8, 16, 24, 32]:
        raise ValueError("dataset_a.ply_counts must be [8, 16, 24, 32]")
    expected_pairs = [[1.0, 1.0], [1.0, 2.0], [2.0, 1.0], [2.0, 2.0]]
    if train["stiffness_scale_pairs"] != expected_pairs:
        raise ValueError(f"dataset_a.stiffness_scale_pairs must be {expected_pairs}")

    solver = cfg["solver"]
    if solver.get("large_deflection") is not True:
        raise ValueError("solver.large_deflection must be true for this study")
    if solver.get("automatic_time_stepping") is not True:
        raise ValueError("solver.automatic_time_stepping must be true for nonlinear static solving")
    n0, nmin, nmax = (int(solver[k]) for k in ("initial_substeps", "minimum_substeps", "maximum_substeps"))
    if not (1 <= nmin <= n0 <= nmax):
        raise ValueError("Solver substeps must satisfy 1 <= minimum <= initial <= maximum")
    if int(solver["maximum_equilibrium_iterations"]) < 1:
        raise ValueError("solver.maximum_equilibrium_iterations must be positive")

    calibration = cfg["calibration"]
    if int(calibration["repeats"]) < 1 or int(calibration["points_per_repeat"]) < 2:
        raise ValueError("Calibration requires at least one repeat and at least two observations per repeat")
    scale_bounds = [float(v) for v in calibration["material_scale_bounds"]]
    if len(scale_bounds) != 2 or not (0 < scale_bounds[0] < scale_bounds[1]):
        raise ValueError("calibration.material_scale_bounds must be two positive increasing values")

    expected_target_sets = {
        "material_interpolation", "material_extrapolation", "layup_transfer",
        "geometry_extrapolation", "material_extrapolation_v2",
    }
    if set(cfg["targets"]) != expected_target_sets:
        raise ValueError(f"targets must define exactly {sorted(expected_target_sets)}")
    resin_model = cfg.get("resin_model", {})
    resin_required = {"baseline_resin_volume_fraction", "matrix_modulus_Pa", "matrix_poisson_ratio"}
    if resin_required - set(resin_model):
        raise ValueError(f"resin_model missing {sorted(resin_required - set(resin_model))}")
    baseline_resin_fraction = float(resin_model["baseline_resin_volume_fraction"])
    if not math.isfinite(baseline_resin_fraction) or not 0 < baseline_resin_fraction < 1:
        raise ValueError("resin_model.baseline_resin_volume_fraction must be between zero and one")
    matrix_e = float(resin_model["matrix_modulus_Pa"])
    if not math.isfinite(matrix_e) or matrix_e <= 0:
        raise ValueError("resin_model.matrix_modulus_Pa must be positive")
    matrix_nu = float(resin_model["matrix_poisson_ratio"])
    if not math.isfinite(matrix_nu) or not -0.99 < matrix_nu < 0.5:
        raise ValueError("resin_model.matrix_poisson_ratio must be between -0.99 and 0.5")
    for name, target in cfg["targets"].items():
        if target["kind"] not in {"tw_glass", "ud_tw_glass"}:
            raise ValueError(f"targets.{name}.kind is unsupported: {target['kind']}")
        for key in ("carbon_scale", "glass_scale"):
            value = float(target[key])
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"targets.{name}.{key} must be positive and finite")
            if not scale_bounds[0] <= value <= scale_bounds[1]:
                raise ValueError(f"targets.{name}.{key} must lie within calibration.material_scale_bounds")
        if name == "material_extrapolation_v2":
            ratio = float(target.get("resin_to_fibre_ratio_multiplier", 0.0))
            if not math.isfinite(ratio) or ratio <= 0:
                raise ValueError("targets.material_extrapolation_v2.resin_to_fibre_ratio_multiplier must be positive and finite")
        elif float(target.get("resin_to_fibre_ratio_multiplier", 1.0)) != 1.0:
            raise ValueError(f"targets.{name} must use the nominal resin-to-fibre ratio")
        if name == "geometry_extrapolation":
            geo_widths = [float(v) for v in target["widths_m"]]
            geo_plies = [int(v) for v in target["ply_counts"]]
            if any(float(raw) != cast for raw, cast in zip(target["ply_counts"], geo_plies)):
                raise ValueError("geometry_extrapolation.ply_counts must be integer values")
            if len(geo_widths) != 5 or any(not math.isfinite(v) or v <= 0 for v in geo_widths):
                raise ValueError("geometry_extrapolation.widths_m must contain five positive finite widths")
            if any(v < 2 or v > 100 or v % 2 for v in geo_plies):
                raise ValueError("geometry_extrapolation.ply_counts must be even integers from 2 to 100")

    for name, card in cfg["materials"].items():
        required = {"ply_thickness_m", "family", *ORTHO_FIELDS}
        missing = required - set(card)
        if missing:
            raise ValueError(f"Material {name} is missing {sorted(missing)}")
        if card["family"] not in {"carbon", "glass"}:
            raise ValueError(f"Material {name}.family must be carbon or glass")
        for key in ("ply_thickness_m", "density_kg_m3", "EX_Pa", "EY_Pa", "EZ_Pa", "GXY_Pa", "GYZ_Pa", "GXZ_Pa"):
            value = float(card[key])
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"Material {name}.{key} must be positive and finite")
        for key in ("PRXY", "PRYZ", "PRXZ"):
            value = float(card[key])
            if not math.isfinite(value) or not -0.99 < value < 0.99:
                raise ValueError(f"Material {name}.{key} must be finite and between -0.99 and 0.99")


def symmetric_stack(kind: str, ply_count: int) -> list[str]:
    """Build a symmetric, bottom-to-top stack from a repeating half-stack."""
    if ply_count <= 0 or ply_count % 2:
        raise ValueError("Symmetric stacks require a positive, even ply count")
    if kind == "tw_glass":
        cycle = ["CF_TWILL", "GF"]
    elif kind == "ud_tw_glass":
        cycle = ["CF_UD", "CF_TWILL", "GF"]
    else:
        raise ValueError(f"Unsupported stack kind: {kind}")
    half = [cycle[i % len(cycle)] for i in range(ply_count // 2)]
    full = half + list(reversed(half))
    if full != list(reversed(full)) or len(full) != ply_count:
        raise AssertionError("Generated ply sequence is not symmetric")
    return full


def make_designs(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return training and target design tables with actual stack metadata."""
    a_rows = []
    cid = 0
    a = cfg["dataset_a"]
    for width, nplies, pair in itertools.product(
        a["widths_m"], a["ply_counts"], a["stiffness_scale_pairs"]
    ):
        stack = symmetric_stack("tw_glass", int(nplies))
        a_rows.append({
            "case_id": cid, "width_m": float(width), "ply_count": int(nplies),
            "carbon_scale": float(pair[0]), "glass_scale": float(pair[1]),
            "target_set": "training", "stack_kind": "tw_glass",
            "ply_sequence_bottom_to_top": "-".join(stack),
        })
        cid += 1

    target_rows = []
    set_specs = cfg["targets"]
    for set_name in ("material_interpolation", "material_extrapolation", "layup_transfer"):
        spec = set_specs[set_name]
        for nplies, width in itertools.product(a["ply_counts"], a["widths_m"]):
            stack = symmetric_stack(spec["kind"], int(nplies))
            target_rows.append({
                "case_id": len(target_rows), "target_set": set_name,
                "width_m": float(width), "ply_count": int(nplies),
                "carbon_scale_truth": float(spec["carbon_scale"]),
                "glass_scale_truth": float(spec["glass_scale"]),
                "resin_to_fibre_ratio_multiplier_truth": float(spec.get("resin_to_fibre_ratio_multiplier", 1.0)),
                "stack_kind": spec["kind"],
                "ply_sequence_bottom_to_top": "-".join(stack),
            })
    spec = set_specs["geometry_extrapolation"]
    training_w = (min(a["widths_m"]), max(a["widths_m"]))
    training_n = (min(a["ply_counts"]), max(a["ply_counts"]))
    for nplies, width in itertools.product(spec["ply_counts"], spec["widths_m"]):
        stack = symmetric_stack(spec["kind"], int(nplies))
        target_rows.append({
            "case_id": len(target_rows), "target_set": "geometry_extrapolation",
            "width_m": float(width), "ply_count": int(nplies),
            "carbon_scale_truth": float(spec["carbon_scale"]),
            "glass_scale_truth": float(spec["glass_scale"]),
            "resin_to_fibre_ratio_multiplier_truth": float(spec.get("resin_to_fibre_ratio_multiplier", 1.0)),
            "stack_kind": spec["kind"],
            "ply_sequence_bottom_to_top": "-".join(stack),
            "width_region": "in_range" if training_w[0] <= width <= training_w[1] else "extrapolation",
            "ply_count_region": "in_range" if training_n[0] <= nplies <= training_n[1] else "extrapolation",
        })
    spec = set_specs["material_extrapolation_v2"]
    for nplies, width in itertools.product(a["ply_counts"], a["widths_m"]):
        stack = symmetric_stack(spec["kind"], int(nplies))
        target_rows.append({
            "case_id": len(target_rows), "target_set": "material_extrapolation_v2",
            "width_m": float(width), "ply_count": int(nplies),
            "carbon_scale_truth": float(spec["carbon_scale"]),
            "glass_scale_truth": float(spec["glass_scale"]),
            "resin_to_fibre_ratio_multiplier_truth": float(spec["resin_to_fibre_ratio_multiplier"]),
            "stack_kind": spec["kind"],
            "ply_sequence_bottom_to_top": "-".join(stack),
        })
    return pd.DataFrame(a_rows), pd.DataFrame(target_rows)


def _adjust_cards_for_resin_ratio(cfg: dict, resin_ratio_multiplier: float) -> dict:
    """Apply an illustrative matrix-rich modulus proxy at fixed ply thickness.

    The configured ply cards are treated as effective properties at the
    baseline resin fraction. Their notional fibre contribution is inferred
    with a Voigt mixture, then remixed at the requested resin:fibre ratio.
    Density and Poisson ratios remain fixed; this is a compliance-study proxy,
    not a physical resin-uptake or mass model.
    """
    factor = float(resin_ratio_multiplier)
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError("resin-to-fibre ratio multiplier must be positive and finite")
    resin_cfg = cfg["resin_model"]
    base_resin = float(resin_cfg["baseline_resin_volume_fraction"])
    new_resin = factor * base_resin / (1.0 - base_resin + factor * base_resin)
    matrix_e = float(resin_cfg["matrix_modulus_Pa"])
    matrix_g = matrix_e / (2.0 * (1.0 + float(resin_cfg["matrix_poisson_ratio"])))
    cards = copy.deepcopy(cfg["materials"])
    for card in cards.values():
        for key in ELASTIC_KEYS:
            matrix_value = matrix_g if key.startswith("G") else matrix_e
            baseline = float(card[key])
            notional_fibre = (baseline - base_resin * matrix_value) / (1.0 - base_resin)
            adjusted = (1.0 - new_resin) * notional_fibre + new_resin * matrix_value
            if not math.isfinite(adjusted) or adjusted <= 0:
                raise ValueError(f"Resin mixture produced invalid {key} for a material card")
            card[key] = adjusted
    return cards


def build_layup(cfg: dict, kind: str, ply_count: int,
                carbon_scale: float = 1.0, glass_scale: float = 1.0,
                resin_ratio_multiplier: float = 1.0) -> dict:
    """Make a generic layup object accepted by the laminate property helper."""
    cards = _adjust_cards_for_resin_ratio(cfg, resin_ratio_multiplier)
    for name, card in cards.items():
        factor = carbon_scale if card["family"] == "carbon" else glass_scale
        for key in ELASTIC_KEYS:
            card[key] = float(card[key]) * float(factor)
    stack = symmetric_stack(kind, int(ply_count))
    return {"material_cards": cards, "plies": [
        {"material": name, "orientation_deg": 0.0} for name in stack
    ]}


def section_data(cfg: dict, row: dict, carbon_scale: float | None = None,
                 glass_scale: float | None = None,
                 resin_ratio_multiplier: float | None = None) -> dict:
    """Calculate total thickness, mass and CLT section stiffnesses."""
    cscale = float(row.get("carbon_scale_truth", row.get("carbon_scale", 1.0)) if carbon_scale is None else carbon_scale)
    gscale = float(row.get("glass_scale_truth", row.get("glass_scale", 1.0)) if glass_scale is None else glass_scale)
    resin_multiplier = float(
        row.get("resin_to_fibre_ratio_multiplier_truth", 1.0)
        if resin_ratio_multiplier is None else resin_ratio_multiplier
    )
    kind = str(row["stack_kind"])
    layup = build_layup(cfg, kind, int(row["ply_count"]), cscale, gscale, resin_multiplier)
    d_eff, s_eff = effective_bending_shear_stiffness(layup)
    cards = layup["material_cards"]
    thickness = sum(float(cards[p["material"]]["ply_thickness_m"]) for p in layup["plies"])
    mass_per_area = sum(
        float(cards[p["material"]]["density_kg_m3"]) * float(cards[p["material"]]["ply_thickness_m"])
        for p in layup["plies"]
    )
    width = float(row["width_m"])
    kappa = float(cfg["beam"]["shear_correction_factor"])
    return {
        "total_thickness_m": thickness,
        "D_eff_Nm_per_m": d_eff,
        "S_eff_N_per_m": s_eff,
        "EI_Nm2": width * d_eff,
        "GA_N": width * s_eff,
        "kGA_N": kappa * width * s_eff,
        "mass_kg": mass_per_area * width * float(cfg["beam"]["length_m"]),
        "carbon_scale": cscale,
        "glass_scale": gscale,
        "resin_to_fibre_ratio_multiplier": resin_multiplier,
    }


def _node_id(i: int, j: int, ny: int) -> int:
    """Return the deterministic MAPDL node number for grid indices ``i, j``."""
    return i * (ny + 1) + j + 1


def configure_large_deflection(mapdl, cfg: dict) -> None:
    """Activate nonlinear geometric static solution controls in MAPDL."""
    solver = cfg["solver"]
    if solver.get("large_deflection") is not True:
        raise ValueError("This runner requires solver.large_deflection=true")
    mapdl.run("NLGEOM,ON")
    mapdl.run("AUTOTS,ON")
    mapdl.run(
        f"NSUBST,{int(solver['initial_substeps'])},"
        f"{int(solver['minimum_substeps'])},{int(solver['maximum_substeps'])}"
    )
    mapdl.run(f"NEQIT,{int(solver['maximum_equilibrium_iterations'])}")
    mapdl.run("NROPT,FULL")


def _case_signature(cfg: dict, row: dict) -> str:
    """Hash the full FEA definition and case design to prevent stale resume."""
    payload = {
        "beam": cfg["beam"], "mesh": cfg["mesh"], "solver": cfg["solver"],
        "resin_model": cfg["resin_model"],
        "materials": cfg["materials"], "case": row,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
        default=lambda value: value.item() if isinstance(value, np.generic) else str(value),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def solve_shell_case(mapdl, cfg: dict, row: dict) -> dict:
    """Solve one generic symmetric layered shell cantilever in MAPDL."""
    started = perf_counter()
    cscale = float(row.get("carbon_scale_truth", row.get("carbon_scale", 1.0)))
    gscale = float(row.get("glass_scale_truth", row.get("glass_scale", 1.0)))
    resin_ratio = float(row.get("resin_to_fibre_ratio_multiplier_truth", 1.0))
    cards = _adjust_cards_for_resin_ratio(cfg, resin_ratio)
    for card in cards.values():
        factor = cscale if card["family"] == "carbon" else gscale
        for key in ELASTIC_KEYS:
            card[key] = float(card[key]) * factor
    stack = str(row["ply_sequence_bottom_to_top"]).split("-")
    names = list(cards)
    mat_ids = {name: i + 1 for i, name in enumerate(names)}
    mesh = cfg["mesh"]
    nx, ny = int(mesh["nx"]), int(mesh["ny"])
    load = float(cfg["beam"]["tip_load_N"])
    length = float(cfg["beam"]["length_m"])
    width = float(row["width_m"])
    mapdl.clear()
    mapdl.prep7()
    for name, mid in mat_ids.items():
        for prop, label in ORTHO_FIELDS.items():
            mapdl.mp(label, mid, float(cards[name][prop]))
    mapdl.et(1, "SHELL181")
    mapdl.keyopt(1, 8, 1)
    mapdl.sectype(1, "SHELL")
    thickness = 0.0
    for material in stack:
        card = cards[material]
        t = float(card["ply_thickness_m"])
        mapdl.secdata(t, mat_ids[material], 0.0, int(mesh["integration_points_per_ply"]))
        thickness += t
    for i in range(nx + 1):
        x = length * i / nx
        for j in range(ny + 1):
            y = -0.5 * width + width * j / ny
            mapdl.n(_node_id(i, j, ny), x, y, 0.0)
    mapdl.type(1)
    mapdl.secnum(1)
    for i in range(nx):
        for j in range(ny):
            mapdl.e(_node_id(i, j, ny), _node_id(i + 1, j, ny),
                    _node_id(i + 1, j + 1, ny), _node_id(i, j + 1, ny))
    roots = [_node_id(0, j, ny) for j in range(ny + 1)]
    tips = [_node_id(nx, j, ny) for j in range(ny + 1)]
    mapdl.finish()
    mapdl.slashsolu()
    mapdl.antype("STATIC")
    configure_large_deflection(mapdl, cfg)
    for node in roots:
        mapdl.d(node, "ALL", 0.0)
    for node in tips:
        mapdl.f(node, "FZ", -load / len(tips))
    mapdl.solve()
    mapdl.finish()
    mapdl.post1()
    mapdl.set(1, 1)
    uz = float(np.mean([float(mapdl.get_value("NODE", n, "U", "Z")) for n in tips]))
    mapdl.finish()
    props = section_data(cfg, row, cscale, gscale, resin_ratio)
    return {
        **props, "tip_deflection_m": uz,
        "compliance_m_per_N": abs(uz) / load,
        "max_deflection_over_length": abs(uz) / length,
        "tip_load_N": load, "nx": nx, "ny": ny,
        "analysis_formulation": "large_deflection_static",
        "fea_runtime_s": perf_counter() - started,
        "ply_count": len(stack),
        "ply_sequence_bottom_to_top": "-".join(stack),
    }


def run_fea(cfg: dict, design: pd.DataFrame, output: Path,
            dataset: str, max_cases: int | None = None) -> None:
    """Run/resume a dataset and persist each completed case immediately."""
    from ansys.mapdl.core import launch_mapdl
    validate_config(cfg)
    output.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if output.exists():
        prior = pd.read_csv(output)
        if "case_id" in prior:
            if "analysis_signature" not in prior or "analysis_formulation" not in prior:
                raise ValueError(
                    f"Refusing to resume {output}: it lacks the large-deflection analysis signature. "
                    "Use a new output directory; old linear-static results cannot be reused."
                )
            expected = {
                int(row["case_id"]): _case_signature(cfg, row)
                for row in design.to_dict(orient="records")
            }
            for prior_row in prior.to_dict(orient="records"):
                cid = int(prior_row["case_id"])
                if (prior_row["analysis_formulation"] != "large_deflection_static"
                        or prior_row["analysis_signature"] != expected.get(cid)):
                    raise ValueError(
                        f"Refusing to resume {output}: case {cid} was generated with a different "
                        "FEA definition. Choose a fresh output directory."
                    )
            done = {int(r["case_id"]): r for r in prior.to_dict(orient="records")}
    mapdl = launch_mapdl(nproc=int(cfg["mesh"].get("nproc", 2)))
    try:
        run_design = design if max_cases is None else design.head(max_cases)
        for row in run_design.to_dict(orient="records"):
            cid = int(row["case_id"])
            if cid in done:
                continue
            solved = solve_shell_case(mapdl, cfg, row)
            record = {**row, **solved, "analysis_signature": _case_signature(cfg, row)}
            done[cid] = record
            pd.DataFrame(done.values()).sort_values("case_id").to_csv(output, index=False)
            print(f"{dataset} case {cid:03d}: width={row['width_m']:.4f} m, "
                  f"plies={row['ply_count']}, h={solved['total_thickness_m']*1e3:.3f} mm, "
                  f"C={solved['compliance_m_per_N']:.6e} m/N, "
                  f"time={solved['fea_runtime_s']:.2f} s")
    finally:
        mapdl.exit()
    print(f"Saved {len(done)} {dataset} cases to {output}")


def mechanics_matrix(frame: pd.DataFrame, length: float) -> np.ndarray:
    """Build bending and shear compliance basis columns for each beam case."""
    return np.column_stack((
        length**3 / (3 * frame["EI_Nm2"].to_numpy(float)),
        length / frame["kGA_N"].to_numpy(float),
    ))


def fit_models(data_a: pd.DataFrame, cfg: dict,
               runtime_seconds: dict[str, float] | None = None) -> dict:
    """Fit all model families on Dataset A and optionally record wall time."""
    y = data_a["compliance_m_per_N"].to_numpy(float)
    if np.any(y <= 0):
        raise ValueError("Compliance must be positive")
    length = float(cfg["beam"]["length_m"])
    started = perf_counter()
    mech = LinearRegression(fit_intercept=False, positive=True).fit(mechanics_matrix(data_a, length), y)
    if runtime_seconds is not None:
        runtime_seconds["mechanics"] = perf_counter() - started
    x = model_inputs(data_a)
    logy = np.log(y)
    started = perf_counter()
    quad = make_pipeline(StandardScaler(), PolynomialFeatures(degree=2, include_bias=False),
                         Ridge(alpha=float(cfg["models"].get("quadratic_ridge_alpha", 1e-8))))
    quad.fit(x, logy)
    if runtime_seconds is not None:
        runtime_seconds["quadratic"] = perf_counter() - started
    started = perf_counter()
    kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(length_scale=np.ones(len(FEATURES)), nu=2.5) + WhiteKernel(noise_level=1e-5)
    gp = make_pipeline(StandardScaler(), GaussianProcessRegressor(
        kernel=kernel, normalize_y=True, n_restarts_optimizer=int(cfg["models"].get("gp_restarts", 2)),
        random_state=int(cfg["calibration"]["seed"])))
    gp.fit(x, logy)
    if runtime_seconds is not None:
        runtime_seconds["gaussian_process"] = perf_counter() - started
    return {"mechanics": mech, "quadratic": quad, "gaussian_process": gp}


def cross_validated_training_metrics(data_a: pd.DataFrame, cfg: dict,
                                     models: dict,
                                     runtime_seconds: dict[str, float] | None = None) -> dict:
    """Five-fold CV predictions on Dataset A; no target data enters model selection."""
    folds = KFold(n_splits=5, shuffle=True, random_state=int(cfg["calibration"]["seed"]))
    y = data_a["compliance_m_per_N"].to_numpy(float)
    results = {}
    for name, prototype in models.items():
        started = perf_counter()
        pred = np.empty_like(y)
        for train_ix, test_ix in folds.split(data_a):
            fitted = clone(prototype)
            train = data_a.iloc[train_ix]
            test = data_a.iloc[test_ix]
            if name == "mechanics":
                fitted.fit(mechanics_matrix(train, float(cfg["beam"]["length_m"])),
                           train["compliance_m_per_N"].to_numpy(float))
                pred[test_ix] = np.maximum(fitted.predict(
                    mechanics_matrix(test, float(cfg["beam"]["length_m"]))), 1e-15)
            else:
                fitted.fit(model_inputs(train),
                           np.log(train["compliance_m_per_N"].to_numpy(float)))
                pred[test_ix] = np.exp(np.clip(fitted.predict(model_inputs(test)), -50, 50))
        results[name] = metrics(y, pred)
        if runtime_seconds is not None:
            runtime_seconds[name] = perf_counter() - started
    return results


def predict(model_name: str, model, frame: pd.DataFrame, cfg: dict) -> np.ndarray:
    """Predict positive compliance using a fitted model and case properties."""
    if model_name == "mechanics":
        return np.maximum(model.predict(mechanics_matrix(frame, float(cfg["beam"]["length_m"]))), 1e-15)
    return np.exp(np.clip(model.predict(model_inputs(frame)), -50, 50))


def metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    """Return sample count, absolute-error scores, MAPE and coefficient of determination."""
    ape = 100 * np.abs(pred-y) / np.maximum(np.abs(y), 1e-15)
    return {
        "n": int(len(y)), "rmse_m_per_N": float(np.sqrt(mean_squared_error(y, pred))),
        "mae_m_per_N": float(mean_absolute_error(y, pred)),
        "mape_percent": float(np.mean(ape)), "max_ape_percent": float(np.max(ape)),
        "r2": float(r2_score(y, pred)) if len(y) > 1 else None,
    }


def generate_features_for_scales(cfg: dict, target_row: dict,
                                 cscale: float, gscale: float) -> dict:
    """Recompute a target case's effective section predictors at candidate scales."""
    row = {**target_row, "carbon_scale_truth": cscale, "glass_scale_truth": gscale}
    row["resin_to_fibre_ratio_multiplier_truth"] = 1.0
    props = section_data(cfg, row, cscale, gscale, resin_ratio_multiplier=1.0)
    return {**target_row, **props}


def fit_evaluate(cfg: dict, a_path: Path, target_path: Path, out_dir: Path) -> dict:
    """Fit on Dataset A and evaluate repeated two-case calibration per target set."""
    data_a = pd.read_csv(a_path)
    target_data = pd.read_csv(target_path)
    fit_runtime = {}
    models = fit_models(data_a, cfg, fit_runtime)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, model in models.items():
        joblib.dump(model, out_dir / f"model_{name}.joblib")
    train_metrics = {}
    for name, model in models.items():
        train_metrics[name] = metrics(data_a["compliance_m_per_N"].to_numpy(float),
                                      predict(name, model, data_a, cfg))
    cv_runtime = {}
    cv_metrics = cross_validated_training_metrics(data_a, cfg, models, cv_runtime)

    rng = np.random.default_rng(int(cfg["calibration"]["seed"]))
    n_repeats = int(cfg["calibration"]["repeats"])
    n_cal = int(cfg["calibration"]["points_per_repeat"])
    bounds = cfg["calibration"]["material_scale_bounds"]
    result_rows, metric_rows, split_rows, target_runtime_rows = [], [], [], []
    for set_name, group in target_data.groupby("target_set", sort=False):
        group = group.sort_values("case_id").reset_index(drop=True)
        if len(group) <= n_cal:
            raise ValueError(f"Target set {set_name} too small for calibration")
        for repeat in range(n_repeats):
            chosen = np.sort(rng.choice(len(group), size=n_cal, replace=False))
            cal_ids = group.iloc[chosen]["case_id"].astype(int).tolist()
            split_rows.append({"target_set": set_name, "repeat": repeat + 1,
                               "calibration_case_ids": ",".join(map(str, cal_ids))})
            y_all = group["compliance_m_per_N"].to_numpy(float)
            cal_mask = np.zeros(len(group), dtype=bool)
            cal_mask[chosen] = True
            for model_name, model in models.items():
                evaluation_started = perf_counter()
                nominal_frame = pd.DataFrame([
                    generate_features_for_scales(cfg, r, 1.0, 1.0)
                    for r in group.to_dict(orient="records")
                ])
                pred_nominal = predict(model_name, model, nominal_frame, cfg)

                def residual(params):
                    xframe = pd.DataFrame([
                        generate_features_for_scales(cfg, group.iloc[idx].to_dict(), params[0], params[1])
                        for idx in chosen
                    ])
                    p = predict(model_name, model, xframe, cfg)
                    return np.log(p) - np.log(y_all[chosen])

                initial = np.array([1.0, 1.0])
                fit = least_squares(residual, initial, bounds=(float(bounds[0]), float(bounds[1])),
                                    x_scale="jac", max_nfev=1000)
                q = fit.x
                jac = np.asarray(fit.jac, float)
                rank = int(np.linalg.matrix_rank(jac))
                singular = np.linalg.svd(jac, compute_uv=False)
                cond = float(singular[0] / singular[-1]) if len(singular) > 1 and singular[-1] > 0 else float("inf")
                fitted_frame = pd.DataFrame([
                    generate_features_for_scales(cfg, r, q[0], q[1])
                    for r in group.to_dict(orient="records")
                ])
                pred_fit = predict(model_name, model, fitted_frame, cfg)
                target_runtime_rows.append({
                    "target_set": set_name, "repeat": repeat + 1, "model": model_name,
                    "calibration_and_prediction_runtime_s": perf_counter() - evaluation_started,
                })
                valid = ~cal_mask
                region_masks = {"all_validation": valid}
                if "width_region" in group:
                    for region_name in group["width_region"].unique():
                        region_masks[f"width_{region_name}"] = valid & (group["width_region"].to_numpy() == region_name)
                if "ply_count_region" in group:
                    for region_name in group["ply_count_region"].unique():
                        region_masks[f"ply_count_{region_name}"] = valid & (group["ply_count_region"].to_numpy() == region_name)
                for mode, pred in (("nominal", pred_nominal), ("calibrated", pred_fit)):
                    for region, region_mask in region_masks.items():
                        if not np.any(region_mask):
                            continue
                        stat = metrics(y_all[region_mask], pred[region_mask])
                        metric_rows.append({
                            "target_set": set_name, "repeat": repeat + 1,
                            "model": model_name, "mode": mode, "evaluation_region": region,
                            "calibration_case_ids": ",".join(map(str, cal_ids)),
                            "carbon_scale_fit": float(q[0]), "glass_scale_fit": float(q[1]),
                            "jacobian_rank": rank, "jacobian_condition_number": cond,
                            "ill_conditioned": bool(rank < 2 or cond > 1e6), **stat,
                        })
                for ix, row in group.iterrows():
                    result_rows.append({
                        "target_set": set_name, "case_id": int(row["case_id"]),
                        "repeat": repeat + 1, "model": model_name,
                        "is_calibration": bool(cal_mask[ix]),
                        "observed_compliance_m_per_N": float(y_all[ix]),
                        "nominal_prediction_m_per_N": float(pred_nominal[ix]),
                        "calibrated_prediction_m_per_N": float(pred_fit[ix]),
                        "carbon_scale_fit": float(q[0]), "glass_scale_fit": float(q[1]),
                        "jacobian_rank": rank, "jacobian_condition_number": cond,
                        "width_region": row.get("width_region", "within_training"),
                        "ply_count_region": row.get("ply_count_region", "within_training"),
                    })

    pd.DataFrame(result_rows).to_csv(out_dir / "predictions.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(out_dir / "metrics_by_repeat.csv", index=False)
    pd.DataFrame(split_rows).to_csv(out_dir / "calibration_splits.csv", index=False)
    target_runtime = pd.DataFrame(target_runtime_rows)
    runtime_summary = target_runtime.groupby("model")["calibration_and_prediction_runtime_s"].agg(
        calibration_prediction_mean_s="mean", calibration_prediction_sd_s="std"
    ).fillna(0.0)
    runtime_table = pd.DataFrame({
        "model": list(models),
        "fit_runtime_s": [fit_runtime[name] for name in models],
        "five_fold_cv_runtime_s": [cv_runtime[name] for name in models],
    }).set_index("model").join(runtime_summary).reset_index()
    runtime_table.to_csv(out_dir / "runtime_by_model.csv", index=False)
    target_runtime.to_csv(out_dir / "runtime_by_repeat.csv", index=False)
    summary = {
        "training_model_metrics": train_metrics,
        "training_fivefold_cv_metrics": cv_metrics,
        "runtime_by_model": runtime_table.to_dict(orient="records"),
        "fea_runtime_column": "fea_runtime_s (seconds per case; recorded in dataset_A.csv and targets.csv)",
        "experiment": {"dataset_a_cases": int(len(data_a)), "target_cases": int(len(target_data)),
                       "calibration_repeats": n_repeats, "calibration_points": n_cal},
        "calibration_warning": "Two points are the algebraic minimum for two factors. Inspect Jacobian rank and condition number before interpreting factors.",
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    """Parse the CLI, validate inputs, and run selected design/FEA/model stages."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.extended_transfer.example.json")
    parser.add_argument("--output-dir", default=str(BASE_DIR))
    parser.add_argument("--stage", choices=("design", "dataset-a", "targets", "fit", "all"), default="all")
    parser.add_argument("--max-cases", type=int, default=None,
                        help="Run only the first N cases for the selected FEA stage (smoke test).")
    args = parser.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    validate_config(cfg)
    if args.max_cases is not None and args.max_cases < 1:
        parser.error("--max-cases must be at least 1")
    out = Path(args.output_dir)
    a_design, target_design = make_designs(cfg)
    a_design_path = out / "dataset_A_design.csv"
    target_design_path = out / "target_design.csv"
    if args.stage in ("design", "all"):
        out.mkdir(parents=True, exist_ok=True)
        a_design.to_csv(a_design_path, index=False)
        target_design.to_csv(target_design_path, index=False)
        print(f"Wrote Dataset A design ({len(a_design)} cases), target design ({len(target_design)} cases).")
    if args.stage in ("dataset-a", "all"):
        run_fea(cfg, a_design, out / "dataset_A.csv", "Dataset A", args.max_cases)
    if args.stage in ("targets", "all"):
        run_fea(cfg, target_design, out / "targets.csv", "Targets", args.max_cases)
    if args.stage in ("fit", "all"):
        result = fit_evaluate(cfg, out / "dataset_A.csv", out / "targets.csv", out / "model")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
