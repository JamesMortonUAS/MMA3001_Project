"""Run the factorial FEA-to-FEA material-transfer experiment.

The script creates an 80-case training grid, four 20-case target sets, and
compares mechanics-informed, quadratic and Gaussian-process compliance models.
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
import itertools
import json
import math
from pathlib import Path

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

from laminate_section import effective_bending_shear_stiffness


BASE_DIR = Path("results/property_transfer_extended")
ORTHO_FIELDS = {
    "EX_Pa": "EX", "EY_Pa": "EY", "EZ_Pa": "EZ",
    "PRXY": "PRXY", "PRYZ": "PRYZ", "PRXZ": "PRXZ",
    "GXY_Pa": "GXY", "GYZ_Pa": "GYZ", "GXZ_Pa": "GXZ",
    "density_kg_m3": "DENS",
}
ELASTIC_KEYS = ("EX_Pa", "EY_Pa", "EZ_Pa", "GXY_Pa", "GYZ_Pa", "GXZ_Pa")
FEATURES = ["width_m", "total_thickness_m", "D_eff_Nm_per_m", "S_eff_N_per_m"]


def model_inputs(frame: pd.DataFrame) -> np.ndarray:
    """Log-transform positive geometry and section-stiffness predictors."""
    x = frame[FEATURES].to_numpy(float)
    if np.any(~np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("Model geometry/section predictors must be positive and finite")
    return np.log(x)


def validate_config(cfg: dict) -> None:
    """Check the factors that determine sweep size and valid stacks."""
    if len(cfg["dataset_a"]["widths_m"]) != 5:
        raise ValueError("Dataset A requires five width levels")
    if cfg["dataset_a"]["ply_counts"] != [8, 16, 24, 32]:
        raise ValueError("Dataset A ply counts must be [8, 16, 24, 32]")
    if len(cfg["dataset_a"]["stiffness_scale_pairs"]) != 4:
        raise ValueError("Dataset A requires four carbon/glass scale pairs")
    for name, card in cfg["materials"].items():
        required = {"ply_thickness_m", "family", *ORTHO_FIELDS}
        missing = required - set(card)
        if missing:
            raise ValueError(f"Material {name} is missing {sorted(missing)}")
        if float(card["ply_thickness_m"]) <= 0:
            raise ValueError(f"Material {name} ply thickness must be positive")


def symmetric_stack(kind: str, ply_count: int) -> list[str]:
    """Build a symmetric, bottom-to-top stack from a repeating half-stack."""
    if ply_count <= 0 or ply_count % 2:
        raise ValueError("Symmetric stacks require a positive, even ply count")
    cycle = ["CF_TWILL", "GF"] if kind == "tw_glass" else ["CF_UD", "CF_TWILL", "GF"]
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
            "stack_kind": spec["kind"],
            "ply_sequence_bottom_to_top": "-".join(stack),
            "width_region": "in_range" if training_w[0] <= width <= training_w[1] else "extrapolation",
            "ply_count_region": "in_range" if training_n[0] <= nplies <= training_n[1] else "extrapolation",
        })
    return pd.DataFrame(a_rows), pd.DataFrame(target_rows)


def build_layup(cfg: dict, kind: str, ply_count: int,
                carbon_scale: float = 1.0, glass_scale: float = 1.0) -> dict:
    """Make a generic layup object accepted by the laminate property helper."""
    cards = copy.deepcopy(cfg["materials"])
    for name, card in cards.items():
        factor = carbon_scale if card["family"] == "carbon" else glass_scale
        for key in ELASTIC_KEYS:
            card[key] = float(card[key]) * float(factor)
    stack = symmetric_stack(kind, int(ply_count))
    return {"material_cards": cards, "plies": [
        {"material": name, "orientation_deg": 0.0} for name in stack
    ]}


def section_data(cfg: dict, row: dict, carbon_scale: float | None = None,
                 glass_scale: float | None = None) -> dict:
    """Calculate total thickness, mass and CLT section stiffnesses."""
    cscale = float(row.get("carbon_scale_truth", row.get("carbon_scale", 1.0)) if carbon_scale is None else carbon_scale)
    gscale = float(row.get("glass_scale_truth", row.get("glass_scale", 1.0)) if glass_scale is None else glass_scale)
    kind = str(row["stack_kind"])
    layup = build_layup(cfg, kind, int(row["ply_count"]), cscale, gscale)
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
    }


def _node_id(i: int, j: int, ny: int) -> int:
    return i * (ny + 1) + j + 1


def solve_shell_case(mapdl, cfg: dict, row: dict) -> dict:
    """Solve one generic symmetric layered shell cantilever in MAPDL."""
    cards = copy.deepcopy(cfg["materials"])
    cscale = float(row.get("carbon_scale_truth", row.get("carbon_scale", 1.0)))
    gscale = float(row.get("glass_scale_truth", row.get("glass_scale", 1.0)))
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
    props = section_data(cfg, row, cscale, gscale)
    return {
        **props, "tip_deflection_m": uz,
        "compliance_m_per_N": abs(uz) / load,
        "max_deflection_over_length": abs(uz) / length,
        "tip_load_N": load, "nx": nx, "ny": ny,
        "ply_count": len(stack),
        "ply_sequence_bottom_to_top": "-".join(stack),
    }


def run_fea(cfg: dict, design: pd.DataFrame, output: Path,
            dataset: str, max_cases: int | None = None) -> None:
    """Run/resume a dataset and persist each completed case immediately."""
    from ansys.mapdl.core import launch_mapdl
    output.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if output.exists():
        prior = pd.read_csv(output)
        if "case_id" in prior:
            done = {int(r["case_id"]): r for r in prior.to_dict(orient="records")}
    mapdl = launch_mapdl(nproc=int(cfg["mesh"].get("nproc", 2)))
    try:
        run_design = design if max_cases is None else design.head(max_cases)
        for row in run_design.to_dict(orient="records"):
            cid = int(row["case_id"])
            if cid in done:
                continue
            solved = solve_shell_case(mapdl, cfg, row)
            record = {**row, **solved}
            done[cid] = record
            pd.DataFrame(done.values()).sort_values("case_id").to_csv(output, index=False)
            print(f"{dataset} case {cid:03d}: width={row['width_m']:.4f} m, "
                  f"plies={row['ply_count']}, h={solved['total_thickness_m']*1e3:.3f} mm, "
                  f"C={solved['compliance_m_per_N']:.6e} m/N")
    finally:
        mapdl.exit()
    print(f"Saved {len(done)} {dataset} cases to {output}")


def mechanics_matrix(frame: pd.DataFrame, length: float) -> np.ndarray:
    return np.column_stack((
        length**3 / (3 * frame["EI_Nm2"].to_numpy(float)),
        length / frame["kGA_N"].to_numpy(float),
    ))


def fit_models(data_a: pd.DataFrame, cfg: dict) -> dict:
    """Fit all three model families using Dataset A only."""
    y = data_a["compliance_m_per_N"].to_numpy(float)
    if np.any(y <= 0):
        raise ValueError("Compliance must be positive")
    length = float(cfg["beam"]["length_m"])
    mech = LinearRegression(fit_intercept=False, positive=True).fit(mechanics_matrix(data_a, length), y)
    x = model_inputs(data_a)
    logy = np.log(y)
    quad = make_pipeline(StandardScaler(), PolynomialFeatures(degree=2, include_bias=False),
                         Ridge(alpha=float(cfg["models"].get("quadratic_ridge_alpha", 1e-8))))
    quad.fit(x, logy)
    kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(length_scale=np.ones(len(FEATURES)), nu=2.5) + WhiteKernel(noise_level=1e-5)
    gp = make_pipeline(StandardScaler(), GaussianProcessRegressor(
        kernel=kernel, normalize_y=True, n_restarts_optimizer=int(cfg["models"].get("gp_restarts", 2)),
        random_state=int(cfg["calibration"]["seed"])))
    gp.fit(x, logy)
    return {"mechanics": mech, "quadratic": quad, "gaussian_process": gp}


def cross_validated_training_metrics(data_a: pd.DataFrame, cfg: dict,
                                     models: dict) -> dict:
    """Five-fold CV predictions on Dataset A; no target data enters model selection."""
    folds = KFold(n_splits=5, shuffle=True, random_state=int(cfg["calibration"]["seed"]))
    y = data_a["compliance_m_per_N"].to_numpy(float)
    results = {}
    for name, prototype in models.items():
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
    return results


def predict(model_name: str, model, frame: pd.DataFrame, cfg: dict) -> np.ndarray:
    if model_name == "mechanics":
        return np.maximum(model.predict(mechanics_matrix(frame, float(cfg["beam"]["length_m"]))), 1e-15)
    return np.exp(np.clip(model.predict(model_inputs(frame)), -50, 50))


def metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    ape = 100 * np.abs(pred-y) / np.maximum(np.abs(y), 1e-15)
    return {
        "n": int(len(y)), "rmse_m_per_N": float(np.sqrt(mean_squared_error(y, pred))),
        "mae_m_per_N": float(mean_absolute_error(y, pred)),
        "mape_percent": float(np.mean(ape)), "max_ape_percent": float(np.max(ape)),
        "r2": float(r2_score(y, pred)) if len(y) > 1 else None,
    }


def generate_features_for_scales(cfg: dict, target_row: dict,
                                 cscale: float, gscale: float) -> dict:
    row = {**target_row, "carbon_scale_truth": cscale, "glass_scale_truth": gscale}
    props = section_data(cfg, row, cscale, gscale)
    return {**target_row, **props}


def fit_evaluate(cfg: dict, a_path: Path, target_path: Path, out_dir: Path) -> dict:
    """Fit on Dataset A and evaluate repeated two-case calibration per target set."""
    data_a = pd.read_csv(a_path)
    target_data = pd.read_csv(target_path)
    models = fit_models(data_a, cfg)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, model in models.items():
        joblib.dump(model, out_dir / f"model_{name}.joblib")
    train_metrics = {}
    for name, model in models.items():
        train_metrics[name] = metrics(data_a["compliance_m_per_N"].to_numpy(float),
                                      predict(name, model, data_a, cfg))
    cv_metrics = cross_validated_training_metrics(data_a, cfg, models)

    rng = np.random.default_rng(int(cfg["calibration"]["seed"]))
    n_repeats = int(cfg["calibration"]["repeats"])
    n_cal = int(cfg["calibration"]["points_per_repeat"])
    bounds = cfg["calibration"]["material_scale_bounds"]
    result_rows, metric_rows, split_rows = [], [], []
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
    summary = {
        "training_model_metrics": train_metrics,
        "training_fivefold_cv_metrics": cv_metrics,
        "experiment": {"dataset_a_cases": int(len(data_a)), "target_cases": int(len(target_data)),
                       "calibration_repeats": n_repeats, "calibration_points": n_cal},
        "calibration_warning": "Two points are the algebraic minimum for two factors. Inspect Jacobian rank and condition number before interpreting factors.",
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
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
