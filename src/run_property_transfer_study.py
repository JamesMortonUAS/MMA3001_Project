"""Generate a tensor-grid, hidden-property FEA transfer study.

Dataset A is a complete linearly spaced grid in width, thickness, in-plane
stiffness and transverse-shear stiffness. Configured laminate property pairs
are excluded from training. Dataset B evaluates one excluded property pair
over in-range and out-of-range geometry. One observation fits one common
compliance correction factor; it cannot identify independent in-plane and
shear stiffness scales. Run the FEA stages with licensed MAPDL.

Example (72 training FEA cases plus 25 target-laminate FEA cases):
    python src/run_property_transfer_study.py --stage all

Stages:
    --stage dataset-a
    --stage dataset-b
    --stage calibrate
"""
from __future__ import annotations

import argparse
import copy
import itertools
import json
from pathlib import Path
from types import SimpleNamespace

import joblib
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict

from laminate_section import effective_bending_shear_stiffness
from run_laminate_shell_sweep import solve_case


MODEL_FEATURES = ["bending_compliance_basis_m_per_N", "shear_compliance_basis_m_per_N"]
OBSERVATION_COLUMNS = ["case_id", "width_m", "thickness_scale", "compliance_m_per_N"]
PROPERTY_NAMES = ["in_plane_stiffness_scale", "transverse_shear_stiffness_scale"]


def scale_material_layup(layup: dict, thickness_scale: float,
                         in_plane_scale: float,
                         transverse_shear_scale: float) -> dict:
    """Copy a layup with separate in-plane and transverse-shear property scales."""
    values = (thickness_scale, in_plane_scale, transverse_shear_scale)
    if any(not np.isfinite(value) or value <= 0 for value in values):
        raise ValueError("Thickness and stiffness scales must be positive and finite")
    out = copy.deepcopy(layup)
    for card in out["material_cards"].values():
        card["ply_thickness_m"] = float(card["ply_thickness_m"]) * thickness_scale
        for key in ("EX_Pa", "EY_Pa", "GXY_Pa"):
            card[key] = float(card[key]) * in_plane_scale
        for key in ("GXZ_Pa", "GYZ_Pa"):
            card[key] = float(card[key]) * transverse_shear_scale
    return out


def section_record(layup: dict, width_m: float, shear_correction: float) -> dict:
    """Return CLT laminate and beam-equivalent bending/shear properties."""
    if width_m <= 0 or shear_correction <= 0:
        raise ValueError("Width and shear correction factor must be positive")
    d_eff, s_eff = effective_bending_shear_stiffness(layup)
    total_thickness = sum(
        float(layup["material_cards"][record["material"]]["ply_thickness_m"])
        for record in layup["plies"]
    )
    areal_mass = sum(
        float(layup["material_cards"][record["material"]]["density_kg_m3"])
        * float(layup["material_cards"][record["material"]]["ply_thickness_m"])
        for record in layup["plies"]
    )
    return {
        "total_thickness_m": total_thickness,
        "D_eff_Nm_per_m": d_eff,
        "S_eff_N_per_m": s_eff,
        "EI_Nm2": width_m * d_eff,
        "GA_N": width_m * s_eff,
        "kGA_N": shear_correction * width_m * s_eff,
        "areal_mass_kg_m2": areal_mass,
    }


def make_dataset_a_design(config: dict) -> pd.DataFrame:
    """Create a full tensor-product grid with linear spacing on each axis."""
    spec = config["dataset_a"]
    levels = int(spec["grid_levels"])
    if levels < 2:
        raise ValueError("dataset_a.grid_levels must be at least 2")
    names = ["width_m", "thickness_scale", *PROPERTY_NAMES]
    axes = []
    for name in names:
        lower, upper = map(float, spec[name])
        if not np.isfinite(lower + upper) or lower <= 0 or upper <= lower:
            raise ValueError(f"Invalid positive bounds for {name}: {spec[name]}")
        axes.append(np.linspace(lower, upper, levels))
    rows = [dict(zip(names, point)) for point in itertools.product(*axes)]
    design = pd.DataFrame(rows)
    design.insert(0, "case_id", np.arange(len(design), dtype=int))
    return design


def filter_training_laminates(data: pd.DataFrame, excluded_pairs: list[dict]) -> pd.DataFrame:
    """Remove configured in-plane/shear scale pairs from Dataset A training."""
    keep = np.ones(len(data), dtype=bool)
    for pair in excluded_pairs:
        ip = float(pair["in_plane_stiffness_scale"])
        shear = float(pair["transverse_shear_stiffness_scale"])
        matches = (
            np.isclose(data["in_plane_stiffness_scale"].to_numpy(float), ip, rtol=0, atol=1e-10)
            & np.isclose(data["transverse_shear_stiffness_scale"].to_numpy(float),
                         shear, rtol=0, atol=1e-10)
        )
        keep &= ~matches
    return data.loc[keep].reset_index(drop=True)


def make_dataset_b_design(config: dict) -> pd.DataFrame:
    """Create target-laminate geometry grid and tag in-range/extrapolation cases."""
    spec = config["dataset_b"]
    a_spec = config["dataset_a"]
    w_lo, w_hi = map(float, a_spec["width_m"])
    t_lo, t_hi = map(float, a_spec["thickness_scale"])
    rows = []
    case_id = 0
    for thickness_scale in spec["thickness_scales"]:
        for width_m in spec["widths_m"]:
            inside = w_lo <= float(width_m) <= w_hi and t_lo <= float(thickness_scale) <= t_hi
            rows.append({
                "case_id": case_id,
                "width_m": float(width_m),
                "thickness_scale": float(thickness_scale),
                "geometry_region": "interpolation" if inside else "extrapolation",
            })
            case_id += 1
    return pd.DataFrame(rows)


def mechanics_basis(section_data: pd.DataFrame, length_m: float) -> np.ndarray:
    """Create bending/shear compliance bases from equivalent section rigidities."""
    ei = section_data["EI_Nm2"].to_numpy(float)
    kga = section_data["kGA_N"].to_numpy(float)
    if length_m <= 0 or np.any(ei <= 0) or np.any(kga <= 0):
        raise ValueError("Length and equivalent rigidities must be positive")
    return np.column_stack((length_m**3 / (3.0 * ei), length_m / kga))


def response_surface():
    """Fit positive bending and shear coefficients to Dataset A FEA."""
    return LinearRegression(fit_intercept=False, positive=True)


def select_calibration_indices(data: pd.DataFrame, n_points: int,
                               seed: int) -> np.ndarray:
    """Select the centre case first, then maximin cases if requested."""
    if not 1 <= n_points < len(data):
        raise ValueError("Calibration must use at least one case and leave validation cases")
    geo = data[["width_m", "thickness_scale"]].to_numpy(float)
    mean = geo.mean(axis=0)
    scale = geo.std(axis=0)
    scale[scale == 0] = 1.0
    z = (geo - mean) / scale
    rng = np.random.default_rng(seed)
    centre_d = np.linalg.norm(z - z.mean(axis=0), axis=1)
    candidates = np.flatnonzero(np.isclose(centre_d, centre_d.min()))
    chosen = [int(rng.choice(candidates))]
    while len(chosen) < n_points:
        remaining = np.asarray([i for i in range(len(data)) if i not in chosen])
        min_d = np.min(np.linalg.norm(z[remaining, None, :] - z[chosen][None, :, :], axis=2), axis=1)
        farthest = remaining[np.isclose(min_d, min_d.max())]
        chosen.append(int(rng.choice(farthest)))
    return np.asarray(chosen, dtype=int)


def scores(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    y_true, y_pred = np.asarray(y_true, float), np.asarray(y_pred, float)
    ape = 100.0 * np.abs(y_pred - y_true) / np.maximum(np.abs(y_true), 1e-15)
    return {
        "n": int(y_true.size),
        "rmse_m_per_N": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "mae_m_per_N": float(mean_absolute_error(y_true, y_pred)),
        "mape_percent": float(np.mean(ape)),
        "max_ape_percent": float(np.max(ape)),
        "r2": float(r2_score(y_true, y_pred)) if y_true.size > 1 else None,
    }


def design_properties(layup: dict, geometry: pd.DataFrame,
                      in_plane_scale: float, shear_scale: float,
                      shear_correction: float) -> pd.DataFrame:
    rows = []
    for row in geometry.itertuples(index=False):
        scaled = scale_material_layup(
            layup, float(row.thickness_scale), in_plane_scale, shear_scale
        )
        rows.append(section_record(scaled, float(row.width_m), shear_correction))
    return pd.DataFrame(rows, index=geometry.index)


def create_fea_data(config: dict, layup: dict, dataset: str,
                    output: Path) -> None:
    """Run MAPDL for training Dataset A or target-laminate Dataset B."""
    if dataset == "A":
        full_design = make_dataset_a_design(config)
        geometry = filter_training_laminates(
            full_design, config["dataset_a"].get("excluded_laminates", [])
        )
        excluded_count = len(full_design) - len(geometry)
        print(f"Dataset A: {len(geometry)} training cases; "
              f"{excluded_count} cases from excluded laminate pairs not run.")
        in_plane = geometry["in_plane_stiffness_scale"].to_numpy(float)
        shear = geometry["transverse_shear_stiffness_scale"].to_numpy(float)
    elif dataset == "B":
        geometry = make_dataset_b_design(config)
        target = config["dataset_b"]
        in_plane = np.full(len(geometry), float(target["hidden_in_plane_stiffness_scale"]))
        shear = np.full(len(geometry), float(target["hidden_transverse_shear_stiffness_scale"]))
    else:
        raise ValueError("dataset must be 'A' or 'B'")

    from ansys.mapdl.core import launch_mapdl

    mapdl = launch_mapdl(nproc=int(config["mesh"].get("nproc", 2)))
    output_rows = []
    beam = config["beam"]
    try:
        for pos, row in enumerate(geometry.itertuples(index=False)):
            case_layup = scale_material_layup(
                layup, float(row.thickness_scale), float(in_plane[pos]), float(shear[pos])
            )
            design_row = SimpleNamespace(
                length_m=float(beam["length_m"]),
                root_width_m=float(row.width_m),
                taper_ratio=1.0,
                tip_load_N=float(beam["tip_load_N"]),
                thickness_scale=1.0,
                effective_stiffness_scale=1.0,
            )
            uz, mass, total_thickness = solve_case(
                mapdl, design_row, case_layup,
                int(config["mesh"]["nx"]), int(config["mesh"]["ny"]),
                int(config["mesh"].get("integration_points_per_ply", 3))
            )
            section = section_record(
                case_layup, float(row.width_m),
                float(beam.get("shear_correction_factor", 5.0 / 6.0))
            )
            result = {
                **row._asdict(), **section,
                "laminate_thickness_from_fea_m": total_thickness,
                "mass_kg": mass,
                "tip_load_N": float(beam["tip_load_N"]),
                "tip_deflection_m": uz,
                "compliance_m_per_N": abs(uz) / float(beam["tip_load_N"]),
                "max_deflection_over_length": abs(uz) / float(beam["length_m"]),
                "nx": int(config["mesh"]["nx"]),
                "ny": int(config["mesh"]["ny"]),
            }
            if dataset == "A":
                result["in_plane_stiffness_scale"] = float(in_plane[pos])
                result["transverse_shear_stiffness_scale"] = float(shear[pos])
            output_rows.append(result)
            print(f"Dataset {dataset} case {int(row.case_id):03d}: "
                  f"width={row.width_m:.4f} m, t-scale={row.thickness_scale:.3f}, "
                  f"C={result['compliance_m_per_N']:.6e} m/N")
    finally:
        mapdl.exit()

    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output_rows).to_csv(output, index=False)
    print(f"Saved Dataset {dataset}: {output}")
    if dataset == "B":
        truth = {
            "description": "Hidden Dataset B material scales; post-fit audit only.",
            "in_plane_stiffness_scale": float(config["dataset_b"]["hidden_in_plane_stiffness_scale"]),
            "transverse_shear_stiffness_scale": float(config["dataset_b"]["hidden_transverse_shear_stiffness_scale"]),
        }
        output.with_name("dataset_B_hidden_truth.json").write_text(json.dumps(truth, indent=2))


def fit_calibrate_validate(config: dict, layup: dict, data_a_path: Path,
                           data_b_path: Path, output_dir: Path,
                           n_calibration: int | None = None) -> dict:
    """Fit Model A, calibrate one common correction from one B case, and validate."""
    data_a = pd.read_csv(data_a_path)
    data_b = pd.read_csv(data_b_path)
    need_a = {"EI_Nm2", "kGA_N", "compliance_m_per_N"}
    need_b = set(OBSERVATION_COLUMNS) | {"geometry_region"}
    if missing := sorted(need_a - set(data_a.columns)):
        raise ValueError(f"Dataset A missing columns: {missing}")
    if missing := sorted(need_b - set(data_b.columns)):
        raise ValueError(f"Dataset B missing columns: {missing}")
    if len(data_a) < 20 or len(data_b) < 3:
        raise ValueError("Need at least 20 Dataset A rows and 3 Dataset B rows")
    if (data_a[["EI_Nm2", "kGA_N", "compliance_m_per_N"]] <= 0).any().any():
        raise ValueError("Dataset A rigidities and compliances must be positive")
    if (data_b[["width_m", "thickness_scale", "compliance_m_per_N"]] <= 0).any().any():
        raise ValueError("Dataset B geometry and compliances must be positive")
    if not set(data_b["geometry_region"].unique()) <= {"interpolation", "extrapolation"}:
        raise ValueError("Dataset B geometry_region must be interpolation or extrapolation")

    length_m = float(config["beam"]["length_m"])
    x_a = mechanics_basis(data_a, length_m)
    y_a = data_a["compliance_m_per_N"].to_numpy(float)
    cv = KFold(n_splits=5, shuffle=True, random_state=int(config["dataset_a"].get("seed", 0)))
    cv_predictions = cross_val_predict(response_surface(), x_a, y_a, cv=cv)
    model = response_surface().fit(x_a, y_a)

    n_cal = int(n_calibration if n_calibration is not None
                else config["calibration"].get("n_points", 1))
    cal_idx = select_calibration_indices(
        data_b, n_cal, int(config["calibration"].get("random_seed", 0))
    )
    test_mask = np.ones(len(data_b), dtype=bool)
    test_mask[cal_idx] = False
    cal = data_b.iloc[cal_idx].reset_index(drop=True)
    validation = data_b.loc[test_mask].reset_index(drop=True)
    beam = config["beam"]
    kappa = float(beam.get("shear_correction_factor", 5.0 / 6.0))

    def predict_at_common_scale(geometry: pd.DataFrame, common_scale: float) -> np.ndarray:
        section = design_properties(layup, geometry, common_scale, common_scale, kappa)
        return model.predict(mechanics_basis(section, length_m))

    initial = float(config["calibration"].get("initial_common_stiffness_scale", 1.0))
    lo, hi = map(float, config["calibration"]["common_stiffness_scale_bounds"])
    if not lo < initial < hi:
        raise ValueError("Initial common stiffness scale must lie strictly inside its bounds")
    observed = cal["compliance_m_per_N"].to_numpy(float)

    def relative_residual(scales: np.ndarray) -> np.ndarray:
        predicted = predict_at_common_scale(cal, float(scales[0]))
        return (predicted - observed) / observed

    fit = least_squares(relative_residual, x0=np.asarray([initial]),
                        bounds=(np.asarray([lo]), np.asarray([hi])), method="trf")
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"Calibration did not converge: {fit.message}")
    fitted_scale = float(fit.x[0])

    nominal_prediction = predict_at_common_scale(validation, initial)
    calibrated_prediction = predict_at_common_scale(validation, fitted_scale)
    validation_result = validation.copy()
    validation_result["nominal_prediction_m_per_N"] = nominal_prediction
    validation_result["calibrated_prediction_m_per_N"] = calibrated_prediction
    validation_result["nominal_error_percent"] = (
        100 * (nominal_prediction - validation_result["compliance_m_per_N"])
        / validation_result["compliance_m_per_N"]
    )
    validation_result["calibrated_error_percent"] = (
        100 * (calibrated_prediction - validation_result["compliance_m_per_N"])
        / validation_result["compliance_m_per_N"]
    )
    validation_result["is_calibration_case"] = False

    cal_result = cal.copy()
    cal_result["nominal_prediction_m_per_N"] = predict_at_common_scale(cal, initial)
    cal_result["calibrated_prediction_m_per_N"] = predict_at_common_scale(cal, fitted_scale)
    cal_result["nominal_error_percent"] = (
        100 * (cal_result["nominal_prediction_m_per_N"] - cal_result["compliance_m_per_N"])
        / cal_result["compliance_m_per_N"]
    )
    cal_result["calibrated_error_percent"] = (
        100 * (cal_result["calibrated_prediction_m_per_N"] - cal_result["compliance_m_per_N"])
        / cal_result["compliance_m_per_N"]
    )
    cal_result["is_calibration_case"] = True
    predictions = pd.concat([cal_result, validation_result], ignore_index=True).sort_values("case_id")

    group_scores = {}
    for region in ("interpolation", "extrapolation"):
        group = validation_result.loc[validation_result["geometry_region"] == region]
        if not group.empty:
            group_scores[region] = {
                "nominal": scores(group["compliance_m_per_N"].to_numpy(),
                                  group["nominal_prediction_m_per_N"].to_numpy()),
                "calibrated": scores(group["compliance_m_per_N"].to_numpy(),
                                     group["calibrated_prediction_m_per_N"].to_numpy()),
            }

    truth_path = data_b_path.with_name("dataset_B_hidden_truth.json")
    hidden_truth = json.loads(truth_path.read_text()) if truth_path.exists() else None
    excluded_pairs = config["dataset_a"].get("excluded_laminates", [])
    metrics = {
        "method": "Dataset A tensor-grid mechanics regression; one-point calibration of a common stiffness scale",
        "data_a": str(data_a_path), "data_b": str(data_b_path),
        "n_dataset_a_training": len(data_a), "n_dataset_b": len(data_b),
        "n_calibration": len(cal), "calibration_case_ids": cal["case_id"].astype(int).tolist(),
        "validation_case_ids": validation["case_id"].astype(int).tolist(),
        "validation_regions": {
            "interpolation_case_ids": validation.loc[
                validation["geometry_region"] == "interpolation", "case_id"
            ].astype(int).tolist(),
            "extrapolation_case_ids": validation.loc[
                validation["geometry_region"] == "extrapolation", "case_id"
            ].astype(int).tolist(),
        },
        "dataset_a_excluded_laminate_pairs": excluded_pairs,
        "model_features": MODEL_FEATURES,
        "surrogate_dataset_a_5fold_cv": scores(y_a, cv_predictions),
        "nominal_common_stiffness_scale": initial,
        "fitted_common_stiffness_scale": fitted_scale,
        "hidden_target_scales_for_post_fit_audit": hidden_truth,
        "calibration_relative_residual_sse": float(np.sum(fit.fun**2)),
        "nominal_validation": scores(validation["compliance_m_per_N"].to_numpy(),
                                     nominal_prediction),
        "calibrated_validation": scores(validation["compliance_m_per_N"].to_numpy(),
                                        calibrated_prediction),
        "validation_by_region": group_scores,
        "interpretation": (
            "One compliance observation fits one common stiffness correction, not "
            "independent in-plane and transverse-shear scales. Interpolation and "
            "extrapolation are classified by whether both geometry inputs lie inside "
            "the Dataset A width/thickness bounds. FEA-to-FEA only."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_dir / "dataset_b_predictions.csv", index=False)
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    joblib.dump({
        "model_a": model,
        "features": MODEL_FEATURES,
        "fitted_common_stiffness_scale": fitted_scale,
    }, output_dir / "model_a_calibrated.joblib")
    print(json.dumps(metrics, indent=2))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("all", "dataset-a", "dataset-b", "calibrate"), default="all")
    parser.add_argument("--config", default="config.hidden_property.example.json")
    parser.add_argument("--layup", default="layup.json")
    parser.add_argument("--dataset-a", default="results/property_transfer_grid/dataset_A.csv")
    parser.add_argument("--dataset-b", default="results/property_transfer_grid/dataset_B.csv")
    parser.add_argument("--output-dir", default="results/property_transfer_grid/model")
    parser.add_argument("--calibration-points", type=int, default=None)
    args = parser.parse_args()

    config = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text())
    data_a_path, data_b_path = Path(args.dataset_a), Path(args.dataset_b)
    if args.stage in ("all", "dataset-a"):
        create_fea_data(config, layup, "A", data_a_path)
    if args.stage in ("all", "dataset-b"):
        create_fea_data(config, layup, "B", data_b_path)
    if args.stage in ("all", "calibrate"):
        fit_calibrate_validate(config, layup, data_a_path, data_b_path,
                               Path(args.output_dir), args.calibration_points)


if __name__ == "__main__":
    main()
