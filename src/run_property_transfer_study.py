"""Generate and evaluate a hidden-property, FEA-to-FEA calibration study.

Dataset A trains a mechanics-informed regression over equivalent beam
rigidities and two laminate property scales. Dataset B uses a fixed, withheld pair of property
scales at several known width/thickness combinations. A small, geometry-spread
subset of B calibrates the two scales; the remaining B cases are held out for
validation. Run with a licensed MAPDL installation for the FEA stages.

Example (run all 132 FEA cases, then fit/calibrate):
    python src/run_property_transfer_study.py --stage all

The separate stages are useful when the FEA datasets have already been run:
    --stage dataset-a
    --stage dataset-b
    --stage calibrate
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import joblib
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import qmc
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
    """Copy a layup with consistent thickness and orthotropic property scales.

    In-plane elastic terms (EX, EY, GXY) and transverse shear terms
    (GXZ, GYZ) are scaled separately. Poisson ratios and densities are fixed.
    The 20-ply sequence itself remains the supplied project sequence.
    """
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
    """Generate Dataset A by Latin hypercube over geometry and property scales."""
    spec = config["dataset_a"]
    n = int(spec["n_cases"])
    if n < 20:
        raise ValueError("Dataset A needs at least 20 cases for stable five-fold fitting")
    unit = qmc.LatinHypercube(d=4, seed=int(spec["seed"])).random(n)
    names = ["width_m", "thickness_scale", *PROPERTY_NAMES]
    bounds = [spec[name] for name in names]
    values = [lo + unit[:, i] * (hi - lo) for i, (lo, hi) in enumerate(bounds)]
    return pd.DataFrame({"case_id": np.arange(n), **{
        name: values[i] for i, name in enumerate(names)
    }})


def make_dataset_b_design(config: dict) -> pd.DataFrame:
    """Create a width/thickness grid at one fixed hidden property setting."""
    spec = config["dataset_b"]
    rows = []
    case_id = 0
    for thickness_scale in spec["thickness_scales"]:
        for width_m in spec["widths_m"]:
            rows.append({"case_id": case_id, "width_m": float(width_m),
                         "thickness_scale": float(thickness_scale)})
            case_id += 1
    return pd.DataFrame(rows)


def mechanics_basis(section_data: pd.DataFrame, length_m: float) -> np.ndarray:
    """Create bending/shear compliance bases from the equivalent rigidities.

    For a uniform cantilever, the Timoshenko terms scale as L^3/(3 EI) and
    L/(kGA). Their coefficients are learned from Dataset A FEA, rather than
    fixed to analytical values. This gives the regression the correct
    asymptotic behavior while allowing the FEA to set the correction factors.
    """
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
    """Choose a reproducible maximin subset across width and thickness."""
    if not 2 <= n_points < len(data):
        raise ValueError("Use at least two calibration cases and leave B holdouts")
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
    in_plane_values = np.broadcast_to(np.asarray(in_plane_scale, float), (len(geometry),))
    shear_values = np.broadcast_to(np.asarray(shear_scale, float), (len(geometry),))
    rows = []
    for position, row in enumerate(geometry.itertuples(index=False)):
        scaled = scale_material_layup(
            layup, float(row.thickness_scale), float(in_plane_values[position]),
            float(shear_values[position]))
        section = section_record(scaled, float(row.width_m), shear_correction)
        rows.append(section)
    return pd.DataFrame(rows, index=geometry.index)


def create_fea_data(config: dict, layup: dict, dataset: str,
                    output: Path) -> None:
    """Run MAPDL for Dataset A or B and save response plus section properties."""
    if dataset == "A":
        geometry = make_dataset_a_design(config)
        in_plane = geometry["in_plane_stiffness_scale"].to_numpy()
        shear = geometry["transverse_shear_stiffness_scale"].to_numpy()
    else:
        geometry = make_dataset_b_design(config)
        target = config["dataset_b"]
        in_plane = np.full(len(geometry), float(target["hidden_in_plane_stiffness_scale"]))
        shear = np.full(len(geometry), float(target["hidden_transverse_shear_stiffness_scale"]))

    from ansys.mapdl.core import launch_mapdl

    mapdl = launch_mapdl(nproc=int(config["mesh"].get("nproc", 2)))
    output_rows = []
    beam = config["beam"]
    try:
        for pos, row in enumerate(geometry.itertuples(index=False)):
            case_layup = scale_material_layup(
                layup, float(row.thickness_scale), float(in_plane[pos]), float(shear[pos]))
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
                int(config["mesh"].get("integration_points_per_ply", 3)))
            section = section_record(
                case_layup, float(row.width_m),
                float(beam.get("shear_correction_factor", 5.0 / 6.0)))
            result = {
                **row._asdict(),
                **section,
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
            "description": "Hidden target FEA material scales; excluded from calibration fitting.",
            "in_plane_stiffness_scale": float(config["dataset_b"]["hidden_in_plane_stiffness_scale"]),
            "transverse_shear_stiffness_scale": float(config["dataset_b"]["hidden_transverse_shear_stiffness_scale"]),
        }
        output.with_name("dataset_B_hidden_truth.json").write_text(json.dumps(truth, indent=2))


def fit_calibrate_validate(config: dict, layup: dict, data_a_path: Path,
                           data_b_path: Path, output_dir: Path,
                           n_calibration: int | None = None) -> dict:
    """Fit Model A, calibrate property scales on B subset, score B holdouts."""
    data_a = pd.read_csv(data_a_path)
    data_b = pd.read_csv(data_b_path)
    need_a = {"EI_Nm2", "kGA_N", "compliance_m_per_N"}
    need_b = set(OBSERVATION_COLUMNS)
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

    length_m = float(config["beam"]["length_m"])
    x_a = mechanics_basis(data_a, length_m)
    y_a = data_a["compliance_m_per_N"].to_numpy(float)
    cv = KFold(n_splits=5, shuffle=True, random_state=int(config["dataset_a"]["seed"]))
    cv_predictions = cross_val_predict(response_surface(), x_a, y_a, cv=cv)
    model = response_surface().fit(x_a, y_a)

    n_cal = int(n_calibration if n_calibration is not None
                else config["calibration"]["n_points"])
    cal_idx = select_calibration_indices(
        data_b, n_cal, int(config["calibration"].get("random_seed", 0)))
    test_mask = np.ones(len(data_b), dtype=bool)
    test_mask[cal_idx] = False
    cal = data_b.iloc[cal_idx].reset_index(drop=True)
    test = data_b.loc[test_mask].reset_index(drop=True)
    beam = config["beam"]
    kappa = float(beam.get("shear_correction_factor", 5.0 / 6.0))

    def predict_at_scales(geometry: pd.DataFrame, scales: np.ndarray) -> np.ndarray:
        section = design_properties(layup, geometry, float(scales[0]),
                                    float(scales[1]), kappa)
        return model.predict(mechanics_basis(section, length_m))

    nominal_scales = np.asarray(config["calibration"].get("initial_scales", [1.0, 1.0]), float)
    cal_observed = cal["compliance_m_per_N"].to_numpy(float)

    def relative_residual(scales: np.ndarray) -> np.ndarray:
        predicted = predict_at_scales(cal, scales)
        return (predicted - cal_observed) / cal_observed

    lo = [float(config["calibration"]["in_plane_scale_bounds"][0]),
          float(config["calibration"]["transverse_shear_scale_bounds"][0])]
    hi = [float(config["calibration"]["in_plane_scale_bounds"][1]),
          float(config["calibration"]["transverse_shear_scale_bounds"][1])]
    if np.any(nominal_scales <= lo) or np.any(nominal_scales >= hi):
        raise ValueError("Initial property scales must lie strictly inside calibration bounds")
    fit = least_squares(relative_residual, nominal_scales, bounds=(lo, hi), method="trf")
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"Calibration did not converge: {fit.message}")
    fitted_scales = fit.x

    nominal_prediction = predict_at_scales(test, nominal_scales)
    calibrated_prediction = predict_at_scales(test, fitted_scales)
    test_result = test.copy()
    test_result["nominal_prediction_m_per_N"] = nominal_prediction
    test_result["calibrated_prediction_m_per_N"] = calibrated_prediction
    test_result["nominal_error_percent"] = 100 * (nominal_prediction - test_result["compliance_m_per_N"]) / test_result["compliance_m_per_N"]
    test_result["calibrated_error_percent"] = 100 * (calibrated_prediction - test_result["compliance_m_per_N"]) / test_result["compliance_m_per_N"]
    test_result["is_calibration_case"] = False

    cal_result = cal.copy()
    cal_result["nominal_prediction_m_per_N"] = predict_at_scales(cal, nominal_scales)
    cal_result["calibrated_prediction_m_per_N"] = predict_at_scales(cal, fitted_scales)
    cal_result["nominal_error_percent"] = 100 * (cal_result["nominal_prediction_m_per_N"] - cal_result["compliance_m_per_N"]) / cal_result["compliance_m_per_N"]
    cal_result["calibrated_error_percent"] = 100 * (cal_result["calibrated_prediction_m_per_N"] - cal_result["compliance_m_per_N"]) / cal_result["compliance_m_per_N"]
    cal_result["is_calibration_case"] = True
    predictions = pd.concat([cal_result, test_result], ignore_index=True).sort_values("case_id")

    truth_path = data_b_path.with_name("dataset_B_hidden_truth.json")
    hidden_truth = json.loads(truth_path.read_text()) if truth_path.exists() else None
    metrics = {
        "method": "Dataset A mechanics-informed regression; Dataset B least-squares calibration of in-plane and transverse-shear stiffness scales",
        "data_a": str(data_a_path), "data_b": str(data_b_path),
        "n_dataset_a": len(data_a), "n_dataset_b": len(data_b),
        "n_calibration": len(cal), "calibration_case_ids": cal["case_id"].astype(int).tolist(),
        "validation_case_ids": test["case_id"].astype(int).tolist(),
        "model_features": MODEL_FEATURES,
        "surrogate_dataset_a_5fold_cv": scores(y_a, cv_predictions),
        "nominal_scales": dict(zip(PROPERTY_NAMES, nominal_scales.tolist())),
        "fitted_scales": dict(zip(PROPERTY_NAMES, fitted_scales.tolist())),
        "hidden_target_scales_for_post_fit_audit": hidden_truth,
        "calibration_relative_residual_sse": float(np.sum(fit.fun**2)),
        "nominal_dataset_b_holdout": scores(test["compliance_m_per_N"].to_numpy(), nominal_prediction),
        "calibrated_dataset_b_holdout": scores(test["compliance_m_per_N"].to_numpy(), calibrated_prediction),
        "interpretation": "B calibration cases are used only to fit the two material scales. B holdouts are excluded from surrogate fitting and scale fitting. Both datasets are FEA; this tests transfer between FEA property settings, not physical accuracy.",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_dir / "dataset_b_predictions.csv", index=False)
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    joblib.dump({"model_a": model, "features": MODEL_FEATURES,
                 "fitted_scales": dict(zip(PROPERTY_NAMES, fitted_scales.tolist()))},
                output_dir / "model_a_calibrated.joblib")
    print(json.dumps(metrics, indent=2))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("all", "dataset-a", "dataset-b", "calibrate"), default="all")
    parser.add_argument("--config", default="config.hidden_property.example.json")
    parser.add_argument("--layup", default="layup.json")
    parser.add_argument("--dataset-a", default="results/property_transfer/dataset_A.csv")
    parser.add_argument("--dataset-b", default="results/property_transfer/dataset_B.csv")
    parser.add_argument("--output-dir", default="results/property_transfer/model")
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
