"""Fit a quadratic FEA response surface and calibrate one stiffness factor."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from scipy.optimize import least_squares

from analytic import tip_compliance_from_stiffness
from laminate_section import effective_bending_shear_stiffness, scaled_layup

FEATURES = ["length_m", "root_width_m", "taper_ratio", "thickness_scale",
            "effective_stiffness_scale"]
GEOMETRY_FEATURES = ["length_m", "root_width_m", "taper_ratio"]


def score(y_true, y_pred) -> dict:
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    rel = np.abs(y_pred - y_true) / np.maximum(np.abs(y_true), 1e-15)
    return {
        "rmse_m_per_N": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "mae_m_per_N": float(mean_absolute_error(y_true, y_pred)),
        "r2": float(r2_score(y_true, y_pred)),
        "mean_absolute_percentage_error_percent": float(100 * np.mean(rel)),
        "max_absolute_percentage_error_percent": float(100 * np.max(rel)),
    }


def quadratic_model() -> object:
    """Full second-order polynomial surface with a numerically stable fit."""
    return make_pipeline(StandardScaler(),
                         PolynomialFeatures(degree=2, include_bias=False),
                         Ridge(alpha=1e-10))


def analytical_reference(data: pd.DataFrame, config: dict, layup: dict,
                         thickness_scale: float,
                         stiffness_scale: float) -> np.ndarray:
    """Calculate reference compliance at a common section-property setting."""
    kappa = float(config["fea"].get("shear_correction_factor_analytic", 5 / 6))
    return np.asarray([
        tip_compliance_from_stiffness(
            row.length_m, row.root_width_m, row.taper_ratio,
            *effective_bending_shear_stiffness(
                scaled_layup(layup, thickness_scale, stiffness_scale)), kappa)
        for row in data.itertuples(index=False)
    ], dtype=float)


def select_calibration_indices(data: pd.DataFrame, train_idx: np.ndarray,
                               n: int) -> np.ndarray:
    """Select training geometries near the centre, then spread them out."""
    geometry = data[GEOMETRY_FEATURES].to_numpy(float)
    centre_scale = StandardScaler().fit(geometry[train_idx])
    z = centre_scale.transform(geometry)
    centre = z[train_idx].mean(axis=0)
    order = train_idx[np.argsort(np.linalg.norm(z[train_idx] - centre, axis=1))]
    if n == 1:
        return order[:1]
    selected = [int(order[0])]
    while len(selected) < n:
        remaining = [int(i) for i in order if int(i) not in selected]
        distances = [min(np.linalg.norm(z[i] - z[j]) for j in selected) for i in remaining]
        selected.append(remaining[int(np.argmax(distances))])
    return np.asarray(selected, dtype=int)


def fit_effective_stiffness_scale(model, x_cal: np.ndarray,
                                  reference_compliance: np.ndarray,
                                  initial_scale: float,
                                  scale_bounds: tuple[float, float]) -> tuple[float, float]:
    """Fit the shared modulus scale by bounded least squares.

    ``x_cal`` contains calibration geometries and thicknesses; its effective
    stiffness-scale column is replaced by the candidate scalar during fitting.
    Returns the fitted scale and the sum of squared residuals.
    """
    x_cal = np.asarray(x_cal, dtype=float)
    reference_compliance = np.asarray(reference_compliance, dtype=float)
    lo, hi = map(float, scale_bounds)
    scale_col = FEATURES.index("effective_stiffness_scale")
    if x_cal.ndim != 2 or x_cal.shape[1] != len(FEATURES):
        raise ValueError("x_cal must contain one column for each model feature")
    if len(reference_compliance) != len(x_cal) or len(x_cal) < 1:
        raise ValueError("Provide one positive reference compliance per calibration case")
    if np.any(reference_compliance <= 0) or not lo < initial_scale < hi:
        raise ValueError("Reference compliances must be positive and initial scale within bounds")

    def residual(parameter):
        candidate = x_cal.copy()
        candidate[:, scale_col] = parameter[0]
        return model.predict(candidate) - reference_compliance

    result = least_squares(residual, x0=np.asarray([initial_scale]),
                           bounds=([lo], [hi]), method="trf")
    if not result.success or not np.isfinite(result.x[0]):
        raise RuntimeError(f"Effective stiffness calibration failed: {result.message}")
    return float(result.x[0]), float(2 * result.cost)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="CSV from run_laminate_shell_sweep.py")
    ap.add_argument("--config", default="config.quadratic.example.json")
    ap.add_argument("--layup", default="layup.json")
    ap.add_argument("--output-dir", default="results/quadratic_model")
    args = ap.parse_args()

    data = pd.read_csv(args.input)
    config = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text())
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    required = FEATURES + ["D_eff_Nm", "S_eff_N_per_m", "compliance_m_per_N"]
    missing = sorted(set(required) - set(data.columns))
    if missing:
        raise ValueError(f"Input CSV missing required columns: {missing}")
    data = data.dropna(subset=required).reset_index(drop=True)
    x = data[FEATURES].to_numpy(float)
    y_fea = data["compliance_m_per_N"].to_numpy(float)
    if len(data) < 42 or np.any(y_fea <= 0):
        raise ValueError("Need at least 42 complete cases and positive compliances")

    mcfg = config.get("model", {})
    seed = int(mcfg.get("seed", 1234))
    train_idx, test_idx = train_test_split(
        np.arange(len(data)), test_size=float(mcfg.get("test_fraction", 0.2)),
        random_state=seed)
    if len(train_idx) < 21:
        raise ValueError("Training set must exceed the 21 coefficients in a 5-input quadratic")
    model = quadratic_model().fit(x[train_idx], y_fea[train_idx])
    y_quad_at_sweep_inputs = model.predict(x[test_idx])

    # Fit the shared effective stiffness scale directly. Calibration and
    # geometry-transfer checks hold thickness fixed to isolate geometry.
    ncal = int(mcfg.get("calibration_points", 1))
    if ncal < 1 or ncal > len(train_idx):
        raise ValueError("calibration_points must be between 1 and n_train")
    cal_idx = select_calibration_indices(data, train_idx, ncal)
    nominal_scale = float(mcfg.get("initial_effective_stiffness_scale", 1.0))
    truth_scale = float(mcfg.get("reference_effective_stiffness_scale", 1.0))
    calibration_thickness = float(mcfg.get("calibration_thickness_scale", 1.0))
    scale_bounds = tuple(map(float, config["sweep"]["effective_stiffness_scale"]))
    if not scale_bounds[0] < scale_bounds[1]:
        raise ValueError("effective_stiffness_scale sweep bounds must increase")
    if not scale_bounds[0] < nominal_scale < scale_bounds[1]:
        raise ValueError("initial_effective_stiffness_scale must lie inside the sweep bounds")
    if not scale_bounds[0] <= truth_scale <= scale_bounds[1]:
        raise ValueError("reference_effective_stiffness_scale must lie inside the sweep bounds")
    if not (config["sweep"]["thickness_scale"][0] <= calibration_thickness
            <= config["sweep"]["thickness_scale"][1]):
        raise ValueError("calibration_thickness_scale must lie inside the thickness sweep")

    y_reference_cal = analytical_reference(
        data.iloc[cal_idx], config, layup, calibration_thickness, truth_scale)
    x_cal = x[cal_idx].copy()
    x_cal[:, FEATURES.index("thickness_scale")] = calibration_thickness

    calibrated_scale, calibration_sse = fit_effective_stiffness_scale(
        model, x_cal, y_reference_cal, nominal_scale, scale_bounds)

    x_test_fixed_section = x[test_idx].copy()
    x_test_fixed_section[:, FEATURES.index("thickness_scale")] = calibration_thickness
    x_test_nominal = x_test_fixed_section.copy()
    x_test_nominal[:, FEATURES.index("effective_stiffness_scale")] = nominal_scale
    x_test_calibrated = x_test_fixed_section.copy()
    x_test_calibrated[:, FEATURES.index("effective_stiffness_scale")] = calibrated_scale
    y_nominal_fixed_section = model.predict(x_test_nominal)
    y_calibrated_fixed_section = model.predict(x_test_calibrated)
    test_geometry = data.iloc[test_idx].copy()
    y_reference_test = analytical_reference(
        test_geometry, config, layup, calibration_thickness, truth_scale)

    metrics = {
        "input_file": str(args.input),
        "data_origin": str(data["data_origin"].iloc[0]) if "data_origin" in data else "PyMAPDL SHELL181 FEA",
        "n_total": int(len(data)), "n_train": int(len(train_idx)), "n_test": int(len(test_idx)),
        "n_analytical_calibration": int(ncal),
        "calibration_case_ids": [int(data.iloc[i]["case_id"]) for i in cal_idx],
        "geometry_transfer_test_case_ids": [int(data.iloc[i]["case_id"]) for i in test_idx],
        "features": FEATURES, "target": "tip compliance |delta|/P [m/N]",
        "response_surface": "full quadratic polynomial (21 coefficients for five inputs)",
        "calibration": {
            "method": "bounded nonlinear least-squares fit of one effective stiffness scale",
            "initial_effective_stiffness_scale": nominal_scale,
            "calibrated_effective_stiffness_scale": calibrated_scale,
            "reference_effective_stiffness_scale": truth_scale,
            "calibration_thickness_scale": calibration_thickness,
            "cost_sum_squared_residuals": calibration_sse,
            "relationship": "select s_E minimizing sum((C_analytical - C_quadratic(geometry, thickness, s_E))^2)",
        },
        "quadratic_vs_fea_test_at_sweep_inputs": score(y_fea[test_idx], y_quad_at_sweep_inputs),
        "nominal_quadratic_vs_analytical_geometry_test": score(y_reference_test, y_nominal_fixed_section),
        "calibrated_quadratic_vs_analytical_geometry_test": score(y_reference_test, y_calibrated_fixed_section),
    }
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    pred = data.iloc[test_idx].copy()
    pred["quadratic_fea_compliance_m_per_N"] = y_quad_at_sweep_inputs
    pred["geometry_test_thickness_scale"] = calibration_thickness
    pred["reference_effective_stiffness_scale"] = truth_scale
    pred["initial_effective_stiffness_scale"] = nominal_scale
    pred["calibrated_effective_stiffness_scale"] = calibrated_scale
    pred["analytical_reference_compliance_m_per_N"] = y_reference_test
    pred["quadratic_nominal_section_compliance_m_per_N"] = y_nominal_fixed_section
    pred["quadratic_calibrated_section_compliance_m_per_N"] = y_calibrated_fixed_section
    pred.to_csv(outdir / "test_predictions.csv", index=False)
    joblib.dump({"quadratic_fea": model, "features": FEATURES,
                 "initial_effective_stiffness_scale": nominal_scale,
                 "calibrated_effective_stiffness_scale": calibrated_scale,
                 "calibration_thickness_scale": calibration_thickness},
                outdir / "quadratic_model.joblib")

    raw = metrics["nominal_quadratic_vs_analytical_geometry_test"]
    corrected = metrics["calibrated_quadratic_vs_analytical_geometry_test"]
    fea = metrics["quadratic_vs_fea_test_at_sweep_inputs"]
    table = r"""% Generated by src/fit_corrected_surrogate.py.
\begin{{table}}[htbp]
\centering
\caption{{Held-out performance for the quadratic beam response surface.}}
\label{{tab:quadratic-validation}}
\begin{{tabular}}{{lrrrr}}
\toprule
Comparison & RMSE & MAE & $R^2$ & MAPE (\%) \\
\midrule
Quadratic vs. FEA at DOE inputs & {f_rmse:.3e} & {f_mae:.3e} & {f_r2:.4f} & {f_mape:.3f} \\
Nominal stiffness vs. analytical & {r_rmse:.3e} & {r_mae:.3e} & {r_r2:.4f} & {r_mape:.3f} \\
Calibrated stiffness vs. analytical & {c_rmse:.3e} & {c_mae:.3e} & {c_r2:.4f} & {c_mape:.3f} \\
\bottomrule
\end{{tabular}}
\end{{table}}
""".format(f_rmse=fea["rmse_m_per_N"], f_mae=fea["mae_m_per_N"], f_r2=fea["r2"],
           f_mape=fea["mean_absolute_percentage_error_percent"],
           r_rmse=raw["rmse_m_per_N"], r_mae=raw["mae_m_per_N"], r_r2=raw["r2"],
           r_mape=raw["mean_absolute_percentage_error_percent"],
           c_rmse=corrected["rmse_m_per_N"], c_mae=corrected["mae_m_per_N"], c_r2=corrected["r2"],
           c_mape=corrected["mean_absolute_percentage_error_percent"])
    (outdir / "results_table.tex").write_text(table)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
