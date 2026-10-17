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

from analytic import tip_compliance_from_stiffness

FEATURES = ["length_m", "root_width_m", "taper_ratio", "thickness_scale",
            "effective_stiffness_scale"]


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


def analytical_reference(data: pd.DataFrame, config: dict) -> np.ndarray:
    kappa = float(config["fea"].get("shear_correction_factor_analytic", 5 / 6))
    return np.asarray([
        tip_compliance_from_stiffness(
            row.length_m, row.root_width_m, row.taper_ratio,
            row.D_eff_Nm, row.S_eff_N_per_m, kappa)
        for row in data.itertuples(index=False)
    ], dtype=float)


def select_calibration_indices(x: np.ndarray, train_idx: np.ndarray, n: int) -> np.ndarray:
    """Select training cases nearest the design-space centre, then spread out."""
    centre_scale = StandardScaler().fit(x[train_idx])
    z = centre_scale.transform(x)
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="CSV from run_laminate_shell_sweep.py")
    ap.add_argument("--config", default="config.quadratic.example.json")
    ap.add_argument("--output-dir", default="results/quadratic_model")
    args = ap.parse_args()

    data = pd.read_csv(args.input)
    config = json.loads(Path(args.config).read_text())
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
    y_quad = model.predict(x[test_idx])

    # The analytic benchmark differs from the response surface by one global
    # stiffness scale: C_corrected = beta*C_quad, lambda_eff = 1/beta.
    # One scalar observation is the mathematical minimum; other test cases
    # remain untouched and measure transfer of that correction.
    ncal = int(mcfg.get("calibration_points", 1))
    if ncal < 1 or ncal > len(train_idx):
        raise ValueError("calibration_points must be between 1 and n_train")
    cal_idx = select_calibration_indices(x, train_idx, ncal)
    y_reference = analytical_reference(data, config)
    q_cal = model.predict(x[cal_idx])
    ref_cal = y_reference[cal_idx]
    beta = float(np.dot(q_cal, ref_cal) / np.dot(q_cal, q_cal))
    if beta <= 0 or not np.isfinite(beta):
        raise ValueError("Least-squares fit returned a non-positive stiffness correction")
    effective_stiffness_factor = 1.0 / beta
    y_corrected = beta * y_quad

    metrics = {
        "input_file": str(args.input),
        "data_origin": str(data["data_origin"].iloc[0]) if "data_origin" in data else "PyMAPDL SHELL181 FEA",
        "n_total": int(len(data)), "n_train": int(len(train_idx)), "n_test": int(len(test_idx)),
        "n_analytical_calibration": int(ncal),
        "calibration_case_ids": [int(data.iloc[i]["case_id"]) for i in cal_idx],
        "features": FEATURES, "target": "tip compliance |delta|/P [m/N]",
        "response_surface": "full quadratic polynomial (21 coefficients for five inputs)",
        "calibration": {
            "method": "one-parameter least-squares fit on training cases",
            "compliance_multiplier_beta": beta,
            "effective_section_stiffness_factor": effective_stiffness_factor,
            "relationship": "C_corrected = beta * C_quadratic; effective stiffness factor = 1/beta",
        },
        "quadratic_vs_fea_test": score(y_fea[test_idx], y_quad),
        "corrected_quadratic_vs_analytical_test": score(y_reference[test_idx], y_corrected),
        "raw_quadratic_vs_analytical_test": score(y_reference[test_idx], y_quad),
    }
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    pred = data.iloc[test_idx].copy()
    pred["analytical_compliance_m_per_N"] = y_reference[test_idx]
    pred["quadratic_fea_compliance_m_per_N"] = y_quad
    pred["corrected_compliance_m_per_N"] = y_corrected
    pred.to_csv(outdir / "test_predictions.csv", index=False)
    joblib.dump({"quadratic_fea": model, "features": FEATURES,
                 "compliance_multiplier_beta": beta,
                 "effective_section_stiffness_factor": effective_stiffness_factor},
                outdir / "quadratic_model.joblib")

    raw, corrected, fea = (metrics["raw_quadratic_vs_analytical_test"],
                           metrics["corrected_quadratic_vs_analytical_test"],
                           metrics["quadratic_vs_fea_test"])
    table = r"""% Generated by src/fit_corrected_surrogate.py.
\begin{{table}}[htbp]
\centering
\caption{{Held-out performance for the quadratic beam response surface.}}
\label{{tab:quadratic-validation}}
\begin{{tabular}}{{lrrrr}}
\toprule
Comparison & RMSE & MAE & $R^2$ & MAPE (\%) \\
\midrule
Quadratic vs. FEA & {f_rmse:.3e} & {f_mae:.3e} & {f_r2:.4f} & {f_mape:.3f} \\
Raw quadratic vs. analytical & {r_rmse:.3e} & {r_mae:.3e} & {r_r2:.4f} & {r_mape:.3f} \\
Corrected quadratic vs. analytical & {c_rmse:.3e} & {c_mae:.3e} & {c_r2:.4f} & {c_mape:.3f} \\
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
