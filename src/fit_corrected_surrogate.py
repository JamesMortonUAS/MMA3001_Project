"""Fit FEA compliance model and sparse analytical discrepancy correction."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from analytic import tip_compliance, tip_compliance_from_stiffness
from laminate_section import effective_bending_shear_stiffness

FEATURES = ["length_m", "root_width_m", "taper_ratio"]


def gp_model(seed: int) -> GaussianProcessRegressor:
    kernel = (ConstantKernel(1.0, (1e-3, 1e3))
              * Matern(length_scale=np.ones(3), length_scale_bounds=(1e-2, 1e2), nu=2.5)
              + WhiteKernel(noise_level=1e-8, noise_level_bounds=(1e-12, 1e-3)))
    return make_pipeline(StandardScaler(), GaussianProcessRegressor(
        kernel=kernel, normalize_y=True, n_restarts_optimizer=3,
        random_state=seed, alpha=1e-10))


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


def analytic_column(df: pd.DataFrame, config: dict, layup: dict | None = None) -> np.ndarray:
    b = config["baseline"]
    kappa = float(config["fea"].get("shear_correction_factor_analytic", 5 / 6))
    if layup is not None:
        d_eff, s_eff = effective_bending_shear_stiffness(layup)
        return np.asarray([tip_compliance_from_stiffness(
            r.length_m, r.root_width_m, r.taper_ratio, d_eff, s_eff, kappa)
            for r in df.itertuples()], dtype=float)
    return np.asarray([tip_compliance(
        r.length_m, r.root_width_m, r.taper_ratio, b["thickness_m"],
        b["youngs_modulus_Pa"], b["poissons_ratio"], kappa)
        for r in df.itertuples()], dtype=float)


def tex_escape(s: str) -> str:
    return s.replace("_", r"\_")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--layup", default=None,
                    help="Optional layup.json; uses laminate CLT effective stiffnesses for analytical reference")
    ap.add_argument("--output-dir", default="results/model")
    args = ap.parse_args()
    data = pd.read_csv(args.input)
    config = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text()) if args.layup else None
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    required = FEATURES + ["compliance_m_per_N"]
    missing = sorted(set(required) - set(data.columns))
    if missing:
        raise ValueError(f"Input CSV missing required columns: {missing}")
    data = data.dropna(subset=required).reset_index(drop=True)
    if len(data) < 16:
        raise ValueError("Need at least 16 complete cases for a train/test assessment")
    if np.any(data["compliance_m_per_N"] <= 0):
        raise ValueError("Compliance values must be positive")

    x = data[FEATURES].to_numpy(float)
    y_fea = data["compliance_m_per_N"].to_numpy(float)
    all_idx = np.arange(len(data))
    mcfg = config.get("model", {})
    train_idx, test_idx = train_test_split(
        all_idx, test_size=float(mcfg.get("test_fraction", 0.25)),
        random_state=int(mcfg.get("seed", 1234)))

    gp = gp_model(int(mcfg.get("seed", 1234)))
    gp.fit(x[train_idx], y_fea[train_idx])
    y_raw_test = gp.predict(x[test_idx])

    # Credible baseline: full quadratic response surface with ridge regularization.
    quad = make_pipeline(StandardScaler(), PolynomialFeatures(degree=2, include_bias=False),
                         StandardScaler(), Ridge(alpha=1e-8))
    quad.fit(x[train_idx], y_fea[train_idx])
    y_quad_test = quad.predict(x[test_idx])

    # Sparse ground-truth calibration points are drawn only from the training set.
    ncal = min(int(mcfg.get("calibration_points", 5)), len(train_idx))
    rng = np.random.default_rng(int(mcfg.get("seed", 1234)) + 1)
    cal_idx = np.sort(rng.choice(train_idx, size=ncal, replace=False))
    y_analytical = analytic_column(data, config, layup)
    y_resid_cal = y_analytical[cal_idx] - y_fea[cal_idx]
    discrepancy = gp_model(int(mcfg.get("seed", 1234)) + 2)
    discrepancy.fit(x[cal_idx], y_resid_cal)
    y_corrected_test = y_raw_test + discrepancy.predict(x[test_idx])

    source_name = (str(data["data_origin"].iloc[0]) if "data_origin" in data
                   else "PyMAPDL FEA (assumed)")
    source_short = "synthetic proxy" if "proxy" in source_name.lower() else "FEA"
    metrics = {
        "input_file": str(args.input),
        "data_origin": source_name,
        "n_total": int(len(data)), "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)), "n_analytical_calibration": int(ncal),
        "features": FEATURES, "target": "tip compliance |delta|/P [m/N]",
        "analytical_reference": ("laminate CLT A/B/D + x-z shear stiffness" if layup else "isotropic rectangular beam"),
        "raw_gp_vs_source_test": score(y_fea[test_idx], y_raw_test),
        "quadratic_vs_source_test": score(y_fea[test_idx], y_quad_test),
        "raw_gp_vs_analytical_test": score(y_analytical[test_idx], y_raw_test),
        "corrected_gp_vs_analytical_test": score(y_analytical[test_idx], y_corrected_test),
    }
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    pred = data.iloc[test_idx].copy()
    pred["split"] = "test"
    pred["analytical_compliance_m_per_N"] = y_analytical[test_idx]
    pred["gp_raw_compliance_m_per_N"] = y_raw_test
    pred["gp_corrected_compliance_m_per_N"] = y_corrected_test
    pred["quadratic_compliance_m_per_N"] = y_quad_test
    pred.to_csv(outdir / "test_predictions.csv", index=False)
    joblib.dump({"fea_gp": gp, "quadratic": quad, "discrepancy_gp": discrepancy,
                 "features": FEATURES}, outdir / "surrogate_models.joblib")

    raw = metrics["raw_gp_vs_analytical_test"]
    corrected = metrics["corrected_gp_vs_analytical_test"]
    baseline = metrics["quadratic_vs_source_test"]
    table = r"""% Generated by src/fit_corrected_surrogate.py; do not edit by hand.
\begin{{table}}[htbp]
\centering
\caption{{Held-out performance for the beam benchmark using {source_short} data. Compliance is $|\delta|/P$.}}
\label{{tab:surrogate-validation}}
\begin{{tabular}}{{lrrrr}}
\toprule
Method & RMSE & MAE & $R^2$ & MAPE (\%) \\
\midrule
Quadratic response surface (against {source_short}) & {b_rmse:.3e} & {b_mae:.3e} & {b_r2:.4f} & {b_mape:.3f} \\
Raw GP surrogate (against analytical reference) & {r_rmse:.3e} & {r_mae:.3e} & {r_r2:.4f} & {r_mape:.3f} \\
Corrected GP (against analytical reference) & {c_rmse:.3e} & {c_mae:.3e} & {c_r2:.4f} & {c_mape:.3f} \\
\bottomrule
\end{{tabular}}
\end{{table}}
""".format(
        b_rmse=baseline["rmse_m_per_N"], b_mae=baseline["mae_m_per_N"],
        b_r2=baseline["r2"], b_mape=baseline["mean_absolute_percentage_error_percent"],
        r_rmse=raw["rmse_m_per_N"], r_mae=raw["mae_m_per_N"], r_r2=raw["r2"],
        r_mape=raw["mean_absolute_percentage_error_percent"],
        c_rmse=corrected["rmse_m_per_N"], c_mae=corrected["mae_m_per_N"],
        c_r2=corrected["r2"], c_mape=corrected["mean_absolute_percentage_error_percent"],
        source_short=source_short)
    (outdir / "results_table.tex").write_text(table)
    print(json.dumps(metrics, indent=2))
    print(f"Wrote fitted models, predictions, metrics, and {outdir / 'results_table.tex'}")


if __name__ == "__main__":
    main()
