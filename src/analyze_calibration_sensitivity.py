"""Assess sensitivity to the number and choice of analytical calibration cases.

The quadratic response surface is fitted once on the configured training split.
For each calibration count, a centred/maximin set and repeated random subsets
are drawn from that training split only. The same held-out geometries are used
for every comparison and are never used to fit or calibrate the model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from fit_corrected_surrogate import (
    FEATURES,
    analytical_reference,
    fit_effective_stiffness_scale,
    quadratic_model,
    score,
    select_calibration_indices,
)


DEFAULT_COUNTS = (1, 2, 3, 5)


def make_calibration_subsets(
    data: pd.DataFrame,
    train_idx: np.ndarray,
    n_calibration: int,
    repetitions: int,
    seed: int,
) -> list[tuple[str, np.ndarray]]:
    """Return one centred/maximin subset and unique random training subsets."""
    train_idx = np.asarray(train_idx, dtype=int)
    if not 1 <= n_calibration <= len(train_idx):
        raise ValueError("n_calibration must be between 1 and the training size")
    if repetitions < 1:
        raise ValueError("repetitions must be at least 1")

    centre = select_calibration_indices(data, train_idx, n_calibration)
    selections: list[tuple[str, np.ndarray]] = [("centred_maximin", centre)]
    seen = {tuple(sorted(map(int, centre)))}
    rng = np.random.default_rng(seed)
    attempts = 0
    max_attempts = max(100, repetitions * 100)
    while len(selections) < repetitions:
        candidate = np.asarray(
            rng.choice(train_idx, size=n_calibration, replace=False), dtype=int
        )
        key = tuple(sorted(map(int, candidate)))
        if key not in seen:
            seen.add(key)
            selections.append(("random_training_subset", candidate))
        attempts += 1
        if attempts > max_attempts:
            raise ValueError(
                "Could not generate enough unique subsets; reduce repetitions "
                "or increase the training set"
            )
    return selections


def sensitivity_plot(runs: pd.DataFrame, path: Path) -> None:
    """Plot held-out error changes and fitted scale by calibration count."""
    counts = sorted(runs["n_calibration"].unique())
    mape_delta = [
        runs.loc[runs.n_calibration == n, "mape_change_percentage_points"].to_numpy()
        for n in counts
    ]
    rmse_delta = [
        runs.loc[runs.n_calibration == n, "rmse_change_percent"].to_numpy()
        for n in counts
    ]
    scale_values = [
        runs.loc[runs.n_calibration == n, "calibrated_effective_stiffness_scale"].to_numpy()
        for n in counts
    ]
    version = tuple(int(part) for part in matplotlib.__version__.split(".")[:2])
    label_argument = "tick_labels" if version >= (3, 9) else "labels"
    box_labels = {label_argument: counts}

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    axes[0].boxplot(mape_delta, showmeans=True, **box_labels)
    axes[0].axhline(0, color="black", linewidth=1, linestyle="--")
    axes[0].set(title="MAPE change", xlabel="Calibration observations",
                ylabel="Calibrated − nominal (percentage points)")
    axes[1].boxplot(rmse_delta, showmeans=True, **box_labels)
    axes[1].axhline(0, color="black", linewidth=1, linestyle="--")
    axes[1].set(title="RMSE change", xlabel="Calibration observations",
                ylabel="Change from nominal (%)")
    axes[2].boxplot(scale_values, showmeans=True, **box_labels)
    axes[2].axhline(1, color="black", linewidth=1, linestyle="--",
                    label="Reference scale")
    axes[2].set(title="Fitted stiffness scale", xlabel="Calibration observations",
                ylabel=r"Estimated $s_E$")
    axes[2].legend(frameon=False)
    fig.suptitle("Calibration sensitivity on a fixed held-out geometry set")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def summarize_runs(runs: pd.DataFrame) -> pd.DataFrame:
    """Aggregate accuracy and parameter sensitivity by calibration count."""
    summary_rows = []
    for n, group in runs.groupby("n_calibration", sort=True):
        row = {
            "n_calibration": int(n),
            "n_subsets": int(len(group)),
            "scale_mean": float(group.calibrated_effective_stiffness_scale.mean()),
            "scale_sd": float(group.calibrated_effective_stiffness_scale.std(ddof=0)),
            "scale_median": float(group.calibrated_effective_stiffness_scale.median()),
            "scale_min": float(group.calibrated_effective_stiffness_scale.min()),
            "scale_max": float(group.calibrated_effective_stiffness_scale.max()),
            "nominal_mape_mean_percent": float(group.nominal_mape_percent.mean()),
            "calibrated_mape_mean_percent": float(group.calibrated_mape_percent.mean()),
            "mape_change_mean_percentage_points": float(
                group.mape_change_percentage_points.mean()),
            "mape_change_median_percentage_points": float(
                group.mape_change_percentage_points.median()),
            "mape_change_sd_percentage_points": float(
                group.mape_change_percentage_points.std(ddof=0)),
            "fraction_mape_improved": float(
                (group.mape_change_percentage_points < 0).mean()),
            "nominal_rmse_mean_m_per_N": float(group.nominal_rmse_m_per_N.mean()),
            "calibrated_rmse_mean_m_per_N": float(group.calibrated_rmse_m_per_N.mean()),
            "rmse_change_mean_percent": float(group.rmse_change_percent.mean()),
            "rmse_change_median_percent": float(group.rmse_change_percent.median()),
            "fraction_rmse_improved": float((group.rmse_change_percent < 0).mean()),
            "nominal_mae_mean_m_per_N": float(group.nominal_mae_m_per_N.mean()),
            "calibrated_mae_mean_m_per_N": float(group.calibrated_mae_m_per_N.mean()),
            "mae_change_mean_percent": float(group.mae_change_percent.mean()),
            "mae_change_median_percent": float(group.mae_change_percent.median()),
            "fraction_mae_improved": float((group.mae_change_percent < 0).mean()),
            "nominal_max_ape_mean_percent": float(group.nominal_max_ape_percent.mean()),
            "calibrated_max_ape_mean_percent": float(group.calibrated_max_ape_percent.mean()),
            "fraction_max_ape_improved": float(
                (group.calibrated_max_ape_percent < group.nominal_max_ape_percent).mean()),
        }
        summary_rows.append(row)
    return pd.DataFrame(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Repeat limited-observation stiffness calibration without touching holdouts."
    )
    parser.add_argument("--input", required=True, help="Existing PyMAPDL sweep CSV")
    parser.add_argument("--config", default="config.quadratic.example.json")
    parser.add_argument("--layup", default="layup.json")
    parser.add_argument("--output-dir", default="results/calibration_sensitivity")
    parser.add_argument("--counts", nargs="+", type=int, default=list(DEFAULT_COUNTS),
                        help="Calibration observation counts (default: 1 2 3 5)")
    parser.add_argument("--repetitions", type=int, default=30,
                        help="Unique subsets per count, including centred/maximin set")
    parser.add_argument("--seed", type=int, default=None,
                        help="Subset seed; defaults to model seed + 1")
    args = parser.parse_args()

    data = pd.read_csv(args.input)
    config = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text())
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    required = FEATURES + ["case_id", "D_eff_Nm", "S_eff_N_per_m", "compliance_m_per_N"]
    missing = sorted(set(required) - set(data.columns))
    if missing:
        raise ValueError(f"Input CSV missing required columns: {missing}")
    data = data.dropna(subset=required).reset_index(drop=True)
    x = data[FEATURES].to_numpy(float)
    y_fea = data["compliance_m_per_N"].to_numpy(float)
    if len(data) < 42 or np.any(y_fea <= 0):
        raise ValueError("Need at least 42 complete cases and positive compliances")

    model_config = config.get("model", {})
    split_seed = int(model_config.get("seed", 1234))
    train_idx, test_idx = train_test_split(
        np.arange(len(data)),
        test_size=float(model_config.get("test_fraction", 0.2)),
        random_state=split_seed,
    )
    if len(train_idx) < 21:
        raise ValueError("Training set must exceed the 21 coefficients in the quadratic")
    model = quadratic_model().fit(x[train_idx], y_fea[train_idx])

    nominal_scale = float(model_config.get("initial_effective_stiffness_scale", 1.0))
    reference_scale = float(model_config.get("reference_effective_stiffness_scale", 1.0))
    calibration_thickness = float(model_config.get("calibration_thickness_scale", 1.0))
    scale_bounds = tuple(map(float, config["sweep"]["effective_stiffness_scale"]))
    if not scale_bounds[0] < nominal_scale < scale_bounds[1]:
        raise ValueError("initial_effective_stiffness_scale must be inside stiffness bounds")
    if not scale_bounds[0] <= reference_scale <= scale_bounds[1]:
        raise ValueError("reference_effective_stiffness_scale must be inside stiffness bounds")
    if not config["sweep"]["thickness_scale"][0] <= calibration_thickness <= config["sweep"]["thickness_scale"][1]:
        raise ValueError("calibration_thickness_scale must be inside thickness sweep")

    heldout_geometry = data.iloc[test_idx].copy()
    y_reference_test = analytical_reference(
        heldout_geometry, config, layup, calibration_thickness, reference_scale)
    x_test = x[test_idx].copy()
    x_test[:, FEATURES.index("thickness_scale")] = calibration_thickness
    x_test[:, FEATURES.index("effective_stiffness_scale")] = nominal_scale
    nominal_predictions = model.predict(x_test)
    nominal_scores = score(y_reference_test, nominal_predictions)
    ncal_counts = sorted(set(args.counts))
    if not ncal_counts or ncal_counts[0] < 1:
        raise ValueError("All calibration counts must be positive")
    if args.repetitions < 1:
        raise ValueError("--repetitions must be at least 1")
    subset_seed = int(args.seed if args.seed is not None else split_seed + 1)
    rows = []

    for n_cal in ncal_counts:
        subsets = make_calibration_subsets(
            data, train_idx, n_cal, args.repetitions, subset_seed + n_cal)
        for repetition, (selection_method, cal_idx) in enumerate(subsets):
            cal_data = data.iloc[cal_idx]
            y_reference_cal = analytical_reference(
                cal_data, config, layup, calibration_thickness, reference_scale)
            x_cal = x[cal_idx].copy()
            x_cal[:, FEATURES.index("thickness_scale")] = calibration_thickness
            fitted_scale, calibration_sse = fit_effective_stiffness_scale(
                model, x_cal, y_reference_cal, nominal_scale, scale_bounds)

            x_test_calibrated = x_test.copy()
            x_test_calibrated[:, FEATURES.index("effective_stiffness_scale")] = fitted_scale
            calibrated_scores = score(y_reference_test, model.predict(x_test_calibrated))

            row = {
                "n_calibration": n_cal,
                "repetition": repetition,
                "selection_method": selection_method,
                "calibration_case_ids": ";".join(
                    str(int(data.iloc[i]["case_id"])) for i in cal_idx),
                "calibrated_effective_stiffness_scale": fitted_scale,
                "calibration_sse_m_per_N_squared": calibration_sse,
                "nominal_mape_percent": nominal_scores["mean_absolute_percentage_error_percent"],
                "calibrated_mape_percent": calibrated_scores["mean_absolute_percentage_error_percent"],
                "mape_change_percentage_points": (
                    calibrated_scores["mean_absolute_percentage_error_percent"]
                    - nominal_scores["mean_absolute_percentage_error_percent"]),
                "nominal_rmse_m_per_N": nominal_scores["rmse_m_per_N"],
                "calibrated_rmse_m_per_N": calibrated_scores["rmse_m_per_N"],
                "rmse_change_percent": 100 * (
                    calibrated_scores["rmse_m_per_N"] / nominal_scores["rmse_m_per_N"] - 1),
                "nominal_mae_m_per_N": nominal_scores["mae_m_per_N"],
                "calibrated_mae_m_per_N": calibrated_scores["mae_m_per_N"],
                "mae_change_percent": 100 * (
                    calibrated_scores["mae_m_per_N"] / nominal_scores["mae_m_per_N"] - 1),
                "nominal_r2": nominal_scores["r2"],
                "calibrated_r2": calibrated_scores["r2"],
                "r2_change": calibrated_scores["r2"] - nominal_scores["r2"],
                "nominal_max_ape_percent": nominal_scores["max_absolute_percentage_error_percent"],
                "calibrated_max_ape_percent": calibrated_scores["max_absolute_percentage_error_percent"],
            }
            rows.append(row)

    runs = pd.DataFrame(rows)
    summary = summarize_runs(runs)
    runs.to_csv(outdir / "calibration_sensitivity_runs.csv", index=False)
    summary.to_csv(outdir / "calibration_sensitivity_summary.csv", index=False)
    sensitivity_plot(runs, outdir / "calibration_sensitivity.png")

    metadata = {
        "input_file": str(args.input),
        "data_origin": str(data["data_origin"].iloc[0]) if "data_origin" in data else "PyMAPDL SHELL181 FEA",
        "n_total": int(len(data)),
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "split_seed": split_seed,
        "subset_seed": subset_seed,
        "calibration_counts": ncal_counts,
        "subsets_per_count": int(args.repetitions),
        "training_case_ids": [int(data.iloc[i]["case_id"]) for i in train_idx],
        "heldout_case_ids": [int(data.iloc[i]["case_id"]) for i in test_idx],
        "heldout_analytical_nominal_metrics": nominal_scores,
        "calibration_thickness_scale": calibration_thickness,
        "initial_effective_stiffness_scale": nominal_scale,
        "reference_effective_stiffness_scale": reference_scale,
        "effective_stiffness_bounds": scale_bounds,
        "interpretation": (
            "Each fitted scale uses analytical references only at training geometries. "
            "The same held-out geometry set scores every subset and is never used for fitting or calibration."
        ),
    }
    (outdir / "calibration_sensitivity_metrics.json").write_text(
        json.dumps(metadata, indent=2))
    print(summary.to_string(index=False))
    print(f"\nWrote sensitivity outputs to {outdir}")


if __name__ == "__main__":
    main()
