"""Generate summary figures for one laminate sweep and its held-out fit."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sweep", required=True, help="Laminate sweep CSV")
    parser.add_argument("--predictions", required=True, help="Held-out predictions CSV")
    parser.add_argument("--output-dir", default="results/figures")
    args = parser.parse_args()

    sweep = pd.read_csv(args.sweep)
    pred = pd.read_csv(args.predictions)
    sweep_required = {"length_m", "root_width_m", "taper_ratio", "compliance_m_per_N",
                      "mass_kg", "tip_deflection_m", "case_id"}
    pred_required = {"analytical_reference_compliance_m_per_N",
                     "quadratic_nominal_section_compliance_m_per_N",
                     "quadratic_calibrated_section_compliance_m_per_N"}
    if missing := sweep_required - set(sweep.columns):
        raise ValueError(f"Sweep CSV is missing columns: {sorted(missing)}")
    if missing := pred_required - set(pred.columns):
        raise ValueError(f"Prediction CSV is missing columns: {sorted(missing)}")

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    ax = axes[0, 0]
    pts = ax.scatter(sweep["length_m"], sweep["compliance_m_per_N"],
                     c=sweep["taper_ratio"], cmap="viridis", s=40, alpha=0.9)
    ax.set(title="Compliance across geometry sweep", xlabel="Length (m)",
           ylabel="Tip compliance (m/N)")
    fig.colorbar(pts, ax=ax, label="Taper ratio")

    ax = axes[0, 1]
    ax.scatter(sweep["mass_kg"], sweep["compliance_m_per_N"],
               c=sweep["root_width_m"], cmap="plasma", s=40, alpha=0.9)
    ax.set(title="Mass–compliance trade-off", xlabel="Beam mass (kg)",
           ylabel="Tip compliance (m/N)")

    ax = axes[1, 0]
    x = pred["analytical_reference_compliance_m_per_N"]
    y_raw = pred["quadratic_nominal_section_compliance_m_per_N"]
    y_corrected = pred["quadratic_calibrated_section_compliance_m_per_N"]
    low = min(x.min(), y_raw.min(), y_corrected.min())
    high = max(x.max(), y_raw.max(), y_corrected.max())
    ax.plot([low, high], [low, high], color="0.3", linestyle="--", label="Perfect agreement")
    ax.scatter(x, y_raw, label="Nominal stiffness", marker="o", s=42)
    ax.scatter(x, y_corrected, label="Calibrated stiffness", marker="x", s=48)
    ax.set(title="Geometry transfer at fixed section properties",
           xlabel="Analytical compliance (m/N)", ylabel="Predicted compliance (m/N)")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ratio = sweep["tip_deflection_m"].abs() / sweep["length_m"]
    ax.scatter(sweep["case_id"], ratio, s=34)
    ax.axhline(0.01, color="tab:green", linestyle="--", label="1% length reference")
    ax.axhline(1.0, color="tab:red", linestyle=":", label="Deflection = length")
    ax.set(title="Linearized deflection / length", xlabel="Sweep case ID",
           ylabel="|Tip deflection| / length")
    ax.legend(frameon=False)

    for ax in axes.flat:
        ax.grid(True, alpha=0.25)
    fig.suptitle("MMA3001 composite beam benchmark", fontsize=15)
    fig.savefig(output / "laminate_sweep_summary.png", dpi=220)
    plt.close(fig)
    print(f"Saved {output / 'laminate_sweep_summary.png'}")


if __name__ == "__main__":
    main()
