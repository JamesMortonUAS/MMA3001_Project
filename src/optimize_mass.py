"""Minimize benchmark beam mass subject to a maximum-deflection limit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
from scipy.optimize import minimize

from analytic import beam_mass_kg, tip_compliance, tip_compliance_from_stiffness
from laminate_section import effective_bending_shear_stiffness, layup_plies


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--model", default="results/model/surrogate_models.joblib")
    ap.add_argument("--layup", default=None,
                    help="Optional layup.json for laminate mass and analytical check")
    ap.add_argument("--target-deflection-m", type=float, required=True,
                    help="Maximum permitted tip deflection under configured tip load")
    ap.add_argument("--minimum-tip-width-m", type=float, default=0.0)
    ap.add_argument("--output", default="results/mass_optimum.json")
    args = ap.parse_args()
    if args.target_deflection_m <= 0 or args.minimum_tip_width_m < 0:
        raise ValueError("Deflection target must be positive; width limit nonnegative")

    cfg = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text()) if args.layup else None
    b = cfg["baseline"]
    sw = cfg["sweep"]
    models = joblib.load(args.model)
    length = float(b["length_m"])
    load = float(b["tip_load_N"])
    thickness = (sum(float(layup["material_cards"][n]["ply_thickness_m"])
                     for n, _ in layup_plies(layup)) if layup else float(b["thickness_m"]))
    areal_mass = (sum(float(layup["material_cards"][n]["density_kg_m3"])
                      * float(layup["material_cards"][n]["ply_thickness_m"])
                      for n, _ in layup_plies(layup)) if layup else None)
    b_bounds = [float(sw["root_width_scale"][0] * b["root_width_m"]),
                float(sw["root_width_scale"][1] * b["root_width_m"])]
    rb_bounds = [float(v) for v in sw["taper_ratio"]]

    def features(z):
        return np.array([[length, float(z[0]), float(z[1])]])

    def corrected_compliance(z):
        x = features(z)
        return float(models["fea_gp"].predict(x)[0]
                     + models["discrepancy_gp"].predict(x)[0])

    def objective(z):
        return (areal_mass * length * z[0] * (1.0 + z[1]) / 2.0 if layup
                else beam_mass_kg(length, z[0], z[1], thickness, float(b["density_kg_m3"])))

    constraints = [
        {"type": "ineq", "fun": lambda z: args.target_deflection_m
         - load * corrected_compliance(z)},
        {"type": "ineq", "fun": lambda z: z[0] * z[1]
         - args.minimum_tip_width_m},
    ]
    z0 = np.array([b["root_width_m"], b["taper_ratio"]], dtype=float)
    result = minimize(objective, z0, method="SLSQP", bounds=[b_bounds, rb_bounds],
                      constraints=constraints, options={"ftol": 1e-12, "maxiter": 500})
    if not result.success:
        raise RuntimeError(f"Optimization failed: {result.message}")

    b0, rb = map(float, result.x)
    if layup:
        d_eff, s_eff = effective_bending_shear_stiffness(layup)
        exact_c = tip_compliance_from_stiffness(
            length, b0, rb, d_eff, s_eff,
            float(cfg["fea"].get("shear_correction_factor_analytic", 5/6)))
    else:
        exact_c = tip_compliance(length, b0, rb, thickness,
                                 float(b["youngs_modulus_Pa"]),
                                 float(b["poissons_ratio"]),
                                 float(cfg["fea"].get("shear_correction_factor_analytic", 5/6)))
    output = {
        "status": "success",
        "warning": "Only valid within the trained benchmark domain and model assumptions.",
        "length_m_fixed": length,
        "root_width_m": b0,
        "tip_width_m": b0 * rb,
        "taper_ratio": rb,
        "laminate_thickness_m": thickness,
        "predicted_deflection_m_corrected_model": load * corrected_compliance(result.x),
        "predicted_deflection_m_analytical_check": load * exact_c,
        "maximum_allowed_deflection_m": args.target_deflection_m,
        "mass_kg": objective(result.x),
        "optimizer_message": str(result.message),
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2))
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
