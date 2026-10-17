"""Minimize laminate-beam mass subject to a corrected deflection limit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
from scipy.optimize import minimize

from analytic import tip_compliance_from_stiffness
from laminate_section import (effective_bending_shear_stiffness, layup_plies,
                              scaled_layup)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.quadratic.example.json")
    ap.add_argument("--model", default="results/quadratic_model/quadratic_model.joblib")
    ap.add_argument("--layup", default="layup.json")
    ap.add_argument("--target-deflection-m", type=float, required=True)
    ap.add_argument("--minimum-tip-width-m", type=float, default=0.0)
    ap.add_argument("--output", default="results/mass_optimum.json")
    args = ap.parse_args()
    if args.target_deflection_m <= 0 or args.minimum_tip_width_m < 0:
        raise ValueError("Deflection target must be positive; width limit nonnegative")

    cfg = json.loads(Path(args.config).read_text())
    layup = json.loads(Path(args.layup).read_text())
    model_data = joblib.load(args.model)
    model = model_data["quadratic_fea"]
    beta = float(model_data["compliance_multiplier_beta"])
    base, sweep = cfg["baseline"], cfg["sweep"]
    length, load = float(base["length_m"]), float(base["tip_load_N"])
    nominal_stiffness_scale = 1.0
    bounds = [
        (float(sweep["root_width_scale"][0]) * float(base["root_width_m"]),
         float(sweep["root_width_scale"][1]) * float(base["root_width_m"])),
        tuple(map(float, sweep["taper_ratio"])),
        tuple(map(float, sweep["thickness_scale"])),
    ]
    nominal_areal_mass = sum(
        float(layup["material_cards"][mat]["density_kg_m3"])
        * float(layup["material_cards"][mat]["ply_thickness_m"])
        for mat, _ in layup_plies(layup))

    def features(z):
        return np.array([[length, float(z[0]), float(z[1]), float(z[2]),
                          nominal_stiffness_scale]])

    def compliance(z):
        return beta * float(model.predict(features(z))[0])

    def mass(z):
        b0, rb, thickness_scale = map(float, z)
        return nominal_areal_mass * thickness_scale * length * b0 * (1 + rb) / 2

    constraints = [
        {"type": "ineq", "fun": lambda z: args.target_deflection_m - load * compliance(z)},
        {"type": "ineq", "fun": lambda z: z[0] * z[1] - args.minimum_tip_width_m},
    ]
    z0 = np.array([base["root_width_m"], base["taper_ratio"], 1.0])
    result = minimize(mass, z0, method="SLSQP", bounds=bounds, constraints=constraints,
                      options={"ftol": 1e-12, "maxiter": 500})
    if not result.success:
        raise RuntimeError(f"Optimization failed: {result.message}")

    b0, rb, thickness_scale = map(float, result.x)
    check_layup = scaled_layup(layup, thickness_scale, nominal_stiffness_scale)
    d_eff, s_eff = effective_bending_shear_stiffness(check_layup)
    analytical_c = tip_compliance_from_stiffness(
        length, b0, rb, d_eff, s_eff,
        float(cfg["fea"].get("shear_correction_factor_analytic", 5 / 6)))
    output = {
        "status": "success",
        "warning": "Valid only within the fitted benchmark domain and model assumptions.",
        "length_m_fixed": length,
        "root_width_m": b0,
        "tip_width_m": b0 * rb,
        "taper_ratio": rb,
        "thickness_scale": thickness_scale,
        "laminate_thickness_m": sum(
            float(check_layup["material_cards"][name]["ply_thickness_m"])
            for name, _ in layup_plies(check_layup)),
        "effective_stiffness_scale": nominal_stiffness_scale,
        "predicted_deflection_m_corrected_model": load * compliance(result.x),
        "predicted_deflection_m_analytical_check": load * analytical_c,
        "maximum_allowed_deflection_m": args.target_deflection_m,
        "mass_kg": mass(result.x),
        "optimizer_message": str(result.message),
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2))
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
