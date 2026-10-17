"""Pipeline demonstration using a low-fidelity numerical proxy, not FEA."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

from analytic import tip_compliance


def midpoint_proxy(length, b0, rb, t, youngs, poisson, kappa, cells=4):
    """Coarse midpoint integration proxy; intentionally not a finite element solve."""
    g = youngs / (2 * (1 + poisson))
    d = youngs * t**3 / 12
    s = g * t
    dx = length / cells
    x = (np.arange(cells) + 0.5) * dx
    width = b0 * (1 - (1 - rb) * x / length)
    terms = (length - x)**2 / (width * d) + 1 / (kappa * width * s)
    return float(dx * np.sum(terms))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.example.json")
    ap.add_argument("--output", default="results/synthetic_proxy.csv")
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    b, sw, fea = cfg["baseline"], cfg["sweep"], cfg["fea"]
    n = int(sw["n_design_points"])
    u = qmc.LatinHypercube(d=3, seed=int(sw["seed"])).random(n)
    rows = []
    for i, a in enumerate(u):
        L = b["length_m"] * (sw["length_scale"][0] + a[0] * (sw["length_scale"][1] - sw["length_scale"][0]))
        w0 = b["root_width_m"] * (sw["root_width_scale"][0] + a[1] * (sw["root_width_scale"][1] - sw["root_width_scale"][0]))
        rb = sw["taper_ratio"][0] + a[2] * (sw["taper_ratio"][1] - sw["taper_ratio"][0])
        low = midpoint_proxy(L, w0, rb, b["thickness_m"], b["youngs_modulus_Pa"],
                             b["poissons_ratio"], fea.get("shear_correction_factor_analytic", 5/6))
        exact = tip_compliance(L, w0, rb, b["thickness_m"], b["youngs_modulus_Pa"],
                               b["poissons_ratio"], fea.get("shear_correction_factor_analytic", 5/6))
        rows.append({"case_id": i, "length_m": L, "root_width_m": w0,
                     "thickness_m": b["thickness_m"], "taper_ratio": rb,
                     "tip_load_N": b["tip_load_N"], "youngs_modulus_Pa": b["youngs_modulus_Pa"],
                     "poissons_ratio": b["poissons_ratio"], "compliance_m_per_N": low,
                     "reference_exact_compliance_m_per_N": exact,
                     "data_origin": "synthetic four-cell midpoint proxy; NOT FEA"})
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Wrote {len(rows)} synthetic proxy cases to {out}; these are not Ansys results.")


if __name__ == "__main__":
    main()
