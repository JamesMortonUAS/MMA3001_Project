"""Analytical Euler-Bernoulli + Timoshenko compliance for width taper."""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss


def tip_compliance(length_m: float, root_width_m: float, taper_ratio: float,
                   thickness_m: float, youngs_modulus_pa: float,
                   poissons_ratio: float, shear_correction: float = 5.0 / 6.0,
                   quadrature_order: int = 64) -> float:
    """Return tip compliance |delta/P| [m/N] for a rectangular cantilever.

    Assumes a clamped root, free tip, transverse tip point force, constant
    thickness and isotropic linear-elastic material. The width varies linearly
    from ``root_width_m`` to ``taper_ratio * root_width_m``.
    """
    vals = (length_m, root_width_m, thickness_m, youngs_modulus_pa,
            shear_correction)
    if any(not np.isfinite(v) or v <= 0 for v in vals):
        raise ValueError("length, width, thickness, E and kappa must be positive")
    if not (0 < taper_ratio <= 1):
        raise ValueError("taper_ratio must be in (0, 1]")
    if not (-1 < poissons_ratio < 0.5):
        raise ValueError("isotropic Poisson ratio must be in (-1, 0.5)")
    if quadrature_order < 2:
        raise ValueError("quadrature_order must be >= 2")

    # For a rectangular isotropic section: D_eff = E*t^3/12 and S_eff = G*t.
    g_pa = youngs_modulus_pa / (2.0 * (1.0 + poissons_ratio))
    d_eff = youngs_modulus_pa * thickness_m**3 / 12.0  # N m, per unit width
    s_eff = g_pa * thickness_m  # N/m, per unit width
    return tip_compliance_from_stiffness(
        length_m, root_width_m, taper_ratio, d_eff, s_eff,
        shear_correction, quadrature_order)


def tip_compliance_from_stiffness(length_m: float, root_width_m: float,
                                  taper_ratio: float, d_eff_Nm: float,
                                  s_eff_N_per_m: float,
                                  shear_correction: float = 5.0 / 6.0,
                                  quadrature_order: int = 64) -> float:
    """Return compliance from laminate stiffnesses per unit width.

    ``d_eff_Nm`` is bending stiffness per unit width (N m); ``s_eff_N_per_m``
    is transverse shear stiffness per unit width (N/m). This form can accept
    laminate-equivalent stiffnesses when those have been independently derived.
    """
    vals = (length_m, root_width_m, d_eff_Nm, s_eff_N_per_m, shear_correction)
    if any(not np.isfinite(v) or v <= 0 for v in vals):
        raise ValueError("length, width, stiffnesses and kappa must be positive")
    if not (0 < taper_ratio <= 1):
        raise ValueError("taper_ratio must be in (0, 1]")
    if quadrature_order < 2:
        raise ValueError("quadrature_order must be >= 2")
    xi, wi = leggauss(quadrature_order)
    x = 0.5 * length_m * (xi + 1.0)
    w = root_width_m * (1.0 - (1.0 - taper_ratio) * x / length_m)
    integrand = ((length_m - x)**2 / (w * d_eff_Nm)
                 + 1.0 / (shear_correction * w * s_eff_N_per_m))
    return float(0.5 * length_m * np.dot(wi, integrand))


def beam_mass_kg(length_m: float, root_width_m: float, taper_ratio: float,
                 thickness_m: float, density_kg_m3: float) -> float:
    """Mass for constant-thickness, uniform-density linearly tapered beam."""
    return float(density_kg_m3 * thickness_m * length_m * root_width_m
                 * (1.0 + taper_ratio) / 2.0)
