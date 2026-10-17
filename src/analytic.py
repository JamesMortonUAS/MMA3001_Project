"""Analytical small-deflection response for a linearly width-tapered beam."""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss


def tip_compliance_from_stiffness(length_m: float, root_width_m: float,
                                  taper_ratio: float, d_eff_Nm: float,
                                  s_eff_N_per_m: float,
                                  shear_correction: float = 5.0 / 6.0,
                                  quadrature_order: int = 64) -> float:
    """Return tip compliance from laminate bending/shear stiffnesses.

    ``d_eff_Nm`` is bending stiffness per unit width (N m); ``s_eff_N_per_m``
    is transverse shear stiffness per unit width (N/m). Assumes a clamped root,
    a free tip, a transverse tip point load, and small-deflection response.
    """
    vals = (length_m, root_width_m, d_eff_Nm, s_eff_N_per_m, shear_correction)
    if any(not np.isfinite(v) or v <= 0 for v in vals):
        raise ValueError("length, width, stiffnesses and shear factor must be positive")
    if not (0 < taper_ratio <= 1):
        raise ValueError("taper_ratio must be in (0, 1]")
    if quadrature_order < 2:
        raise ValueError("quadrature_order must be >= 2")
    xi, wi = leggauss(quadrature_order)
    x = 0.5 * length_m * (xi + 1.0)
    width = root_width_m * (1.0 - (1.0 - taper_ratio) * x / length_m)
    integrand = ((length_m - x)**2 / (width * d_eff_Nm)
                 + 1.0 / (shear_correction * width * s_eff_N_per_m))
    return float(0.5 * length_m * np.dot(wi, integrand))
