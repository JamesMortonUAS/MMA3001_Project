"""Create the specified 20-ply shell section in a PyMAPDL session.

The function intentionally requires cured-ply properties as inputs: textile
areal mass alone is not enough to define a structural material card.
"""
from __future__ import annotations

import math
import numpy as np


PLY_SEQUENCE = tuple(
    ("CF_TWILL_200GSM", 0.0) if ply % 2 == 1 else ("GF_100GSM", 0.0)
    for ply in range(1, 21)
)

MATERIAL_IDS = {"CF_TWILL_200GSM": 1, "GF_100GSM": 2}
ORTHOTROPIC_FIELDS = {
    "EX_Pa": "EX", "EY_Pa": "EY", "EZ_Pa": "EZ",
    "PRXY": "PRXY", "PRYZ": "PRYZ", "PRXZ": "PRXZ",
    "GXY_Pa": "GXY", "GYZ_Pa": "GYZ", "GXZ_Pa": "GXZ",
    "density_kg_m3": "DENS",
}


def layup_plies(layup: dict) -> list[tuple[str, float]]:
    """Normalize project JSON ply records to (material, angle) tuples."""
    out = []
    for record in layup["plies"]:
        if isinstance(record, dict):
            out.append((str(record["material"]), float(record["orientation_deg"])))
        else:
            name, angle = record
            out.append((str(name), float(angle)))
    return out


def laminate_thickness(layup: dict) -> float:
    """Sum cured ply thicknesses from the layup definition, in metres."""
    cards = layup["material_cards"]
    return sum(float(cards[name]["ply_thickness_m"])
               for name, _angle in layup_plies(layup))


def define_twenty_ply_shell_section(mapdl, material_cards: dict,
                                    section_id: int = 1,
                                    integration_points_per_ply: int = 3) -> float:
    """Define material cards and an SHELL181 section; return total thickness.

    ``material_cards`` must have entries ``CF_TWILL_200GSM`` and ``GF_100GSM``.
    Each entry requires ``ply_thickness_m`` plus the orthotropic constants
    EX/EY/EZ, PRXY/PRYZ/PRXZ, GXY/GYZ/GXZ and cured-composite density.
    Angles in ``PLY_SEQUENCE`` are relative to the shell's local in-plane axis 1.
    Confirm that axis 1 is aligned with the beam axis for the mesh in use.
    """
    if integration_points_per_ply not in (1, 3, 5, 7, 9):
        raise ValueError("Use 1, 3, 5, 7 or 9 through-thickness integration points")
    missing_materials = set(MATERIAL_IDS) - set(material_cards)
    if missing_materials:
        raise ValueError(f"Missing material cards: {sorted(missing_materials)}")

    for name, mat_id in MATERIAL_IDS.items():
        card = material_cards[name]
        required = {"ply_thickness_m", *ORTHOTROPIC_FIELDS}
        missing = required - set(card)
        if missing:
            raise ValueError(f"{name} missing properties: {sorted(missing)}")
        vals = [float(card[k]) for k in required]
        if not all(math.isfinite(v) for v in vals):
            raise ValueError(f"{name} has missing/non-finite material values")
        if float(card["ply_thickness_m"]) <= 0:
            raise ValueError(f"{name} ply_thickness_m must be positive")
        for key in ("EX_Pa", "EY_Pa", "EZ_Pa", "GXY_Pa", "GYZ_Pa", "GXZ_Pa", "density_kg_m3"):
            if float(card[key]) <= 0:
                raise ValueError(f"{name} {key} must be positive")
        # Orthotropic Poisson ratios do not share the isotropic (-1, 0.5)
        # bound; ensure finiteness here and let MAPDL check matrix stability.
        for key, label in ORTHOTROPIC_FIELDS.items():
            mapdl.mp(label, mat_id, float(card[key]))

    mapdl.et(1, "SHELL181")
    mapdl.keyopt(1, 8, 1)  # store layer output for postprocessing
    mapdl.sectype(section_id, "SHELL")
    total_thickness = 0.0
    for material_name, angle_deg in PLY_SEQUENCE:
        card = material_cards[material_name]
        tply = float(card["ply_thickness_m"])
        mapdl.secdata(tply, MATERIAL_IDS[material_name], angle_deg,
                      integration_points_per_ply)
        total_thickness += tply
    return total_thickness


def effective_bending_shear_stiffness(layup: dict) -> tuple[float, float]:
    """Return CLT effective D11 [N m] and x-z shear stiffness [N/m].

    For the all-zero-degree stack, this assembles the A/B/D matrices using
    each ply's in-plane orthotropic constants. D_eff is obtained by inverting
    the full ABD matrix with zero membrane resultants and zero non-x moments,
    so extension-bending coupling in the unsymmetric stack is retained.
    S_eff is the through-thickness sum of G13*t for the beam's x-z shear.
    """
    cards = layup["material_cards"]
    plies = layup_plies(layup)
    h = sum(float(cards[name]["ply_thickness_m"]) for name, _ in plies)
    z = -0.5 * h
    abd = np.zeros((6, 6), dtype=float)
    shear_xz = 0.0
    for name, angle in plies:
        if abs(float(angle)) > 1e-12:
            raise NotImplementedError("CLT stiffness helper currently assumes all plies at 0 degrees")
        p = cards[name]
        t = float(p["ply_thickness_m"])
        e1, e2, nu12, g12 = (float(p[k]) for k in ("EX_Pa", "EY_Pa", "PRXY", "GXY_Pa"))
        nu21 = nu12 * e2 / e1
        den = 1.0 - nu12 * nu21
        if den <= 0:
            raise ValueError(f"Invalid in-plane elastic constants for {name}")
        q = np.zeros((3, 3), dtype=float)
        q[0, 0] = e1 / den
        q[1, 1] = e2 / den
        q[0, 1] = q[1, 0] = nu12 * e2 / den
        q[2, 2] = g12
        z_next = z + t
        abd[:3, :3] += q * t
        abd[:3, 3:] += q * (z_next**2 - z**2) / 2.0
        abd[3:, :3] += q * (z_next**2 - z**2) / 2.0
        abd[3:, 3:] += q * (z_next**3 - z**3) / 3.0
        shear_xz += float(p["GXZ_Pa"]) * t
        z = z_next
    compliance_kx_per_moment = float(np.linalg.inv(abd)[3, 3])
    if compliance_kx_per_moment <= 0 or shear_xz <= 0:
        raise ValueError("Laminate stiffness matrix is not positive definite")
    return 1.0 / compliance_kx_per_moment, shear_xz
