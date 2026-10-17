import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from laminate_section import (effective_bending_shear_stiffness, laminate_thickness,
                              layup_plies, scaled_layup)


def load_layup():
    return json.loads((ROOT / "layup.json").read_text())


def test_layup_has_twenty_alternating_zero_degree_plies():
    plies = layup_plies(load_layup())
    assert len(plies) == 20
    assert plies[0][0] == "CF_TWILL_200GSM"
    assert plies[1][0] == "GF_100GSM"
    assert all(angle == 0.0 for _, angle in plies)


def test_laminate_thickness_matches_specification():
    assert abs(laminate_thickness(load_layup()) - 0.0035) < 1e-12


def test_laminate_effective_stiffness_is_positive():
    d_eff, s_eff = effective_bending_shear_stiffness(load_layup())
    assert d_eff > 0
    assert s_eff > 0


def test_thickness_and_effective_stiffness_scales_update_section_properties():
    baseline = load_layup()
    changed = scaled_layup(baseline, thickness_scale=1.2, stiffness_scale=0.9)
    d0, s0 = effective_bending_shear_stiffness(baseline)
    d1, s1 = effective_bending_shear_stiffness(changed)
    assert abs(laminate_thickness(changed) / laminate_thickness(baseline) - 1.2) < 1e-12
    assert d1 > 0 and s1 > 0
    assert d1 != d0 and s1 != s0


def test_new_design_contains_thickness_and_effective_stiffness_sweep():
    from run_laminate_shell_sweep import make_design

    cfg = json.loads((ROOT / "config.quadratic.example.json").read_text())
    design = make_design(cfg)
    assert len(design) == 120
    assert set(("thickness_scale", "effective_stiffness_scale")) <= set(design.columns)
    assert design["tip_load_N"].eq(10.0).all()


def test_least_squares_recovers_shared_effective_stiffness_scale():
    import numpy as np
    from fit_corrected_surrogate import fit_effective_stiffness_scale

    class InverseStiffnessResponse:
        def predict(self, x):
            geometry_response = (0.01 * x[:, 0] + 0.02 * x[:, 1] + 0.03 * x[:, 2])
            return geometry_response / x[:, 4]

    x_cal = np.array([
        [0.4, 0.04, 0.5, 1.0, 0.8],
        [0.6, 0.06, 0.8, 1.0, 1.2],
    ])
    truth_scale = 1.08
    target = InverseStiffnessResponse().predict(
        np.column_stack([x_cal[:, :4], np.full(len(x_cal), truth_scale)]))
    fitted, sse = fit_effective_stiffness_scale(
        InverseStiffnessResponse(), x_cal, target, 1.0, (0.8, 1.2))
    assert abs(fitted - truth_scale) < 1e-6
    assert sse < 1e-15
