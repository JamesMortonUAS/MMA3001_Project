import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from laminate_section import effective_bending_shear_stiffness, laminate_thickness, layup_plies


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
