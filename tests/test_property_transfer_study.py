import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from laminate_section import effective_bending_shear_stiffness
from run_property_transfer_study import (
    make_dataset_a_design,
    make_dataset_b_design,
    scale_material_layup,
    section_record,
    select_calibration_indices,
)


def load_json(name):
    return json.loads((ROOT / name).read_text())


def test_dataset_a_varies_geometry_and_both_material_property_groups():
    config = load_json("config.hidden_property.example.json")
    design = make_dataset_a_design(config)
    assert len(design) == 120
    assert set(("width_m", "thickness_scale", "in_plane_stiffness_scale",
                "transverse_shear_stiffness_scale")) <= set(design.columns)
    assert design["width_m"].between(0.03, 0.07).all()
    assert design["thickness_scale"].between(0.7, 1.3).all()
    assert design["in_plane_stiffness_scale"].between(0.7, 1.3).all()
    assert design["transverse_shear_stiffness_scale"].between(0.7, 1.3).all()


def test_dataset_b_contains_known_geometry_but_no_hidden_material_scales():
    config = load_json("config.hidden_property.example.json")
    design = make_dataset_b_design(config)
    assert len(design) == 12
    assert set(("width_m", "thickness_scale")) <= set(design.columns)
    assert "in_plane_stiffness_scale" not in design.columns
    assert "transverse_shear_stiffness_scale" not in design.columns


def test_layup_property_scales_change_bending_and_shear_separately():
    layup = load_json("layup.json")
    baseline_d, baseline_s = effective_bending_shear_stiffness(layup)
    in_plane = scale_material_layup(layup, 1.0, 1.2, 1.0)
    shear = scale_material_layup(layup, 1.0, 1.0, 1.2)
    d_in, s_in = effective_bending_shear_stiffness(in_plane)
    d_sh, s_sh = effective_bending_shear_stiffness(shear)
    assert np.isclose(d_in / baseline_d, 1.2, rtol=1e-10)
    assert np.isclose(s_in / baseline_s, 1.0, rtol=1e-10)
    assert np.isclose(d_sh / baseline_d, 1.0, rtol=1e-10)
    assert np.isclose(s_sh / baseline_s, 1.2, rtol=1e-10)


def test_section_record_tracks_total_rigidities_and_mass():
    layup = load_json("layup.json")
    section = section_record(layup, width_m=0.05, shear_correction=5 / 6)
    assert np.isclose(section["total_thickness_m"], 0.0035)
    assert section["EI_Nm2"] > 0
    assert section["kGA_N"] > 0
    assert section["areal_mass_kg_m2"] > 0


def test_calibration_subset_is_reproducible_and_spans_geometry():
    import pandas as pd

    data = pd.DataFrame({
        "width_m": [0.035, 0.05, 0.065] * 4,
        "thickness_scale": np.repeat([0.7, 0.9, 1.1, 1.3], 3),
    })
    first = select_calibration_indices(data, n_points=4, seed=123)
    second = select_calibration_indices(data, n_points=4, seed=123)
    assert np.array_equal(first, second)
    assert len(set(first)) == 4
    assert len(np.unique(data.iloc[first]["thickness_scale"])) >= 3
