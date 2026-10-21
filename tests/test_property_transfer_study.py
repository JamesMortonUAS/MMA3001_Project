import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from laminate_section import effective_bending_shear_stiffness
from run_property_transfer_study import (
    design_properties,
    filter_training_laminates,
    fit_calibrate_validate,
    mechanics_basis,
    make_dataset_a_design,
    make_dataset_b_design,
    scale_material_layup,
    section_record,
    select_calibration_indices,
)


def load_json(name):
    return json.loads((ROOT / name).read_text())


def test_dataset_a_is_full_linearly_spaced_tensor_grid_and_excludes_target_pair():
    config = load_json("config.hidden_property.example.json")
    design = make_dataset_a_design(config)
    assert len(design) == 3**4
    names = ("width_m", "thickness_scale", "in_plane_stiffness_scale",
             "transverse_shear_stiffness_scale")
    for name in names:
        levels = np.sort(design[name].unique())
        assert len(levels) == 3
        assert np.allclose(np.diff(levels), np.diff(levels)[0])
    training = filter_training_laminates(
        design, config["dataset_a"]["excluded_laminates"]
    )
    assert len(training) == 72
    assert not (
        np.isclose(training.in_plane_stiffness_scale, 0.9)
        & np.isclose(training.transverse_shear_stiffness_scale, 1.15)
    ).any()
    assert set(names) <= set(training.columns)


def test_dataset_b_contains_interpolation_and_extrapolation_geometries():
    config = load_json("config.hidden_property.example.json")
    design = make_dataset_b_design(config)
    assert len(design) == 25
    assert (design.geometry_region == "interpolation").sum() == 9
    assert (design.geometry_region == "extrapolation").sum() == 16
    assert set(("width_m", "thickness_scale", "geometry_region")) <= set(design.columns)
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


def test_one_point_calibration_selection_is_centre_case_and_reproducible():
    data = make_dataset_b_design(load_json("config.hidden_property.example.json"))
    first = select_calibration_indices(data, n_points=1, seed=123)
    second = select_calibration_indices(data, n_points=1, seed=123)
    assert np.array_equal(first, second)
    chosen = data.iloc[first[0]]
    assert np.isclose(chosen.width_m, 0.05)
    assert np.isclose(chosen.thickness_scale, 1.0)


def test_one_point_common_scale_calibration_reports_separate_validation_bins(tmp_path):
    config = load_json("config.hidden_property.example.json")
    layup = load_json("layup.json")
    length = config["beam"]["length_m"]
    kappa = config["beam"]["shear_correction_factor"]
    coefficients = np.array([0.98, 1.02])

    training = filter_training_laminates(
        make_dataset_a_design(config), config["dataset_a"]["excluded_laminates"]
    )
    section_rows = []
    for row in training.itertuples(index=False):
        scaled = scale_material_layup(
            layup, row.thickness_scale,
            row.in_plane_stiffness_scale,
            row.transverse_shear_stiffness_scale,
        )
        section_rows.append(section_record(scaled, row.width_m, kappa))
    data_a = pd.concat(
        [training.reset_index(drop=True), pd.DataFrame(section_rows)], axis=1
    )
    data_a["compliance_m_per_N"] = mechanics_basis(data_a, length) @ coefficients

    data_b = make_dataset_b_design(config)
    section_b = design_properties(layup, data_b, 0.85, 0.85, kappa)
    data_b["compliance_m_per_N"] = mechanics_basis(section_b, length) @ coefficients
    path_a, path_b = tmp_path / "a.csv", tmp_path / "b.csv"
    data_a.to_csv(path_a, index=False)
    data_b.to_csv(path_b, index=False)

    metrics = fit_calibrate_validate(
        config, layup, path_a, path_b, tmp_path / "model"
    )
    assert metrics["n_calibration"] == 1
    assert np.isclose(metrics["fitted_common_stiffness_scale"], 0.85, rtol=1e-5)
    grouped = metrics["validation_by_region"]
    assert grouped["interpolation"]["calibrated"]["n"] == 8
    assert grouped["extrapolation"]["calibrated"]["n"] == 16
    assert grouped["interpolation"]["calibrated"]["mape_percent"] < 1e-4
    assert grouped["extrapolation"]["calibrated"]["mape_percent"] < 1e-4
