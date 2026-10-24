import json
from pathlib import Path
import sys
import copy
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from run_extended_transfer_study import (
    _case_signature,
    configure_large_deflection,
    cross_validated_training_metrics,
    fit_models,
    make_designs,
    model_inputs,
    section_data,
    symmetric_stack,
    validate_config,
)


CONFIG = Path(__file__).parents[1] / "config.extended_transfer.example.json"


class ExtendedTransferDesignTests(unittest.TestCase):
  def test_large_deflection_solver_controls_are_set(self):
    cfg = json.loads(CONFIG.read_text())

    class FakeMapdl:
        commands = []
        def run(self, command):
            self.commands.append(command)

    mapdl = FakeMapdl()
    configure_large_deflection(mapdl, cfg)
    self.assertEqual(mapdl.commands, [
        "NLGEOM,ON", "AUTOTS,ON", "NSUBST,10,100,1", "NEQIT,50", "NROPT,FULL"
    ])

  def test_config_validation_rejects_unsupported_small_deflection_or_invalid_geometry(self):
    cfg = json.loads(CONFIG.read_text())
    validate_config(cfg)
    bad_solver = copy.deepcopy(cfg)
    bad_solver["solver"]["large_deflection"] = False
    with self.assertRaisesRegex(ValueError, "large_deflection must be true"):
        validate_config(bad_solver)
    bad_geometry = copy.deepcopy(cfg)
    bad_geometry["targets"]["geometry_extrapolation"]["ply_counts"] = [8, 16, 31, 100]
    with self.assertRaisesRegex(ValueError, "even integers"):
        validate_config(bad_geometry)

  def test_stack_builder_rejects_unknown_family_and_odd_ply_count(self):
    with self.assertRaisesRegex(ValueError, "Unsupported stack kind"):
        symmetric_stack("unidirectional_only", 8)
    with self.assertRaisesRegex(ValueError, "positive, even"):
        symmetric_stack("tw_glass", 7)

  def test_model_inputs_reject_nonpositive_or_nonfinite_values(self):
    import numpy as np
    import pandas as pd
    frame = pd.DataFrame({
        "width_m": [0.05], "total_thickness_m": [0.003],
        "D_eff_Nm_per_m": [100.0], "S_eff_N_per_m": [1.0e6],
    })
    self.assertTrue(np.isfinite(model_inputs(frame)).all())
    frame.loc[0, "width_m"] = 0.0
    with self.assertRaisesRegex(ValueError, "positive and finite"):
        model_inputs(frame)

  def test_factorial_case_counts_and_levels(self):
    cfg = json.loads(CONFIG.read_text())
    a, targets = make_designs(cfg)
    self.assertEqual(len(a), 80)
    self.assertEqual(len(targets), 100)
    self.assertEqual(a["width_m"].nunique(), 5)
    self.assertEqual(a["ply_count"].nunique(), 4)
    self.assertEqual(len(a[["carbon_scale", "glass_scale"]].drop_duplicates()), 4)
    self.assertEqual(targets.groupby("target_set").size().to_dict(), {
        "material_interpolation": 20,
        "material_extrapolation": 20,
        "layup_transfer": 20,
        "geometry_extrapolation": 20,
        "material_extrapolation_v2": 20,
    })
    v2 = targets[targets.target_set == "material_extrapolation_v2"]
    self.assertTrue((v2.resin_to_fibre_ratio_multiplier_truth == 2.0).all())


  def test_symmetric_stacks_and_nominal_thicknesses(self):
    for nplies, expected_mm in ((8, 1.2), (16, 2.4), (24, 3.6), (32, 4.8), (100, 15.0)):
        stack = symmetric_stack("tw_glass", nplies)
        self.assertEqual(stack, stack[::-1])
        self.assertEqual(stack.count("CF_TWILL"), nplies // 2)
        self.assertEqual(stack.count("GF"), nplies // 2)
        self.assertAlmostEqual(stack.count("CF_TWILL") * 0.2 + stack.count("GF") * 0.1, expected_mm)


  def test_design_records_actual_stack_and_extrapolation_axes(self):
    cfg = json.loads(CONFIG.read_text())
    _, targets = make_designs(cfg)
    geometry = targets[targets.target_set == "geometry_extrapolation"]
    self.assertEqual(set(geometry.width_region), {"in_range", "extrapolation"})
    self.assertEqual(set(geometry.ply_count_region), {"in_range", "extrapolation"})
    row = targets[targets.target_set == "layup_transfer"].iloc[0].to_dict()
    props = section_data(cfg, row)
    self.assertEqual(len(row["ply_sequence_bottom_to_top"].split("-")), row["ply_count"])
    self.assertGreater(props["total_thickness_m"], 0)
    self.assertGreater(props["D_eff_Nm_per_m"], 0)
    self.assertGreater(props["S_eff_N_per_m"], 0)

  def test_fea_case_signature_changes_with_analysis_settings(self):
    cfg = json.loads(CONFIG.read_text())
    a, _ = make_designs(cfg)
    row = a.iloc[0].to_dict()
    signature = _case_signature(cfg, row)
    changed = copy.deepcopy(cfg)
    changed["solver"]["maximum_substeps"] = 200
    self.assertNotEqual(signature, _case_signature(changed, row))
    changed = copy.deepcopy(cfg)
    changed["resin_model"]["matrix_modulus_Pa"] *= 1.1
    self.assertNotEqual(signature, _case_signature(changed, row))

  def test_resin_ratio_proxy_changes_truth_stiffness_but_prediction_uses_nominal(self):
    cfg = json.loads(CONFIG.read_text())
    _, targets = make_designs(cfg)
    row = targets[targets.target_set == "material_extrapolation_v2"].iloc[0].to_dict()
    truth = section_data(cfg, row)
    nominal = section_data(cfg, row, carbon_scale=2.25, glass_scale=0.75,
                           resin_ratio_multiplier=1.0)
    self.assertEqual(truth["resin_to_fibre_ratio_multiplier"], 2.0)
    self.assertLess(truth["D_eff_Nm_per_m"], nominal["D_eff_Nm_per_m"])

  def test_model_fit_and_cross_validation_record_runtime(self):
    import numpy as np
    import pandas as pd
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    cfg = json.loads(CONFIG.read_text())
    cfg["models"]["gp_restarts"] = 0
    n = 20
    width = np.linspace(0.02, 0.10, n)
    thickness = np.linspace(0.0012, 0.0048, n)
    d_eff = np.linspace(1.0e3, 8.0e3, n)
    s_eff = np.linspace(1.0e6, 5.0e6, n)
    ei = width * d_eff
    kga = cfg["beam"]["shear_correction_factor"] * width * s_eff
    compliance = 0.7 * cfg["beam"]["length_m"]**3 / (3 * ei) + 0.3 * cfg["beam"]["length_m"] / kga
    frame = pd.DataFrame({
        "width_m": width, "total_thickness_m": thickness,
        "D_eff_Nm_per_m": d_eff, "S_eff_N_per_m": s_eff,
        "EI_Nm2": ei, "kGA_N": kga, "compliance_m_per_N": compliance,
    })
    fit_times = {}
    cv_times = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        models = fit_models(frame, cfg, fit_times)
        cross_validated_training_metrics(frame, cfg, models, cv_times)
    self.assertEqual(set(fit_times), set(models))
    self.assertEqual(set(cv_times), set(models))
    self.assertTrue(all(np.isfinite(v) and v >= 0 for v in [*fit_times.values(), *cv_times.values()]))


if __name__ == "__main__":
    unittest.main()
