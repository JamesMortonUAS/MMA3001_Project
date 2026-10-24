# Extended transfer experiment: run guide

This study contains 80 Dataset A simulations and 100 target simulations across five sets. Three calibration repeats reuse target cases; they do not create additional FEA jobs. The current runner uses a nonlinear large-deflection static formulation. Earlier CSVs in `results/property_transfer_extended/` were generated with a linear-static formulation and four target sets; they are historical baseline data only. The new default output directory is `results/property_transfer_large_deflection_5sets/`, and the runner refuses to resume output whose case signature does not match the current FEA inputs.

## Design

Dataset A is a full 5 x 4 x 4 factorial grid:

- Width: 0.02, 0.04, 0.06, 0.08 and 0.10 m.
- Total ply count: 8, 16, 24 and 32.
- Carbon/glass modulus scale pairs: (1,1), (1,2), (2,1) and (2,2).

The base half-stack repeats carbon twill (T) and glass (G), then mirrors around the mid-plane. Carbon twill plies are 0.20 mm; glass plies are 0.10 mm. These thicknesses produce the following symmetric sequences and totals:

| Plies | Bottom-to-top sequence | T/G counts | Total thickness |
|---:|---|---:|---:|
| 8 | T-G-T-G-G-T-G-T | 4/4 | 1.2 mm |
| 16 | (T-G repeated 4 times) + reverse | 8/8 | 2.4 mm |
| 24 | (T-G repeated 6 times) + reverse | 12/12 | 3.6 mm |
| 32 | (T-G repeated 8 times) + reverse | 16/16 | 4.8 mm |

The carbon scale multiplies all six elastic/shear moduli in each carbon ply material card; the glass scale does the same for glass. Density, Poisson ratios and cured ply thickness do not change. The response model receives calculated laminate bending and transverse-shear properties, not just a single flexural modulus.

There are five 20-case target sets. Sets 1–3 and 5 use the same 5 widths and 4 ply counts as Dataset A, except for the geometry-extrapolation set:

| Target set | Material / geometry change |
|---|---|
| `material_interpolation` | Twill/glass stack at carbon/glass factors (1.5, 1.5). |
| `material_extrapolation` | Twill/glass stack at factors (2.25, 0.75). |
| `layup_transfer` | Symmetric repeating CF unidirectional (U), CF twill (T), E-glass (G) half-stack, mirrored; factors (1,1). |
| `geometry_extrapolation` | Twill/glass stack; widths 0.015, 0.04, 0.06, 0.08, 0.15 m and ply counts 8, 16, 32, 100. |
| `material_extrapolation_v2` | Twill/glass stack; factors (2.25, 0.75); double the nominal resin-to-fibre volume ratio. |

For layup transfer, the half-stack repeats U-T-G and is truncated to the required half-ply count before mirroring. Carbon UD and twill plies are 0.20 mm; glass plies are 0.10 mm. The UD elastic card in the example configuration is a representative numerical assumption for this method experiment, not a sourced or measured project material specification. Replace it with the desired card before treating results as material predictions.

The 100-ply alternating twill/glass geometry has 50 plies of each type and 15 mm total thickness. Geometry cases are marked separately by whether width and ply count are inside the Dataset A ranges. Thus the extreme 0.015/0.15 m widths and 100-ply cases are explicit extrapolation points.

### Set 5 resin-ratio proxy

The example treats each supplied ply card as an effective lamina at a 50% baseline resin volume fraction. This is an illustrative baseline: 50% is the midpoint of the cited woven-carbon laminate fibre-volume range, and the same baseline is assumed for the glass card. `resin_to_fibre_ratio_multiplier=2` doubles the resin:fibre volume ratio, then renormalizes the two fractions; it changes the resin fraction from 0.50 to 0.667. The matrix proxy uses epoxy tensile modulus 3.0 GPa and Poisson ratio 0.35, giving shear modulus 1.11 GPa. The modulus is taken from a Westlake EPIKOTE resin technical data sheet; Poisson ratio is an explicit assumption.

For each elastic or shear modulus $p_0$ in the configured effective ply card, the code infers a notional fibre contribution $p_f=(p_0-V_{r,0}p_m)/(1-V_{r,0})$ and recombines it as $p=(1-V_r)p_f+V_rp_m$. This simple Voigt mixture is a transparent stiffness perturbation, not a full micromechanics model. Ply thickness, density and Poisson ratios stay fixed; therefore this set tests compliance transfer only and must not be used for mass optimization. The source material cards and the resin proxy are not validated physical properties.

In Set 5, the target FEA uses the 2× resin:fibre ratio and the configured carbon/glass factors. The model is deliberately not told that ratio changed: its nominal section-property features assume a ratio multiplier of 1, and calibration adjusts only the two carbon/glass factors using the two observed compliances. Set 5 therefore tests how the two-factor correction handles an additional unmodelled composition shift; it does not identify or recover the resin ratio.

Sources for illustrative resin assumptions: [Westlake Epoxy EPIKOTE MGS LR285 technical data sheet](https://shop.swiss-composite.ch/shop/resources/downloads/110_1051/t-Epoxyd-Harz-L-285-LF-e.pdf) (typical cured tensile modulus 3.0 GPa, density 1.19 g/cm$^3$); [Gurit SE84LV prepreg data sheet](https://www.revchem.com/document/15pc135966a_tds-pdf/) (woven carbon laminate fibre volume fraction 47–53%).

## FEA definition

The model is a 0.25 m straight cantilever, with a 10 N total transverse force distributed over the free tip edge, a fully fixed root edge, and SHELL181 layered shell elements on a structured 40 x 6 mesh. Three through-thickness integration points are assigned to each ply. MAPDL is set to large-deflection static analysis (`NLGEOM,ON`) with automatic time stepping, 10 initial substeps, 1 to 100 permitted substeps, up to 50 equilibrium iterations per substep, and full Newton-Raphson. This accounts for geometric nonlinearity from large rotations; it does not add material nonlinearity, damage, failure, contact, or dynamics.

The output is mean vertical tip displacement and compliance `abs(delta)/P`. Because the formulation is geometrically nonlinear, compliance can depend on load. The current experiment fixes the force at 10 N and does not claim predictions for another load level; a new load requires regenerating both datasets and refitting. Each record also includes complete ply sequence, count, total thickness, mass, `D_eff`, `S_eff`, `EI`, `kGA`, normalized displacement, mesh counts, and per-case MAPDL wall time `fea_runtime_s`.

This is a 2D shell discretization with a layered section, not a 3D solid analysis. It is a nonlinear static compliance experiment, not an impact, energy absorption, failure or physical validation study. Check mesh convergence, reactions, shell local axes, nonlinear convergence and load-step histories before interpreting the sweep. A 100-ply case takes substantially more section data than the shorter stacks and should first be tried as a single smoke-test case if MAPDL runtime or memory is limited.

## Models, calibration and validation

All three models are fit using Dataset A only:

1. **Mechanics-informed compliance regression:** nonnegative no-intercept least squares in the bending basis `L^3/(3 EI)` and shear basis `L/(k GA)`.
2. **Quadratic response surface:** second-order polynomial in standardized log width, log total thickness, log `D_eff` and log `S_eff`, fitted to log compliance with ridge regularization.
3. **Gaussian process:** Matérn 5/2 kernel with a white-noise term, using the same log predictors and log compliance.

The quadratic and GP models exponentiate their predicted log compliance so their outputs remain positive. Five-fold cross-validation is run on Dataset A. No target-set FEA response enters training or tuning.

For each of the five target sets, the script selects two distinct observations randomly and repeats the selection three times. The same case IDs are used for all three model types. Nonlinear least squares estimates two effective property multipliers (carbon and glass) from those two target responses. The remaining 18 cases in that set are held out for validation. Set 5's resin ratio is intentionally omitted from model features and is not fitted. The fitted multipliers are not guaranteed to be identifiable: each repeat reports the calibration Jacobian rank and condition number, and flags rank deficiency or a condition number over 1e6. Treat unstable multipliers as evidence about the design, not as recovered physical material constants.

Metrics are reported before calibration and after calibration, by target set and repeat. In the geometry-extrapolation set they are additionally separated by width and ply-count regions. Validation excludes both calibration cases. Dataset A cross-validation and target validation remain separate results.

## Requirements and commands

Install the repository's Python dependencies, PyMAPDL and a licensed local MAPDL installation. From PowerShell at the repository root, first inspect the deterministic design tables without starting MAPDL:

```powershell
python src/run_extended_transfer_study.py --stage design
```

Then run the FEA stages and model evaluation:

```powershell
python src/run_extended_transfer_study.py --stage dataset-a
python src/run_extended_transfer_study.py --stage targets
python src/run_extended_transfer_study.py --stage fit
```

Or run the entire sequence with `--stage all`. If a simulation run stops, the FEA stages read completed case IDs from their CSV and skip those cases when resumed. Keep the config and solver logs with the output. To change output location, pass `--output-dir results/my_run`.

The runner includes the effective bending/shear calculation it needs and does not import the legacy `laminate_section.py` module. The archived baseline `run_property_transfer_study.py` is not invoked. Smoke-test one Dataset A case and one Set 5 case before the full run:

```powershell
python src/run_extended_transfer_study.py --stage dataset-a --max-cases 1
python src/run_extended_transfer_study.py --stage targets --target-set material_extrapolation_v2 --max-cases 1
```

Then run `python src/run_extended_transfer_study.py --stage all`. The runner resumes completed case IDs from the same output folder, so it will keep the smoke-tested Set 5 case and run the remaining cases. Do not use `--max-cases` with `--stage all`, because model fitting requires the complete datasets.

### Accepted input domain and failure handling

The JSON configuration uses SI units. The supported design has a positive finite beam length and fixed positive tip load; five unique positive Dataset A widths; Dataset A ply counts exactly 8, 16, 24, and 32; and the four carbon/glass scale pairs `(1,1)`, `(1,2)`, `(2,1)`, `(2,2)`. The geometry target accepts five positive widths and even ply counts from 2 through 100. The current example uses 0.25 m length and 10 N. Material cards require positive finite ply thickness, density, elastic and shear moduli, carbon/glass family labels, and bounded finite Poisson ratios. The resin proxy requires a baseline fraction between zero and one, positive matrix modulus, and valid matrix Poisson ratio. Only the `tw_glass` and `ud_tw_glass` stack families and zero-degree plies are supported.

Large-deflection mode and automatic time stepping must be enabled. Calibration requires at least one repeat, at least two observations, and positive increasing scale bounds; target factors must be within those bounds. Invalid or missing required values fail before MAPDL is launched. Unsupported layup family, odd ply count, or nonzero orientation raises a clear error. Reusing a CSV from another FEA definition raises an error and requires a fresh output directory. In particular, the previous linear-static CSVs cannot be resumed as nonlinear results.

Outputs by default under `results/property_transfer_large_deflection_5sets/`:

- `dataset_A_design.csv` and `target_design.csv`: planned cases, scale factors and full bottom-to-top layup strings.
- `dataset_A.csv` and `targets.csv`: FEA responses and effective section/mass properties.
- `model/model_*.joblib`: fitted model objects.
- `model/calibration_splits.csv`: the three paired calibration selections per target set.
- `model/predictions.csv`: per-case nominal and corrected predictions for each model/repeat.
- `model/metrics_by_repeat.csv`: held-out scores and calibration diagnostics.
- `model/runtime_by_model.csv`: model fit time, five-fold CV time, and mean calibration/prediction time.
- `model/runtime_by_repeat.csv`: calibration/prediction wall time for each target set and repeat.
- `model/summary.json`: Dataset A cross-validation scores and run counts.

Runtime records are wall-clock measurements from the machine and MAPDL installation used for the run. Compare model timings on the same computer and process state. The model runtime output allows an accuracy-versus-cost comparison; per-case `fea_runtime_s` allows the FEA cost to be compared with surrogate fitting and prediction. Do not substitute timings from the historical linear run for the nonlinear study.

Expected size: 80 training + 100 target FEA cases, or 180 unique analyses. Each target-set calibration/validation cycle has two calibration cases and 18 held-out cases per repeat. A complete fit evaluates 5 sets x 3 repeats x 3 model families.

## Interpretation guardrails

- Dataset A and all targets share the same FE modeling framework; success is FEA-to-FEA transfer only.
- The mechanics-informed bending/shear basis is derived from small-deflection beam theory; with NLGEOM targets it is a comparison baseline, not an exact nonlinear mechanics solution.
- Compliance is load-dependent in a geometrically nonlinear analysis. The model's domain is the single configured load used for both datasets; load variation is unsupported without a new sweep and model fit.
- Two observations are the algebraic minimum for two fitted property multipliers. A poor Jacobian condition means the pair cannot robustly distinguish their effects.
- Material interpolation/extrapolation labels refer to the configured stiffness factors. Geometry extrapolation labels are based on width and ply count separately.
- The new UD card and the other representative material cards are computational assumptions, not validated material data.
- A good result does not validate physical landing gear response or curved-gear performance. Physical observations, model verification, and the curved structural geometry remain future work.
