# OOP Hex27/Hex64 commands (new independent project)

From `HermiteStructural_CubicHigherOrder_OOP` use standard Python, without
creating a virtual environment:

```powershell
python -m tools.demo_nodal_oop --order 3
python -m unittest tests.test_nodal_oop tests.test_p_order_study -v
python -m unittest discover -s tests -v
python -m tools.p_order_study --mode verify --order 3
python -m tools.p_order_study --mode plan --order 3 --nx 3 --ny 5 --nz 25
# Only the small 16-cell solved example, not the large 192-cell mesh:
python -m tools.p_order_study --mode solve --order 3 --nx 3 --ny 3 --nz 5
# Optional separate-process small-assembly RAM/time comparison:
python -m tools.benchmark_nodal_assembly --order 3
```

Open `HermiteStructural_OOP.code-workspace` in VS Code for the new OOP tasks.
See [the OOP tutorial](OOP_TUTORIAL.md) and [TODO checklist](OOP_TODO.md).
The original Hermite and viewer commands remain below.

---

# Reproducible commands (Windows PowerShell)

Run from the repository root. Use the normal globally installed Python interpreter, not a virtual environment.

```powershell
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python -m tests.checks.element_modes
python -m tests.checks.source_provenance
```

Then run the three *separate* demonstrations on the same modest grid:

```powershell
python -m tutorials.linear --nx 3 --ny 3 --nz 5
python -m tutorials.quadratic --nx 3 --ny 3 --nz 5
python -m tutorials.cubic --nx 3 --ny 3 --nz 5
python -m tools.compare_results
```

To run a **new, computationally expensive** higher-resolution quadratic case at the **current** modulus (the historical 0.483% result is stored separately):

```powershell
python -m tutorials.quadratic --nx 4 --ny 5 --nz 17
Get-Content .\outputs\quadratic.json
```

Convenience wrapper (defaults to quadratic; **check the selected formulation in the JSON**):

```powershell
python .\run_higher_order_hermite.py --nx 4 --ny 5 --nz 17
Get-Content .\hermite_result.json
```

WARNING: the following optional cubic full-grid research can be **very expensive** (this run was NOT completed in the current revision); run it deliberately rather than as part of a default unit test:

```powershell
python -m tutorials.cubic --nx 4 --ny 5 --nz 17
```

The runner CLI defaults to 4×5×17; **all published small-grid examples specify `--nx 3 --ny 3 --nz 5` explicitly**. The three example files save to different JSON paths so running one never overwrites another. `--out` explicitly overrides the destination. `--weight-density` changes gravity for *both* the FE loading and the analytical comparison. No `E` adjustment is performed.

## VS Code

Double-click `HermiteStructural.code-workspace`. Run **Terminal → Run Task** and choose the linear, quadratic or cubic task, or one of the three checks. For very large grids, consider using your appropriate SciPy/PARDISO installation; the published runner currently uses the portable SciPy backend.

## PyVista viewers (each original Hermite demo)

```powershell
python -m tutorials.linear --nx 3 --ny 3 --nz 5 --view --export-vtu
python -m tutorials.quadratic --nx 3 --ny 3 --nz 5 --view --export-vtu
python -m tutorials.cubic --nx 3 --ny 3 --nz 5 --view --export-vtu
```

Open `outputs/quadratic.vtu` after the fact with `python -m core.visualization.pyvista_viewer outputs/quadratic.vtu`, without re-solving. For full-screen, use `--full-screen`. For PNG, use `--screenshot`. See [VISUALIZATION.md](VISUALIZATION.md).

## Grade-aligned cubic stress study

```powershell
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --subdivisions 3 --export-vtu --view --full-screen --field "Sigma zz (psi)"
Get-Content outputs/cubic_stress_audit.json
python -m core.visualization.pyvista_viewer samples/cubic_3x3x9.vtu --field "Cubic minus beam sigma zz (psi)"
```

Default `E_L=1,100,000 psi` (published Redwood reference grade). Historical old `E_L=2,000,000 psi` numbers are *not* the new active examples.

## A posteriori readiness and cubic refinement (NOT yet certified)

```powershell
python -m tests.checks.certification_readiness
# Deliberate negative test; expected to exit nonzero until real certification:
python -m tests.checks.certification_readiness --require-certified
# One at a time, on your workstation. May be expensive:
python -m tutorials.cubic --nx 3 --ny 3 --nz 17 --subdivisions 3 --out outputs/cubic_3x3x17.json
python -m tools.cubic_refinement samples/cubic_3x3x9.json outputs/cubic_3x3x17.json --json-out outputs/cubic_axial_comparison.json
```

Existing 3×3×9 result already uses 3 sample subdivisions, matching the proposed
new study. The comparator reports beam-disagreement and scalar intermesh differences,
**not** a guaranteed a posteriori bound. See [the complete protocol](CUBIC_REFINEMENT_AND_APOSTERIORI.md).

## Higher-density section-resultant tests

```powershell
python -m tutorials.cubic --nx 5 --ny 5 --nz 129 --section-only --out outputs/cubic_5x5x129.json
python -m tutorials.cubic --nx 5 --ny 5 --nz 65 --section-only --section-axial-positions -0.5 0 0.5 --out outputs/cubic_5x5x65_three_slices.json
python -m tools.section_refinement samples/cubic_5x5x17.json samples/cubic_5x5x33.json samples/cubic_5x5x65.json samples/cubic_5x5x129.json
```

Run these solves sequentially. The [higher-density study](DENSIFICATION_STUDY_20261026.md)
documents the measured midpoint trend and the separate section-location check.

## Independent Mechanics of Materials checks (cubic)

```powershell
python -m unittest tests.test_mechanics_benchmarks -v
python -m tools.mechanics_benchmarks --out outputs/mechanics_benchmarks.json
python -m tutorials.cubic --nx 3 --ny 3 --nz 5 --section-only --mechanics-audit --out outputs/beer_validation_cubic_3x3x5.json
Get-Content outputs/beer_validation_cubic_3x3x5_mechanics_audit.json
```

The suite verifies classical 1D references, signed sigma/tau cube mechanics,
source-cubic polynomial patch modes and opt-in solved work/energy consistency.
It does not certify three-dimensional stresses or an a posteriori bound. See
[verification details](MECHANICS_OF_MATERIALS_VERIFICATION.md).

## Original-source nodal stress audit (new default sampling)

```powershell
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --nodal-stress-audit --out outputs/cubic_original_corners.json
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --nodal-stress-audit --stress-sample-location inset_corners --stress-inset-fraction 0.95 --out outputs/cubic_inset_corners.json
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --stress-audit --equilibrium-diagnostics --out outputs/cubic_corner_screen.json
```

The `_source_vertex_stress.json` files contain **all eight one-sided source
corner stresses per element**, separately from the beam reference. Current
`_stress_audit.json` uses raw subcell-corner samples, not legacy centroids.
Section-resultant default `zeta=-0.9,0.9` evaluates near element ends;
`--section-axial-positions 0` **explicitly re-enables historical midpoint
sections** for like-for-like comparisons with the published tables. Gaussian
quadrature is retained for cross-section/volume integration.
