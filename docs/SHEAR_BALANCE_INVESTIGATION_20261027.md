# Cubic Hermite sectional-shear investigation — 2026-10-27

**Status: RESEARCH DIAGNOSTIC — NOT CERTIFIED. No solver correction is claimed.**

## What was changed and what was not

The `HermiteStructural(4).zip` and `(5).zip` Python-source comparison found
**no changes** in `core/elements/{base,cubic,quadratic,linear}.py`, materials,
loads, constraints, or `core/solvers/section_equilibrium.py`. The changed
`tutorials/cantilever.py` adds optional checks and output. There is therefore
no evidence that the new Mechanics of Materials references altered the
underlying solver, although they exposed existing inaccuracies.

This update **does not edit the stiffness, displacement interpolation,
constitutive matrix, stress operator, applied loads, or constraints**. It adds
an optional raw-stress diagnostic `core/solvers/shear_balance.py`, CLI flag
`--shear-balance`, four regression tests, and reports.

## Independently reproduced 3×3×9 baseline

All comparisons below retain the same cubic formulation, 4×8×120-inch Redwood
geometry, signed end traction -12.8 lbf and gravitational weight density
0.0174 lbf/in³. The default clamp fixes nodal translations and rotations.

- Tip displacement: **0.465124%** from the AWC beam reference.
- Worst interior **midpoint-only** section shear resultant: **20.565722%**.
- Worst interior midpoint section bending moment: **3.226%**.
- Sampled interior longitudinal stress relative L2 from user's run: **1.163%**.

All of these are beam-theory comparisons, not 3D certified stress errors.

## Why a single sectional measurement is insufficient

For each source axial cell, `--shear-balance` computes four *separate* Gauss
cross sections, and integrates each raw FE stress tensor over the cross
section. The FE stress is always `D @ (B @ U)`.

It reports three **different** quantities, without replacing one by another:

1. `raw_FE_shear_y_lbf = integral_A tau_yz dA` at individual axial points;
2. the **axial Gauss average** of that sectional shear within a cell; and
3. `dM_FE/dz`, calculated by differentiating the polynomial interpolant of
   the independently integrated `M_FE = integral_A y sigma_zz dA` at the
   axial Gauss points. This differentiation is a diagnostic approximation;
   it is not an exact pointwise derivative or an equilibrated stress.

For the default 3×3×9 solved case, the code reports, **within fully contained
middle-half cells**, the following extrema:

| Diagnostic | 3×3×9 default clamp |
|---|---:|
| Max sampled local-section shear disagreement with beam force | **18.814%** |
| Max disagreement in axial cell-*averaged* sectional shear | **8.986%** |
| Max mismatch of FE shear vs interpolated derivative of **FE** moment, scaled by beam shear | **118.216%** |

**Do not interpret the 8.986% axial average as a repaired stress error.**
Averaging smooths oscillations and can erase important local errors. The
peak values are *sampled* extrema, not known continuous maxima. On the same
mesh the existing element-midpoint-only metric remains **20.566%** because
that samples different axial coordinates. Higher-order axial differentiation
checks (3, 4, 5, 6, 8 samples) indicated substantial FE shear/moment mismatch
throughout; individual percentages depend on the diagnostic interpolation.

### Clamp experiment (not a solver fix)

Leaving nodal rotations free while imposing displacement constraints changed
3×3×9 midpoint shear discrepancy from **20.57% → 12.66%**, but tip
error worsened from **0.465% → 0.905%**. The new four-point axial diagnostic
still reported **20.79%** local sampled shear discrepancy, **8.03%** axial
cell-average discrepancy, and **105.77%** FE shear–moment balance mismatch.
Thus, merely relaxing clamp rotations is not an acceptable universal fix.

### Context from the earlier refinement study

The pre-existing 5×5×65 stored **three-slice** source-section report has
worst quarter-to-three-quarter-span sampled shear disagreement of **0.02296%**,
compared with **0.00448%** for the midpoint alone. This suggests that some
high-refinement interior section resultants really improve across multiple
axial positions. But it is **not** a guarantee of local stress accuracy or a
rigorous convergence rate. The 5×5×129 midpoint-only **0.0029%** result must
not be generalized to the full continuous stress field. Also, the finer
5×5×129 tip displacement disagreement remained **1.774%**.

A new attempted 5×5×17 solve in the current constrained environment was
terminated before completion; **no new 5×5×17 result is claimed**. Run larger
meshes on the user's workstation.

## Mechanics behind the remaining issue

For a traction-free lateral surface and static classical Cauchy equilibrium,
local balance is `div(sigma)+b=0`. Integrating the appropriate moment balance
across the section gives the expected relation `dM_x/dz = V_y`, with the
signs used in this repository, when lateral terms vanish. The FE solution
satisfies a *weak/discrete* equilibrium equation, but unmodified FE
constitutive stress need not satisfy the strong form pointwise.

The data implicate coarse-mesh oscillations, recovered-derivative / stress
accuracy, and/or boundary/lateral traction defects. The diagnostic does **not**
prove which mechanism dominates. The simple beam reference itself is
approximate for 3D orthotropic Redwood, especially near load and support.

## Run the diagnostic

```powershell
python -m tutorials.cubic `
    --nx 3 --ny 3 --nz 9 `
    --section-only --shear-balance `
    --out outputs/cubic_3x3x9_shear_research.json
```

Compare the original source-section JSON to the new
`*_shear_balance.json`. The output contains *all* sample coordinates,
FE shear, beam reference, FE moment, beam moment, derivative-of-FE-moment,
local discrepancies, and axial averages.

```powershell
python -m unittest tests.test_shear_balance -v
python -m unittest discover -s tests -v
python -m tests.checks.source_provenance
```

## Proposed controlled next engineering steps

1. Run the same **four-Gauss-point axial shear audit** at 3×3×9,
   5×5×17, 5×5×33, and 5×5×65. Compare interior *sampled local* metrics,
   cell-average metrics, and FE shear–moment balance separately.
2. Investigate on the coarse mesh whether the large FE `V_y - dM_x/dz`
   mismatch correlates with raw `div(sigma)+b` and top/side traction
   residuals; do not hide lateral traction effects.
3. Examine recovery stencils and nodal rotations near the clamp **without**
   changing material properties or imposing calibration. Add a separately
   named optional recovery variant only after exact patch and solved-load tests.
4. Compare with an independent conforming high-order solid solution with
   *identical* orthotropic constants, loads, and physical boundary conditions.
5. For certification, pursue admissible H(div) equilibrated stress with proper
   body-force and boundary-face constraints, and rigorous goal-oriented
   bounds. Raw stress smoothing and beam comparisons do not certify the solve.

No claim of a newly repaired local shear field or general a posteriori
certificate is made by this package.
