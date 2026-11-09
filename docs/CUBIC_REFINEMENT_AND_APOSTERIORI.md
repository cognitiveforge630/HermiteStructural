# Cubic refinement and a posteriori verification — research plan

**Status (2026-10-26): a posteriori preparation exists; no mathematically certified Hermite error bound has yet been demonstrated.** The current fail-closed checklist is `docs/CERTIFICATION_STATUS.json`. `python -m tests.checks.certification_readiness --require-certified` intentionally exits with an error. Keep it that way until independently verified proofs and numerical checks satisfy the listed obligations. The correct technical term is *a posteriori* error estimation.

## New measured section-QoI refinement study — 2026-10-26

Direct Gauss section integration now produces resultants from the *original*
`D @ B @ U` strain/stress field without using visualization subdivisions.
It is available through `core.solvers.section_equilibrium.integrate_sections`,
and runs automatically for all cubic examples, even with `--section-only`.
The tables in this section report the nonlocal cubic Hermite formulation on
Hex8 connectivity; conventional Hex27/Hex64 data are reported separately in
the [higher-node study](HIGHER_NODE_HEX27_HEX64_STUDY_20261027.md).

| Node mesh | Tip vs AWC (%) | Worst midspan shear vs beam (%) | Worst midspan moment vs beam (%) |
|:--|--:|--:|--:|
| 3×3×5 | 1.1646 | 18.8384 | 9.3508 |
| 3×3×9 | 0.4651 | 20.5657 | 3.2255 |
| 3×3×17 | 0.4079 | 8.0528 | 2.2964 |
| 3×3×33 | 1.1081 | 6.1162 | 0.3158 |
| 5×5×9 | 0.1864 | 15.0603 | 2.5935 |

**Non-monotonic warning:** 3×3×33 has a larger beam-tip discrepancy than
3×3×17 despite smaller sectional shear/moment discrepancies. This is a
credible obstacle, not an improvement to hide. The original mesh family
is only 3 nodes across x and y, and purely axial refinement is not a uniform
h-convergence experiment. The true 3D elasticity solution need not equal the
approximate beam-theory formula either. Integrated sectional resultants are
additional research QoIs, not a bound or verified admissible stress.

Subsequent section-only runs completed the 5×5×17, 5×5×33 and 5×5×65
series. Their worst interior midpoint shear-resultant discrepancies were
2.0419%, 0.3905% and 0.0045%, respectively. Tip disagreement increased along
this series. See the [higher-density study](DENSIFICATION_STUDY_20261026.md)
for the extended mesh test, section-location sensitivity and saved artifacts.
Inspect solver/KKT constraints and the source-Gauss shear resultants before
paying for visual resampling.

```powershell
python -m tutorials.cubic --nx 5 --ny 5 --nz 17 --section-only --out outputs/cubic_5x5x17.json
python -m tools.section_refinement outputs/cubic_3x3x9.json outputs/cubic_5x5x17.json
```

The tool above needs files written with the **new** runner, including paired
`*_source_sections.json` reports; the archived `samples/cubic_3x3x9.json`
contains only a visualization audit and is not interchangeable. For independent
reproduction, rerun the 3×3×9 case using `--section-only` first.

**Boundary-condition sensitivity:** on 3×3×17, changing from the published
fixed-nodal-rotation clamp to a displacement-only clamp produced tip-reference
disagreement 1.4192% (versus 0.4079%) and midspan section-shear disagreement
7.0156% (versus 8.0528%). These are *different boundary problems*. A
classical 3D Cauchy displacement-Dirichlet specification does not intrinsically
prescribe six vertex DOFs; a genuine model review must settle which rotational
conditions are physically appropriate rather than choose whichever improves
one numerical QoI. For a diagnostic use `--displacement-only-clamp` and keep
its results separate from the published fixed-spin series.

## Evidence already obtained

| Recovery | Nodes | Tip `uy` (in) | Tip vs analytical beam | Interior `σzz` relative L² vs beam | Interior maximum shear-resultant discrepancy |
|---|---|---:|---:|---:|---:|
| Cubic | 3×3×5 | −0.114796439 | 1.165% | 9.12% (`--subdivisions 2`) | 28.40% (`--subdivisions 2`) |
| Cubic | 3×3×9 | −0.115608854 | 0.465% | 2.16% (`--subdivisions 3`) | 19.82% (`--subdivisions 3`) |

Beam-theory `uy = −0.116149091 in`; load and Redwood reference E are identical. These are *reference-model discrepancies*, **not true-solution error estimates**. Because 3×3×5 used 2 subcells per element axis and 3×3×9 used 3, the stress comparison mixes discretization change with stress-sampling change. The actual 3×3×9 stress metric changes from ~1.16% at subdivision 2 to ~2.16% at subdivision 3 **without changing the FE solution**. Never claim these data prove monotonic stress convergence, global accuracy, or a certified error below 10%.

## Controlled next experiments

Run **additional** solves on your own machine, *one at a time*. The new direct-section mesh runs above were completed, but their full visualization stress audits were not recomputed. Large meshes are **not** in the automated test suite. Stencil nonlocality, global sparse solves, and a dense KKT constraint-rank test can become expensive; do not presume cubic mesh refinement is cheap, stable, or monotonically better.

First, redo the existing 3×3×9 case with fixed sampling and a distinct output filename, to create fresh comparable metadata:

```powershell
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --subdivisions 3 --out outputs/cubic_3x3x9.json
```

Then increase **only the beam-direction resolution** (this is an axial refinement study, not a uniform 3D h-study):

```powershell
python -m tutorials.cubic --nx 3 --ny 3 --nz 17 --subdivisions 3 --out outputs/cubic_3x3x17.json
python -m tools.cubic_refinement outputs/cubic_3x3x9.json outputs/cubic_3x3x17.json --json-out outputs/cubic_axial_comparison.json
```

**Optional, more computationally demanding**, refine all axes relative to the 3×3×9 case:

```powershell
python -m tutorials.cubic --nx 5 --ny 5 --nz 17 --subdivisions 3 --out outputs/cubic_5x5x17.json
python -m tools.cubic_refinement outputs/cubic_3x3x9.json outputs/cubic_5x5x17.json --json-out outputs/cubic_uniform_h_comparison.json
```

3×3×9 → 5×5×17 halves *all three element edge lengths*. By contrast 3×3×9 → 3×3×17 halves only the axial element length; it does not constitute a uniform spatial refinement. Do not run 9×9×33 or larger until the recovery/KKT assembly, memory budget, solve residuals, and sampling costs are independently profiled. Keep `E_L=1,100,000 psi`, `P=−12.8 lbf`, gravity, clamp conditions, quadrature, `--subdivisions 3`, and stress extraction definitions unchanged.

The offline comparator `python -m tools.cubic_refinement` **does not solve** and **does not certify**: it validates material/load/boundary consistency, displays each beam-reference discrepancy, measures the change in tip deflection and in *sampled* peak longitudinal tension between meshes, and flags mismatched sampling or non-uniform refinement. It reads the `*_stress_audit.json` files beside their paired result JSONs. Older shipped audit files do not contain `sample_subdivisions`; the tool infers it from section counts. For a fair stress claim, rerun both grids using the same sampling and check sampling sensitivity separately (`--subdivisions 3` versus 4 on a fixed solution when supported).

### Numerical acceptance targets (research only)

Aim for **≤0.5% tip difference** versus independent analytical beam theory (already met by 3×3×9), **≤1% interior longitudinal-stress beam comparison**, and **≤5% changes** in carefully defined scalar QoIs between *successive refined* meshes, with a separate goal of reducing section shear discrepancies. These are *targets*, not results or proof of a true-solution error. Record tip `uy`, integrated section moment and shear, signed bending stress, any support-excluded stress QoIs, matrix/constraint residuals and timings. Do not use clamp-corner raw stress maxima as a convergence certificate: they can be sensitive to local singularities and sampling.

## What the a posteriori certification pipeline still needs

1. **Define the exact variational problem**: 3D orthotropic Cauchy elasticity versus the recovered-rotation formulation, energy norm, load functionals, essential/natural faces, and how nodal rotations relate to displacement. Prove a stable, conforming approximation or derive a nonconforming estimator with the required jump/consistency terms. Patch tests alone do not provide this.
2. **Admissible displacement field**: verify the actual imposed fixed-face trace continuously (not merely at a finite set of KKT samples), interior continuity, and the conditions needed for the bound.
3. **Statically admissible stress reconstruction**: construct an equilibrated stress field satisfying body force, external traction, internal-face traction continuity and any required interface conditions. The present beam stress reference and sampled `D B U` stress are **not** an equilibrated reconstruction.
4. **Bound derivation and numerical enclosure**: for linear elastic problems, a complementary-energy/hypercircle-type construction may bound the error *after* prerequisites are proved; quantify load approximation, quadrature, algebraic solve residual, constitutive tensor uncertainty and rounding where needed. Never substitute an asymptotic mesh-difference or comparison with beam theory for a guaranteed bound.
5. **Independent acceptance**: manufactured-solution tests with analytic errors; effectivity on at least one certified example; energy-error and chosen stress/displacement QoI bounds (dual-weighted estimates require their own justification); rigorous automated fail-closed status before any `certified_bound=true` claim.

Current project support: element-mode tests, load-resultant checks, finite-element stress sampling, section checks, automated unit tests, source-identity checks, and a fail-closed certification readiness JSON. **The actual equilibrated-stress estimator, proved reliability constant and guaranteed bounds are not yet implemented.** That is precisely the difference between being prepared to *develop* a posteriori certification and having achieved it.
