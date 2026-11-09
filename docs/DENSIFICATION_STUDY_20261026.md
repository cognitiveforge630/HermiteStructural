# Higher-density cubic densification study — 26 October 2026

Mesh densification greatly improved the measured interior **section-shear
resultant** agreement with the beam reference. The extended 5×5×129 test
reduces the midpoint discrepancy to **0.002908%**, compared with **20.565722%**
on the original 3×3×9 mesh: **99.9859% smaller**, approximately **7,072-fold**.
A separate three-slices-per-element test on 5×5×65 gives **0.022961%**.
These measurements concern integrated section forces and moments; they do not
provide a certified three-dimensional stress-field error percentage.

## Controlled axial densification

The 5×5 series holds transverse resolution fixed and doubles the axial
element count on each step. The same cubic recovered-rotation formulation,
Redwood material (`E_L=1,100,000 psi`), gravity (`0.0174 lbf/in³`), end force
(`−12.8 lbf`), displacement-trace clamp plus fixed nodal rotations, and SciPy
backend were used. No modulus or load was fitted to the reference.

Source-element stress is evaluated as `D @ B @ U` and integrated using
four-point Gaussian quadrature in both cross-section directions. The table
uses one section at each element's axial midpoint (`ζ=0`). Maxima are taken
over `30 ≤ z ≤ 90 in`, the middle half of the 120-inch span.

| Element formulation | Node mesh | Axial elements | Axial element length (in) | Tip vs beam (%) | Maximum midpoint shear-resultant discrepancy (%) | Maximum midpoint moment discrepancy (%) |
|:--|:--|--:|--:|--:|--:|--:|
| Cubic Hermite Hex8 (nonlocal recovery) | 5×5×9 | 8 | 15.0000 | 0.186354 | 15.060272 | 2.593506 |
| Cubic Hermite Hex8 (nonlocal recovery) | 5×5×17 | 16 | 7.5000 | 0.670870 | 2.041864 | 1.672473 |
| Cubic Hermite Hex8 (nonlocal recovery) | 5×5×33 | 32 | 3.7500 | 1.143716 | 0.390481 | 0.193469 |
| Cubic Hermite Hex8 (nonlocal recovery) | 5×5×65 | 64 | 1.8750 | 1.540615 | 0.004483 | 0.245697 |
| Cubic Hermite Hex8 (nonlocal recovery) | 5×5×129 | 128 | 0.9375 | 1.773562 | 0.002908 | 0.366085 |

![Measured 5×5 axial series and separate 7×7×17 transverse-refinement point](figures/cubic_densification.svg)

The open-diamond **7×7×17** observation is deliberately unconnected to the axial curve: it increases the transverse node count while returning to 17 axial nodes. Its recorded tip discrepancy is **1.062607%**; maximum midpoint section-shear discrepancy **2.079855%**; section-moment discrepancy **0.994642%**. Relative to **5×5×129**, displacement is *closer* to classical beam theory (1.063% versus 1.774%) while the measured section-shear and section-moment discrepancies are *larger* (2.080% versus 0.0029%, and 0.995% versus 0.366%). These are section integrals, not local stress norms; mesh and measurement locations change.

Both plot series are read directly from the saved solver JSON and corresponding source-section reports. Rebuild the SVG from the repository root without a new solve:

```powershell
python -m tools.plot_cubic_densification
```


The sampled midpoint shear discrepancy decreases on every step of the **5×5 axial series**, not over the separate 7×7×17 point. Doubling from 65 to 129 nodes reduces it a further **35.1%**.
Tip-to-beam disagreement increases, and moment discrepancy reaches its
smallest measured value at 33 nodes before increasing. Section locations
change with the mesh and the reduction factors vary substantially, so this
study does not establish a single asymptotic convergence rate or uniform
improvement across all quantities.

## Section-location sensitivity

The new `--section-axial-positions` option allows measurements away from
element midpoints. A fresh 5×5×65 solve used `ζ = −0.5, 0, +0.5`, corresponding
to the quarter, middle and three-quarter positions of each axial element.
This produced 192 sections, including 96 in the middle half of the span.

| Parent section position | Maximum interior shear-resultant discrepancy (%) |
|:--|--:|
| −0.5 | 0.021624 |
| 0.0 | 0.004483 |
| +0.5 | 0.022961 |
| Maximum over all three positions | **0.022961** |

The worst sampled section was at `z=89.53125 in`. The repeated solve's tip
and midpoint shear metric matched the saved 5×5×65 case. The expanded
section set still gives strong resultant agreement, while showing that
midpoint sampling produces a more favorable maximum. The 5×5×129 result
was evaluated at midpoints only; its off-midpoint maximum remains unmeasured.
No continuous-in-z maximum is enclosed by either finite section set.

## What has converged in these measurements

The positive result is a large reduction in the **sampled section-shear
resultant discrepancy**. Integrating a stress distribution can cancel local
errors, and the Galerkin equations enforce weak force balance. Accurate
resultants can therefore coexist with errors in the local stress distribution
or strong-form volume equilibrium. This study has not recomputed full-field
stress-map norms on the denser meshes, reconstructed an admissible equilibrated
stress, or supplied the exact 3D elasticity reference.

The 5×5×129 KKT system had 19,355 unknowns and 23,820,671 nonzeros. Its recorded
normwise relative backward error was `4.20e−17`; the raw infinity-norm residual
was `2.52e−3`, and the RHS-relative residual was `2.74e−3`. These diagnostics
use the system's equation scaling. A small backward error does not provide
a forward-error enclosure for this potentially ill-conditioned system, and
the section percentages have no solver/roundoff certification.

The next stress-accuracy work remains local-equilibrium verification,
equilibrated-stress reconstruction and review of the fixed-face rotational
constraints. The observed densification improvement supplies useful evidence
for that work while the certification gate remains false.

## Reproduce the higher-density tests

Run sequentially from the project root in Windows PowerShell:

```powershell
python -m tutorials.cubic --nx 5 --ny 5 --nz 129 --section-only --out outputs/cubic_5x5x129.json
python -m tutorials.cubic --nx 5 --ny 5 --nz 65 --section-only --section-axial-positions -0.5 0 0.5 --out outputs/cubic_5x5x65_three_slices.json
python -m tools.section_refinement samples/cubic_5x5x17.json samples/cubic_5x5x33.json samples/cubic_5x5x65.json samples/cubic_5x5x129.json
```

The section-only runs avoid visualization sampling. The comparator rejects
mixing the midpoint and three-slice cases because their section definitions
differ. All source result JSONs and paired `*_source_sections.json` files are
saved under `samples/`. The complete 5×5 midpoint series is summarized in
`samples/cubic_5x5_axial_comparison.json`.

Validation after adding the section-location option: **47 tests ran
successfully, with one offscreen graphics check skipped**. The tests check
invalid section positions and preservation of the previous midpoint metric.
The large solves are explicit research benchmarks outside the unit suite.
