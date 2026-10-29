# Cubic Hermite: transverse warping completeness and resolved shear benchmark

**2026-10-27 — Research evidence, not 3-D a posteriori certification.**

## Root cause identified

The cubic recovery stencil builds total-degree-three displacement Taylor terms only where there are enough distinct nodal layers along an axis. For `ny=3`, the term `(0,3,0)` is absent: there are only three distinct `y` coordinates from which to identify `y^3`. Consequently, the interpolation is *not cubic-complete across beam depth*, despite the name `cubic_rotation_hermite`. That is especially relevant to transverse shear of a bent beam: the displacement field may require a longitudinal **cubic warping in y**, whose derivative generates a **quadratic shear strain profile**.

The prescribed kinematic patch `u_z = y^3`, nodal physical spin `omega_x = 1.5*y^2` has exact engineering shear `gamma_yz = 3*y^2` and zero other small-strain components. Evaluating the *existing source B and T*, including interior and exact-corner points:

| Node mesh | `(0,3,0)` recovery term | Relative sampled 6-component strain L2 |
|---|---|---:|
| 3×3×9 | Absent | ~33% |
| 3×4×9 | Present | approximately numerical roundoff |

This is a **polynomial reproduction test**. It does not impose actual shear equilibrium, model free-surface traction or prove accuracy in a solved elasticity problem. Source-element stiffness, material, load integration and constraints were not changed.

## Actual solved cantilever refinement (original unmodified FE stresses)

All rows use the same prescribed Redwood-like orthotropic material (longitudinal E = 1,100,000 psi, without calibration), end-face integrated downward force (-12.8 lbf), gravitational weight density (0.0174 lbf/in³), KKT-clamped translations and rotations, 4-point Gaussian sectional quadrature, parent axial section positions `-0.9,+0.9`, and **exact source-element-corner stress sampling**. Stress comparison uses beam-theory reference as a screen, not a 3-D exact solution.

| Nodal mesh nx×ny×nz | Tip displacement difference vs beam | Interior nodal sigma_zz relative L2 vs beam | Interior nodal tau_yz relative L2 vs beam | Max interior sectional shear difference vs beam | Max interior sectional moment difference vs beam |
|---|---:|---:|---:|---:|---:|
| 3×3×9 | 0.465% | 5.451% | **73.875%** | 16.670% | 10.430% |
| 3×4×9 | 0.475% | 5.141% | 27.759% | 12.256% | 10.130% |
| 3×5×9 | 0.076% | 4.532% | 29.768% | 12.940% | 9.727% |
| 3×4×17 | 0.389% | 1.645% | 10.065% | 1.952% | 3.588% |
| 3×5×17 | 0.670% | 1.465% | 9.473% | 1.811% | 2.999% |
| **3×4×25** | **0.842%** | **0.603%** | **3.570%** | **0.302%** | 1.212% |
| **3×5×25** | **0.989%** | **0.398%** | **1.669%** | **0.604%** | **0.998%** |
| 3×5×33 | 1.189% | 0.123% | 1.074% | 0.117% | 0.305% |
| **3×5×41** | 1.321% | 0.094% | **0.918%** | 0.051% | 0.026% |
| **3×5×49** | 1.413% | 0.142% | **0.912%** | 0.022% | 0.189% |

The **1.669%** value is `100 * ||tau_FE - tau_beam||_2 / ||tau_beam||_2`
over equally weighted, element-sided source-corner samples with
`L/4 <= z <= 3L/4`. Shared vertices contribute one value per owning element;
no cross-element averaging is applied. This measures a discrete comparison
with classical beam theory, rather than a maximum pointwise or
volume-integrated 3D error. It does not establish stress equilibrium.

The same axial-refinement series at fixed 3×5 cross-section produced 1.074%
at 33 axial nodes, 0.918% at 41, and 0.912% at 49. So the reported corner-sample
discrepancy does fall below 1% by 41 axial nodes. The tip comparison moves the
other way over these runs: it grows from 0.989% at 25 nodes to 1.321% at 41
and 1.413% at 49. This is measured behavior for the stated discrete samples,
not proof that the true 3D shear error converges below 1%.

For direction, define the fitted amplitude `alpha = sum(tau_FE*tau_beam) /
sum(tau_beam^2)` over the same equally weighted samples. It is 1.0066 for
3×5×25 and 1.0028 for 3×5×41, so the FE samples have a slightly larger
overall fitted magnitude than the beam reference on those meshes. The 3×5×9
and 3×5×17 values are 0.9172 and 0.9860, respectively: the aggregate direction
changes as the mesh is refined. Individual FE values exceed and fall below the
reference on the same mesh, so “larger” does not hold pointwise everywhere.

**We can now demonstrate less than 5% shear disagreement at the *actual original source vertices* using raw FE constitutive stress** for both 3×4×25 and 3×5×25. This is not achieved by replacing FE shear with a theoretical parabola, by moving sampling points to centroids or by calibrating E. No new numerical method was applied to the stress field. The important change is **mesh resolution that supports the required cubic transverse displacement mode and resolves axial variations**.

The FE shear at nominally free top/bottom surfaces is still not exactly zero. In the quarter-to-three-quarter-span source-corner samples, the maximum absolute `tau_yz` at `y=±4` dropped from ~1.874 psi (3×3×9) to ~0.110 psi (3×5×25). **This does not establish full surface-traction admissibility.** Full face stresses, other traction components and residuals still need scrutiny.

## Additional local-balance warning

An extended 3×5×25 sectional audit (4 axial Gaussian samples per source cell) reports
**0.586% worst sampled interior shear-resultant difference** versus the
classical reference, but **31.360% worst sampled shear-vs-interpolated-FE-`dM/dz`
residual**, normalized by theoretical shear. The derivative uses a local
polynomial fitted to FE moments at Gaussian points; this is an approximation,
not an exact derivative or rigorous 3-D equilibrium residual. It nevertheless
shows why excellent one-dimensional beam-stress agreement does **not** establish
full local equilibrium. The preliminary full volume-residual sampling did not
complete within this environment's limit. No claims are made about it.

`outputs/shear_3x5x25_axial_balance_research.json` preserves this diagnostic,
including per-cell samples, rather than concealing it behind the favorable
source-vertex L2 percentages.

## Reproduction on Windows PowerShell

```powershell
python -m unittest tests.test_cubic_shear_warping -v
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --nodal-stress-audit --shear-warping-audit --out outputs/shear_coarse_3x3x9.json
python -m tutorials.cubic --nx 3 --ny 4 --nz 25 --section-only --nodal-stress-audit --shear-warping-audit --out outputs/shear_improved_3x4x25.json
python -m tutorials.cubic --nx 3 --ny 5 --nz 25 --section-only --nodal-stress-audit --shear-warping-audit --out outputs/shear_improved_3x5x25.json
python -m tutorials.cubic --nx 3 --ny 5 --nz 33 --section-only --nodal-stress-audit --out outputs/shear_depth_3x5x33.json
python -m tutorials.cubic --nx 3 --ny 5 --nz 41 --section-only --nodal-stress-audit --out outputs/shear_depth_3x5x41.json
python -m tutorials.cubic --nx 3 --ny 5 --nz 49 --section-only --nodal-stress-audit --out outputs/shear_depth_3x5x49.json
```

The `--shear-warping-audit` switch separately writes `*_shear_warping_patch.json`. A coarse cubic mesh with only three `y` node layers now prints an explicit warning and records its missing transverse cubic term in the main run JSON.

## Important limitations and remaining steps

1. Percentages are **comparisons with classical beam theory, not certified upper bounds** on 3D elastic stress error. The orthotropic solid may differ from the one-dimensional idealization.
2. The nodal L2 is a *discrete, element-sided, unweighted stress-sampling norm* over quarter to three-quarter span. Meshes have different sets/counts of corners; therefore the table is evidence of improved agreement, **not a strict mesh-convergence study on an identical integration norm**.
3. The 3×5 source-corner shear metric decreases across the listed 9-to-49-node axial series, but do not treat this as proof of general solution convergence: the equally weighted samples change with the mesh, other quantities (including tip displacement) do not improve monotonically, and a common spatial-integration norm plus a tighter 3D reference are still needed.
4. Neither the 3×4×25 nor 3×5×25 example has been shown to provide a globally equilibrated stress tensor satisfying all natural tractions, nor a mathematical a posteriori error bound.
5. Preserve the raw FE stresses. A separate smooth shear profile reconstructed from FE sectional moment gradients tested favorably, but it was **not adopted as a substitute for FE shear**, because it assumes beam-like through-thickness shape and could conceal the real stresses.
6. An experimental least-squares penalty on free-surface `tau_yz` was also tested and **rejected**: it could reduce a sectional discrepancy for some penalties but did not resolve the nodal shear and degraded other quantities. The production stiffness remains unchanged.
7. Next: source-Gauss volume `div(sigma)+b` residual, face tractions at all free sides, common-volume-norm stress refinement, and eventually equilibrated three-dimensional stress reconstruction for certification.
