# Stress recovery at original Hermite source vertices (27 October 2026)

**Classification: RESEARCH / NOT A POSTERIORI CERTIFIED.** This change
alters *diagnostic sampling*, not element stiffness, recovered interpolation,
nonlocal stencils, applied loads, or boundary conditions.

## Matrix construction and singularity checks

- The source displacement/strain matrix `primitive_B(xi, eta, zeta, coords)`
  maps recovered source DOFs to six engineering strains in order `xx, yy,
  zz, xy, yz, xz`. Its geometrical Jacobian is obtained from the Hex8
  coordinate mapping and checked for positive determinant.
- Cubic Hermite nodal derivatives are reconstructed by a polynomial patch fit
  across neighboring source vertices. The patch-rank test prevents an
  underdetermined reconstruction. This is different from **evaluating** `B`
  at a coordinate; recovery is not re-fitted for every evaluation point.
- Element stiffness, consistent volume and surface loads, and FE energy
  continue to use interior Gaussian integration points, where suitable.
  Replacing quadrature by eight unweighted corners would not provide a
  correct general stiffness/force integral and is deliberately avoided.
- A 3×3×5 cubic source-mesh corner test verified finite B arrays and positive
  Jacobian determinants at all eight exact source corners in multiple cells.
  This does **not** prove the absence of degenerate Jacobians or high
  conditioning on arbitrary distorted meshes. Exact corner B is allowed
  only while geometric mapping remains regular.

## What changed

1. `core/verification/nodal_stress.py`: raw stress is computed at each
   original *source Hex8* vertex with `sigma = D @ (B @ T @ U)`. Reports every
   source-element-side stress and a separate analytical beam comparison;
   no cross-element averaging or Gauss-point extrapolation occurs. Optional
   inset evaluation `(-f,+f)^3` with `0<f<1` is available.
2. `core/visualization/sampling.py`: stress is computed at eight corners
   of each visualization subcell, not at the subcell center. Physical corner
   coordinates and six-component raw stress values are retained in memory.
   The displayed cell value is explicitly the eight-corner arithmetic mean;
   **do not mistake it for a nodal stress or a maximum**.
3. `core/visualization/cubic_stress.py`: the independent analytical reference
   is sampled at exactly the same physical corner coordinates; L2 and
   utilization/peak calculations use the individual source-element-sided
   subcell-corner stresses, NOT the cell means. Historical centroid metrics
   have a different sampling norm and cannot be compared directly.
4. The tutorial default cross-section positions are now inset toward each
   element end (`zeta=-0.9,+0.9`) instead of `zeta=0` (which remains opt-in
   for historical reproduction). Cross-sectional force/moment integrations
   **still require Gauss quadrature** in x/y.
5. For volume residual sampling, the default is the eight-point interior
   Gaussian rule (`--volume-residual-order 2`) rather than the single
   mid-element diagnostic.

## Measured 3×3×9 case (same solver, no recalibration)

| Diagnostic | Result |
|---|---:|
| Tip displacement difference vs AWC | 0.465124% |
| Original source **exact-corner** longitudinal σzz discrete interior L2 | 5.451% |
| Original source **exact-corner** shear τyz discrete interior L2 | 73.875% |
| Original source **95%-inset-corner** longitudinal σzz discrete interior L2 | 5.247% |
| Original source **95%-inset-corner** shear τyz discrete interior L2 | 68.865% |
| Display-subcell raw exact-corner longitudinal σzz interior L2 (`--subdivisions 2`) | 4.776% |
| Near-face (`zeta=±0.9`) interior sectional shear discrepancy | 16.670% |
| Near-face (`zeta=±0.9`) interior sectional bending-moment discrepancy | 10.430% |

**Do not call these values a numerical convergence trend.** The source
vertex norm weights source-element-sided corners equally (repeating mesh
vertices), rather than approximating a certified 3D energy norm. Stress can
vary by element side, be sensitive to clamp and traction boundary effects,
and behave poorly at surfaces even when midpoint errors are small. The
approximate rectangular-beam shear reference is not exact orthotropic 3D
elasticity. The large shear numbers warrant independent analysis; switching
where the stress is evaluated is a disclosure/diagnostic improvement,
**not a remedy for the underlying shear imbalance**.

## Commands

```powershell
python -m unittest tests.test_corner_stress_audit -v
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --nodal-stress-audit --out outputs/original_source_vertex.json
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --nodal-stress-audit --stress-sample-location inset_corners --stress-inset-fraction 0.95 --out outputs/inset_source_vertex.json
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --stress-audit --equilibrium-diagnostics --out outputs/cubic_source_corner_full.json
```

Use `--section-axial-positions 0` only when intentionally repeating an
**older midpoint benchmark**. Do not combine figures from sampling variants
without labeling and checking matched positions and norms.

## Still missing for certification

The recovered stress does not yet satisfy proven strong-form volume equilibrium,
interface traction equilibrium on general geometries, or natural-boundary
admissibility to verified tolerances. No fully computable guaranteed energy
or QoI error bound has been established. Retain `certified_error_bound=false`.
