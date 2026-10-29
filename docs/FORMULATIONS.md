# What changes with linear, quadratic and cubic Hermite recovery?

**Important:** The terms refer to *recovery polynomial order* rather than the polynomial degree of a standard displacement-only Cauchy hexahedron. All three methods use the **same eight mesh vertices** and **six DOFs per vertex**: `u_x, u_y, u_z, theta_x, theta_y, theta_z`. The rotations represent physical vertex spin; the symmetric part of the displacement gradient is derived from surrounding translations.

| Feature | Linear | Quadratic | Cubic |
|---|---|---|---|
| Implementation | `core/elements/linear.py` | `core/elements/quadratic.py` | `core/elements/cubic.py` |
| Recovery patch | One-ring | Two-ring | Three-ring |
| Polynomial fitted to neighboring translations | First degree | Second degree | Up to third degree, subject to grid rank |
| Recovered quantities | First derivatives | First + selected mixed second derivatives | First + mixed second + xyz third derivative |
| Primitive unknown layout per Hex8 | 96 | 168 | 192 |
| Status | Very stiff in current bending benchmark | Historical 4×5×17 beam displacement check only (E=2m) | Experimental; measured Hex8 densification and stress studies, not certified |

The cubic results and densification plot refer to this **nonlocal Hermite
formulation on Hex8 connectivity**. The separate conventional C0 Lagrange
families now have both 3×3×5 and matched 3×5×25 results. At 3×5×25, Hex64
has smaller section-shear and section-moment discrepancies than Hex27, while
Hex27 is slightly closer on tip displacement. See the
[Hex27/Hex64 comparison](HIGHER_NODE_HEX27_HEX64_STUDY_20261027.md); neither
the comparison nor the beam reference establishes certified 3D accuracy.

## Shape functions and why the degree matters

The Hermite basis in `core/elements/base.py` is not replaced by a conventional Hex27 or Hex64 Lagrange shape family. Higher recovery orders **enrich the coefficients and mixed-derivative Hermite terms used to construct the displacement operator**. The resulting `primitive_matrices` and `element_transform` differ for the three formulations. A recoverable derivative at a vertex is computed from *global neighboring translations*, so these are nonlocal formulations, even though the mesh topology is Hex8.

In symbolic shorthand:

```text
Linear:    recovered 1st derivatives + physical vertex rotations
Quadratic: linear ingredients + xy, yz, xz mixed 2nd derivatives
Cubic:     quadratic ingredients + xyz mixed 3rd derivatives
```

Quadratic completeness depends on each two-ring recovery patch having full rank 9. The code has a linear-fit fallback for deficient quadratic patches. Cubic uses a rank-checked fit of the available total-degree-three polynomial terms; a small number of coordinate layers can omit unsupported terms. Accordingly, do **not** assume the names mean globally complete degree-2 or degree-3 finite-element approximation spaces on arbitrary unstructured meshes.

## What is constant across comparisons?

* Beam is 4 in × 8 in × 120 in; beam axis `z`, transverse force in global `-y`.
* Prescribed `E_z = 1,100,000 psi` for the new grade-reference benchmark, with identical research orthotropic ratios.
* End-face total force `P = -12.8 lbf`, integrated as traction with Gaussian quadrature.
* Specific weight `0.0174 lbf/in³`, integrated over the volume (total `-66.816 lbf`).
* Total applied y force `-79.616 lbf`.
* AWC-style analytical signed tip deflection `-0.116149 in` under the same point end-load and self-weight.
* Face displacement constraints via KKT with fixed face rotations, as implemented by the canonical runner.

A reference discrepancy is **not** the same thing as an a posteriori discretization-error estimator. Stress QoIs, mesh convergence, sensitivity to clamp modeling, and rigorous equilibrated stress bounds remain separate tasks.

## Revised Redwood source-grade baseline

As of 2026-10-26 the **new active benchmark** uses prescribed **E_L=1,100,000 psi**, not 2,000,000 psi, from Redwood Select Structural dimension-lumber material-reference data. End load −12.8 lbf and self-weight 0.0174 lbf/in³ are unchanged since they pass the preliminary grade/displacement screen. The previous quadratic 0.483% benchmark is historical at 2 million psi, not a fresh run of this material. See [Cubic stress research](CUBIC_STRESS_AND_REDWOOD.md).
