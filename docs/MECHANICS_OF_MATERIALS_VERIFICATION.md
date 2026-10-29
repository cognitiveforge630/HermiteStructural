# Mechanics of Materials verification benchmarks (26 October 2026)

**Scope:** independent classical mechanics references and stress-tensor/element
regression tests for the active **cubic recovered-rotation Hermite** model.
This is a **verification programme for a research prototype**, not a material
failure model, design-code approval, a rigorous a posteriori error estimator,
or a certified three-dimensional finite-element solution.

Reference: Ferdinand P. Beer, E. Russell Johnston Jr., John T. DeWolf and
David F. Mazurek, *Mechanics of Materials*, 6th ed. (McGraw-Hill, 2012):
- Chapter 1, especially Sections 1.11--1.12, stress on inclined planes and
  the small stress cube (Fig. 1.34, **printed page 29**);
- Chapters 4--6, bending, shear and their stress distributions;
- Chapter 7, principal stresses and stress transformation;
- Chapter 9, deflection of beams;
- Chapter 11, elastic energy and work.

No book artwork or copyrighted text is reproduced in this project.

## Quick verification, Windows PowerShell (global Python)

```powershell
python -m pip install -r requirements.txt
python -m unittest tests.test_mechanics_benchmarks -v
python -m tools.mechanics_benchmarks --out outputs/mechanics_benchmarks.json
```

The **quick** commands construct a 3×3×5 cubic source mesh and perform tests
of independently prescribed polynomial displacement/stress fields. **They do
not run a structural load solve.** The CLI exits nonzero if any of its nine
analytical or polynomial-patch checks fails. Output JSON records the numerical
residuals, inputs, test locations, and `certified_error_bound: false`.

### Optional solved-model audit (more expensive)

```powershell
python -m tutorials.cubic --nx 3 --ny 3 --nz 5 --section-only --mechanics-audit --out outputs/beer_validation_cubic_3x3x5.json
```

This additionally writes:

- `outputs/beer_validation_cubic_3x3x5_mechanics_audit.json`: direct source-FE
  signed centerline displacements and all six stresses at three interior
  sections; **separate** textbook bending/shear references; principal stresses;
  traction on the +z cut; FE Gauss-integrated strain energy and load work;
- `outputs/beer_validation_cubic_3x3x5_source_sections.json`: already
  supported source-element section force and moment integrals; and
- `outputs/beer_validation_cubic_3x3x5.json`: existing runner metadata plus
  the new audit location and Clapeyron work/energy residual.

For further refinement, use `--nx 5 --ny 5 --nz 17` or other sizes **on your
own machine**; cubic nonlocal assembly and KKT solves may be expensive. The
new audit performs an extra source-element energy integration; it is opt-in,
including when `--section-only` is used. Existing default behavior is intact.

## Sign conventions and formulas

Global z runs along the beam from **fixed z=0 to free z=L**; x is width, y is
depth and +y is upward. P (lbf) and w (lbf/in) are **signed**, so downward
end traction and weight use negative values. The book's pictured cantilever
has its x-coordinate running in the opposite direction. All deflections here
are signed, not absolute magnitudes.

Define \(r=L-z\), \(I=b h^3/12\). For elementary Euler--Bernoulli bending:

\[
V_y(z)=P+w r,\qquad M_x(z)=-Pr-\frac12 w r^2.
\]
\[
u_y(z)=\frac{Pz^2(3L-z)}{6EI}
 +\frac{wz^2(6L^2-4Lz+z^2)}{24EI}.
\]
\[
\sigma_{zz}^{\rm beam}=\frac{M_x y}{I},\qquad
\tau_{yz}^{\rm beam}=\frac{3V_y}{2bh}\left(1-\frac{4y^2}{h^2}\right).
\]

The transverse shear formula is the elementary rectangular-section beam
approximation, **not an exact pointwise 3D elasticity solution**. In particular,
we do not equate a local 3D shear residual near the clamp to an FE error bound.

The independent beam bending energy and linear-elastic work relation are

\[
U_{\rm beam}=\frac{1}{2EI}\left(\frac{P^2L^3}{3}
 +\frac{PwL^4}{4}+\frac{w^2L^5}{20}\right),\quad
P u_y(L)+\int_0^Lw u_y(z)dz=2 U_{\rm beam}.
\]

All beam formulas exclude shear deformation and localized 3D effects.
The code tests point and self-weight effects individually before superposition.
The independent 1D energy identity uses five-point Gauss integration, exact
for the polynomial integrand, with no fitted elastic modulus.

### The Cauchy stress cube

Three normal and three independent shear components correspond to the book's
\(\sigma\), \(\tau\) notation. **Voigt order in the actual source is**
`xx, yy, zz, xy, yz, xz`. Stress uses physical \(\tau_{ij}\), whereas strains
use *engineering shear* \(\gamma_{ij}=2\epsilon_{ij}\). Consequently,
\(\tfrac12\boldsymbol\sigma_{\rm Voigt}\cdot\boldsymbol\epsilon_{\rm eng}\)
is the correct energy density: do **not** double shear energy again.

`core/verification/stress_tensor.py` independently checks:
- symmetric Cauchy stress matrices; normal/shear tensor component mapping;
- arbitrary-plane traction \(\mathbf t=\boldsymbol\sigma\mathbf n\) (unit normal);
- proper coordinate rotations \(\sigma'=R\sigma R^T\);
- principal stress eigenvalue invariance;
- a 45-degree transformation of uniaxial 100 psi tension gives normal and
  complementary shear of 50 psi;
- resultant forces/moments on a *constant*-stress cube cancel; and
- tensor versus engineering-Voigt strain energy density.

Constant-stress cube balance is an algebraic identity, **not** verification of
\(\nabla\cdot\boldsymbol\sigma+\mathbf b=0\) for an arbitrary computed stress
field. The separate `--equilibrium-diagnostics` flag remains available for
sampled numerical 3D balance defects.

### Source-element polynomial stress tests

`core/verification/patch_benchmarks.py` feeds **independent prescribed
kinematic data** into the actual cubic recovered nonlocal element transform,
B operator and material matrix. It samples three different source elements at
two non-center Gauss-like interior coordinates in each.

- Pure longitudinal tension and compression: include **both lateral Poisson
  contractions** so that the analytical stress is *uniaxial* rather than a
  triaxially constrained state.
- Pure \(yz\) shear: expected \(\tau_{yz}=G_{yz}\gamma_{yz}\).
- Pure bending: an explicit quadratic three-dimensional displacement field
  includes lateral Poisson terms, giving \(\sigma_{zz}=E_z\kappa y\) and zero
  remaining 3D stress components at sampled locations.

These tests **pass at near floating-point precision** on 3×3×5. They do not
prove completeness or convergence on arbitrary geometries, and they prescribe
fields instead of solving a loaded boundary-value problem.

## Reproduced local result: 3×3×5 cubic solved research beam

The following numbers were measured on the accompanying archive using SciPy,
not copied from the beam reference:

| Quantity | Direct result | Interpretation |
| --- | ---: | --- |
| Signed FE tip y displacement | −0.114796439 in | Actual solved source displacement |
| Classical beam tip reference | −0.116149091 in | End force plus self-weight |
| Tip agreement discrepancy | 1.164582% | Beam-theory reference difference, not exact 3D error |
| Maximum interior section-shear discrepancy | 18.838% | Source-element sectional integral versus 1D shear |
| Maximum interior section-moment discrepancy | 9.351% | Source-element moment versus 1D moment |
| Direct integrated FE strain energy | 2.224044024 lbf·in | 3D constitutive energy integral |
| Discrete work of applied nodal loads | 4.448088049 lbf·in | Should equal twice FE energy for this linear solve |
| Relative Clapeyron work-energy residual | 1.30×10⁻¹¹ | Excellent **algebraic** consistency, not physical accuracy |
| 1D Euler–Bernoulli reference bending energy | 2.262677132 lbf·in | Different model; not an exact 3D energy solution |

The three interior section stress-cube samples reveal non-negligible FE-vs-1D
shear differences. The program preserves source stresses unchanged instead of
replacing them with smooth analytical beam stresses.

## What remains UNVERIFIED / NOT CERTIFIED

1. No general conforming-space proof across arbitrary distorted/curved meshes.
2. No globally equilibrated stress reconstruction or validated natural-face
   tractions; sampled strong form and face-jump tests remain diagnostics.
3. No rigorous upper bound for full 3D energy, stress, displacement or QoIs,
   including quadrature, roundoff or linear-solve enclosure.
4. No demonstration of convergence in all critical QoIs below a specified
   target on a refined grid; in particular, sectional shear remains weak in
   the 3×3×5 example.
5. No verified design adequacy for real redwood or a replacement for applicable
   codes, grade adjustments, stability/connection checks or professional review.

**Conclusion:** textbook-style mechanics tests significantly strengthen
implementation verification. They do **not** establish equivalence between
1D beam theory and exact 3D elasticity, and they do not confer certification.
