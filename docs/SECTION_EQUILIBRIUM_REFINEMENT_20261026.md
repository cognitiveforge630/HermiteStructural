# Source-element sectional equilibrium — October 26, 2026

**Research status: not a certified a posteriori stress/error estimate.**

This study addresses an ambiguity in the previous `cubic_stress_audit`: the
section resultants were sums of sampled visualization subcell-center stresses,
so they depended on `--subdivisions` even when the underlying FE solution was
unchanged. The new `core.solvers.section_equilibrium.integrate_sections()`
does **not** use those samples. It integrates the original cubic displacement
solution using its element B matrix and the material's Cauchy stiffness D.

For each interior axial section \(z\), with parent \(\zeta=0\), source
cross-section integration calculates

\[
V_y^{FE}(z) = \int_{A_z} \sigma_{yz}^{FE}\,dA,
\qquad M_x^{FE}(z)=\int_{A_z} y\sigma_{zz}^{FE}\,dA.
\]

Fourth-order tensor Gauss integration computes both integrals on the source
FE patches, including the exact physical surface Jacobian. The signed beam
reference (downward force \(P=-12.8\,\mathrm{lbf}\), specific weight
\(\gamma=0.0174\,\mathrm{lbf/in^3}\), cross-section area
\(A=32\,\mathrm{in^2}\), length \(L=120\,\mathrm{in}\)) is

\[
V_y^{beam}(z)=P-\gamma A(L-z), \qquad
M_x^{beam}(z)=-P(L-z)+\tfrac{1}{2}\gamma A(L-z)^2.
\]

The normalized percentage is \(100 |Q_{FE}-Q_{beam}|/|Q_{beam}|\);
the study quotes the **maximum among midspan sections with**
\(L/4\le z\le 3L/4\). This is neither a 3D stress-norm error nor an
error bound. Section locations change with the mesh, which affects maxima.

## Completed, reproducible cubic runs

All rows below are the nonlocal cubic Hermite formulation on Hex8 connectivity.

| Node mesh | FE tip Uy (in) | Tip vs beam (%) | Worst interior shear vs beam (%) | Worst interior moment vs beam (%) |
|---|---:|---:|---:|---:|
| 3×3×5 | −0.114796439 | 1.1646 | 18.8384 | 9.3508 |
| 3×3×9 | −0.115608854 | 0.4651 | 20.5657 | 3.2255 |
| 3×3×17 | −0.116622812 | 0.4079 | 8.0528 | 2.2964 |
| 3×3×33 | −0.117436153 | 1.1081 | 6.1162 | 0.3158 |
| 5×5×9 | −0.116365539 | 0.1864 | 15.0603 | 2.5935 |
| 5×5×17 | −0.116928301 | 0.6709 | 2.0419 | 1.6725 |
| 5×5×33 | −0.117477507 | 1.1437 | 0.3905 | 0.1935 |
| 5×5×65 | −0.117938502 | 1.5406 | 0.0045 | 0.2457 |
| 5×5×129 | −0.118209067 | 1.7736 | 0.0029 | 0.3661 |
| 7×7×17 | −0.117383300 | 1.0626 | 2.0800 | 0.9950 |

All ten used fixed KKT displacement trace **and** fixed nodal rotations,
identical modulus/gravity/traction, SciPy, and fourth-order section
quadrature. The four new section-only runs took about 29 s (5×5×17), 85 s
(5×5×33), 142 s (5×5×65), and 101 s (7×7×17) on the workstation used for
this study. Numerical comparisons are against AWC/beam theory.

**Interpretation:** The 5×5 axial series reduced the worst midspan shear
resultant discrepancy from 2.04% at 17 nodes to 0.39% at 33 and 0.0045% at
65; its moment discrepancy was 0.25% at 65. Refining x and y from 5×5 to
7×7 at 17 axial nodes left shear near 2.08%. Tip-to-beam discrepancy, however,
rose from 0.67% at 5×5×17 to 1.54% at 5×5×65. **These results are
non-monotone across quantities and do not establish stress-field convergence.**
The near-zero 5×5×65 value measures the integrated section force
\(\int_A\sigma_{yz}\,dA\) against a beam-theory resultant. A Galerkin solution
enforces weak force balance, so a good section resultant can coexist with
incorrect pointwise stress and appreciable local \(\nabla\cdot\sigma+b\).
The true 3D orthotropic solution is not supplied; beam agreement is not a
stress error estimate or certification.

The new run artifacts are saved in `samples/cubic_5x5x17*`,
`samples/cubic_5x5x33*`, `samples/cubic_5x5x65*`, and
`samples/cubic_7x7x17*`. Re-run and compare the axial series with:

```powershell
python -m tutorials.cubic --nx 5 --ny 5 --nz 17 --section-only --out outputs/cubic_5x5x17.json
python -m tutorials.cubic --nx 5 --ny 5 --nz 33 --section-only --out outputs/cubic_5x5x33.json
python -m tutorials.cubic --nx 5 --ny 5 --nz 65 --section-only --out outputs/cubic_5x5x65.json
python -m tools.section_refinement outputs/cubic_5x5x17.json outputs/cubic_5x5x33.json outputs/cubic_5x5x65.json --json-out outputs/cubic_5x5_axial_comparison.json
```

The 5×5×65 result is a particularly favorable **section-resultant** comparison,
not evidence that stress has 0.0045% error. A physically admissible stress
reconstruction, validation of local equilibrium and a specified 3D reference
remain necessary for claims about actual stress accuracy.

The [higher-density follow-up](DENSIFICATION_STUDY_20261026.md) extends the
5×5 axial series to 129 nodes: midpoint shear discrepancy falls further to
0.0029%, while tip and moment discrepancies increase to 1.7736% and 0.3661%.
It also checks the effect of sampling sections away from element midpoints.

Cross-section integration was independently checked at Gaussian orders 3,
4 and 5 for the 3×3×5 solution; shear and moment integral differences were
at floating-point noise. An independent cubic constant-strain affine
patch test verifies exact cross-sectional area, axial force and zero
bending/shear for that imposed displacement, again **not** full certification.

### Reproducing with Windows PowerShell

```powershell
python -m unittest discover -s tests -v
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --out outputs/cubic_3x3x9.json
python -m tutorials.cubic --nx 3 --ny 3 --nz 17 --section-only --out outputs/cubic_3x3x17.json
python -m tools.section_refinement outputs/cubic_3x3x9.json outputs/cubic_3x3x17.json --json-out outputs/section_comparison.json
```

The `--section-only` switch avoids creating a 3D visualization sample and
**does not perform the visualization-based Redwood sampled-stress screen**.
The runner still enforces the preliminary beam/reference load screen; neither
screen is a design approval. For detailed stress maps, omit `--section-only`
and choose `--subdivisions` explicitly. A separate `*_source_sections.json`
is written automatically regardless of that selection.

For uniform h-refinement relative to 3×3×9, the completed 5×5×17 case halves
all three element edge lengths. For larger workstation runs, prefer
`--section-only` initially. `--linear-solver pypardiso` is supported if the
user separately installed the optional backend; SciPy remains the default
and required reproducible fallback.

### Boundary conditions are a material numerical uncertainty

The standard case imposes a zero displacement trace on the fixed face **and**
fixes nodal rotations. Research comparisons with free nodal rotation at the
same zero face trace (`--displacement-only-clamp`) gave:

| Mesh | Clamp constraint | Tip vs beam (%) | Worst interior shear vs beam (%) |
|---|---|---:|---:|
| 3×3×5 | Fixed spin | 1.1646 | 18.8384 |
| 3×3×5 | Free spin | 0.1048 | 25.8599 |
| 3×3×9 | Fixed spin | 0.4651 | 20.5657 |
| 3×3×9 | Free spin | 0.9047 | 12.6568 |
| 3×3×17 | Fixed spin | 0.4079 | 8.0528 |
| 3×3×17 | Free spin | 1.4192 | 7.0156 |

These choices define different discrete boundary problems, rather than
interchangeable mesh solutions. True Cauchy displacement-Dirichlet data do
not by themselves prescribe independent nodal rotations. The research
recovery method adds rotational unknowns, so a strong formulation and trace
admissibility analysis is needed before deciding which constraint treatment
faithfully represents the desired boundary condition.

## What this actually closes in the a posteriori programme

The new integrated sections remove the **visualization sampling choice** from
section-resultant diagnostics and quantify how strongly the raw discrete
stress agrees with an independent resultant reference. They do **not**
equilibrate stress. Global reaction equilibrium, interface traction jumps,
body-force residuals, error norm estimation and rigorous upper bounds still
require independent proof/implementation. Do not set `certified_bound=true`
in `docs/CERTIFICATION_STATUS.json` based on these tests.

## Raw 3D local equilibrium residual (new research diagnostic)

The optional flag `--equilibrium-diagnostics` writes a separate report of
sampled `div(sigma) + b` in units of lbf/in³, evaluated by central differences
of the original polynomial element stress. Internal-face `[(sigma)n]` is
integrated directly from both sides using tensor face Gaussian quadrature.

On the present shared-vertex structured cubic elements the internal-face
traction jumps are approximately **10⁻¹³ psi** (floating-point noise); the
raw stress is numerically continuous at these tested interfaces. This does
**not** imply local equilibrium or prove global H(div) admissibility under
natural traction conditions. The sampled **volume** residual is nonzero:

| Mesh | One cell-center residual RMS (lbf/in³) | Eight Gauss points/cell residual RMS (lbf/in³) |
|---|---:|---:|
| 3×3×9 | 0.8521 | 1.8281 |
| 3×3×17 | 0.7033 | 1.5139 |
| 3×3×33 | 0.6627 | 1.3685 |
| 5×5×9 | 1.2780 | Not calculated |

The *large sampling dependence* makes clear that midpoint-only residuals are
not a verified norm. The denser calculation is a **numerical approximation**,
not a validated integral, reliability guarantee, or true-solution error
percentage. Small central-difference steps of 1e−4 versus 5e−5 agreed closely
for 3×3×9 and 3×3×17 at fixed residual sampling locations, but this is not
a rigorous truncation/roundoff enclosure. The finite-element continuity
checks and quadratic manufactured `u_y=α z²`, `theta_x=−α z` strong-form
balance test pass without presuming that stress is equilibrated for the
loaded beam.

```powershell
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --equilibrium-diagnostics --out outputs/cubic_3x3x9.json
Get-Content outputs/cubic_3x3x9_equilibrium_diagnostics.json
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --section-only --equilibrium-diagnostics --volume-residual-order 2 --out outputs/cubic_3x3x9_8point.json
```

**Next mathematical obstacle:** obtain a stress field with verified
`div σ = −b` and the prescribed exterior tractions, or construct admissible
corrections to the existing continuous but locally unbalanced stress. Only
then can a complementary-energy/hypercircle-type guarantee even be considered,
subject to the required displacement admissibility, functional framework,
exact integration/roundoff treatment, and proof. Do not treat any sampled
residual as a certification bound.
