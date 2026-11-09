# Cubic Hermite stress recovery and Redwood screening

**Research status: longitudinal bending stress promising; complete 3D stress distribution NOT certified.**
All distances in inches, stresses in psi, forces in lbf. The demo is an *idealized*
4 × 8 × 120 solid cantilever, grain parallel to global +z, distributed downward
end traction `P = −12.8 lbf`, and self-weight `γ = 0.0174 lbf/in³`. Gravity is
integrated over the beam; the end force is integrated over the end face.
Deformations in PyVista are **visually magnified** unless `--deform-scale 1` is used.

## Material basis — revised to published Redwood properties

The original 2,000,000 psi longitudinal modulus was a prescribed *hypothetical*
Redwood-like parameter, not tied to a dimension-lumber grade. This release instead
uses the published **Redwood Select Structural, dimension lumber** reference:

| Property | Unadjusted published value |
|---|---:|
| Longitudinal modulus `E_L` | **1,100,000 psi** |
| Bending reference `F_b` | **1,100 psi** |
| Tension parallel `F_t` | **625 psi** |
| Shear parallel `F_v` | **160 psi** |
| Compression parallel `F_c` | **1,100 psi** |
| Compression perpendicular `F_c⊥` | **425 psi** |

Source: *Redwood Technical Data Sheet, Structural Design Values*, table for
2–4 in thick dimension lumber, Select Structural,
https://getredwood.com/wp-content/uploads/2024/07/2022-Structural-Design-Values.pdf .
These are **reference design values, not maximum physical failure strengths**.
Actual timbers require a species/grade/size and all relevant design adjustments;
the model uses idealized full 4 × 8 in geometry, not actual dressed dimensions.
The other orthotropic elastic constants retain the project's *research ratios*
from `core/materials/redwood.py`: they are not a grade-certified full constitutive
set. Tension **perpendicular** to grain is NOT covered by this reference table;
this must not be interpreted as a passed safety check. No connection, bearing,
buckling, creep, stability, duration, moisture or full NDS design check is made.

## Pre-solve screening and saved FE stress screening

The baseline load passes preliminary theoretical maximum bending, tension,
compression, shear, and a separately declared **L/360 serviceability example**.
**L/360 is NOT a Redwood fracture/displacement limit.** The preload screen
fails closed for larger input loads that exceed any one of those reference
checks (`--end-force` supports exploring other *downward* end forces). The cubic
run additionally checks *sampled* source-FE longitudinal and shear components
against the quoted reference stresses. Peak stresses between samples and clamp
singularities are NOT bounded; these checks are not a structural approval.
Von Mises is retained for generic visualization only, **not** as a valid wood
failure criterion.

## Why this cubic stress field is a meaningful improvement

Earlier viewers produced actual `σ = D ε` at **subcell centers**. The current
viewer evaluates at subcell corners and computes display-cell means; smooth color
bands could falsely suggest accuracy. In the revised viewer, raw FE stress is
kept separate from the beam-theory stress and the **signed FE-minus-reference**
field. No analytical/reference values overwrite the finite-element stress.

At each source element's material point `(ξ,η,ζ)`:

    ε_H = B_H(ξ,η,ζ) T_element U_support
    σ_H = D_orthotropic ε_H

The nonlocal recovery matrix `T_element` contains cubic Taylor recovery,
including the mixed `xyz` derivative; this is *cubic recovery*, not an ordinary
cubic Lagrange Hex64 or a mathematically proven C¹ Hermite element. For the
independent Saint-Venant bending comparison, with `q = γ b h`:

    M_x(z) = |P|(L-z) + (q/2)(L-z)²
    σ_zz,beam(y,z) = M_x(z)y / I_x,  I_x = b h³ / 12
    V_y(z) = P - q(L-z)
    τ_yz,beam(y,z) ≈ (3/2)[V_y(z)/(bh)] [1-(2y/h)²]

The shear expression is an approximate *rectangular-beam* reference, NOT a
3D traction solution for all sections. Comparisons near the clamp/loaded end
are vulnerable to 3D end effects and are excluded from the quoted interior L²
stress discrepancy; they are **not removed from the stress visualization or
material screening**. Per-section FE bending moment and shear resultants are
reported independently, including signed axial force.

## Numerical diagnostics (identical new material and loads)

| Cubic node mesh | FE tip `u_y` | Beam tip | Tip discrepancy | Interior `σ_zz` relative L² discrepancy |
|---|---:|---:|---:|---:|
| 3×3×5 | −0.114796 in | −0.116149 in | 1.165% | 9.12% (subdivisions 2) |
| 3×3×9 | −0.115609 in | −0.116149 in | 0.465% | 2.16% (subdivisions 3) |

**Historical metrics in the preceding table** were evaluated from raw source-FE
stresses at display subcell centers, spanning 1/4 to 3/4 of beam length.
The current default uses raw individual corner-side samples with the same span
window; do not put new corner-norm results in that old table without relabeling.
This is **not** a
certified a posteriori error bound, and is sampling-resolution dependent.
At 3×3×9 and subdivisions 2, the interior error was about 1.16%; changing
sampling to 3 subdivisions gives about 2.16%. These differences are a reason
for caution, NOT grounds to select whichever number looks best.

At 3×3×9 with subdivisions 2, the maximum interior sectional bending moment
comparison error was about 7.3%; the corresponding shear-force resultant
comparison was about **12.2%**. Thus, even though `σ_zz` broadly resembles
classical bending, **the complete stress state does not yet meet a 5% target**.
The current work supports a useful **cubic bending-stress research viewer**,
not an accurate/certified full 3D stress solution.

## References in the supplied books

* **Zienkiewicz, Taylor & Govindjee**, *The Finite Element Method: Its Basis and
  Fundamentals*, 8th ed., **§10.8** (reporting displacement/strain/stress),
  **§§11.3–11.8** (patch tests and higher-order robustness), **§12.3** (mixed
  methods/locking), **Appendix I, pp. 575–587** (stress recovery, optimal
  sampling, superconvergent patch recovery, boundary/interface recovery).
  Appendix I's polynomial least-squares fit is a research lead for *future*
  equilibrated stress reconstruction, not proof that the present cubic
  displacement recovery already provides superconvergent stresses.
* **Sadd**, *Elasticity: Theory, Applications, and Numerics*, 5th ed.,
  **§§3.6–3.7** (stress contours/equilibrium), **§9.8–9.9** (flexure), **§11**
  (anisotropic elasticity), **§16** (numerical finite elements).
* **Timoshenko & Goodier**, *Theory of Elasticity*, 2nd ed., **Ch. 12**
  (bending of prismatic bars): a reference for physical stress distributions.
* **Atanackovic & Guran**, *Theory of Elasticity for Scientists and Engineers*,
  **§§5.7–5.8** (pure bending and terminal load), **§3.3** (orthotropy),
  **§7** (energy and variational principles).
* **Eringen**, *Microcontinuum Field Theories I*, **Ch. 5** (micropolar
  elasticity/stress/couple stress and balance laws). Important distinction:
  independently rotating microrotation DOFs in true Eringen/Cosserat mechanics
  imply different field equations and couple-stresses. The current
  *recovered-rotation Hermite* beam is a classical orthotropic displacement-
  energy model with recovery DOFs, not automatically an Eringen material.

## Reproduce the cubic research result and viewer

```powershell
python -m pip install -r requirements.txt
python -m tutorials.cubic --nx 3 --ny 3 --nz 9 --subdivisions 3 --export-vtu --view --full-screen --field "Sigma zz (psi)"
python -m core.visualization.pyvista_viewer outputs/cubic.vtu --field "Cubic minus beam sigma zz (psi)"
Get-Content outputs/cubic_stress_audit.json
```

PyVista: `6` displays raw cubic **σzz**; `0` shows **beam reference σzz**;
`B` (lowercase `b` key) shows **cubic-minus-beam σzz**;
`T` displays **longitudinal tension utilization**. `C` toggles a midspan
cross-section; `D` toggles deformed shape; `S` screenshot. Some shortcuts require
an exported cubic VTU with the corresponding extra arrays. Keep magnification
in the HUD in mind and use `--deform-scale 1` for physical deformation.

Saved example: `samples/cubic_3x3x9.vtu` can be opened *without a solve*;
`samples/cubic_3x3x9_stress_audit.json` and the matching `_sections.csv`
contain the section-by-section verification data. This viewer does NOT convert
FE stresses into analytical stresses and does not claim certification.
