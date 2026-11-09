# Publication-readiness audit — 29 October 2026

**Scope:** `HermiteStructural_CubicHigherOrder_OOP(8).zip`, repaired here only
where supported by review and fast regression testing. This is a technical
code review, **not an independent certification of a finite-element method**.
All findings below concern the actual supplied archive, not a prior branch.

## Bottom line

**Acceptable to publish as experimental/research/tutorial code if the research
limitations are prominent. NOT ready for a public claim of 1% accurate local
3D stress, verified peak stress, or independently validated structural design.**
Theoretical calculations are correct for the stated idealized loads and
geometry; gravity is included in both solver load vectors and beam theory.
The strongest historical 0.002908% result is a *section shear resultant*.
The corresponding local shear-stress field was not saved or solved for the
5×5×129 case. The 3.732 psi theoretical maximum has not been matched against
5×5×129 local FE shear stress.

## Prioritized findings

| Priority | Finding and evidence | Disposition |
|---|---|---|
| **P0 — claims** | No 5×5×129 local-stress solution archive or local-shear discrepancy exists in the supplied ZIP. No 1% local-stress result can be claimed. | Clearly labeled **NOT MEASURED**; new solve and mesh study required. |
| **P0 — numerical reliability** | Saved 5×5×129 KKT solve: RHS-relative infinity residual **0.002739826** (0.274%), normwise backward error **4.20×10⁻¹⁷**. A tiny backward error relative to matrix norm can coexist with a large forward error for a poorly conditioned system. | Original saved diagnostic retained, but publication no longer implies that backward error certifies the 0.002908% section metric. Evaluate scaling, condition sensitivity, independent solver refinement. |
| **P0 — reference mismatch** | Classical 3.732 psi is maximum *rectangular beam shear* at clamp centerline. Saved FEM point samples are inset and generally elsewhere. 129.96 psi is nominal outer-fiber bending, not max 3D von Mises. | README makes both differences explicit. Do not compare unlike extrema to claim 1%. |
| **P1 — centroid diagnostic loophole (fixed)** | The optional strong-form stress-divergence diagnostic previously permitted a one-point center-only evaluation, inconsistent with the no-centroid research-output policy. | **Fixed:** residual sampling now uses even Gauss orders 2 or 4; explicit centroid stress-derivative evaluations rejected by regression. Assembly and section integration were not altered. |
| **P1 — reproducibility defect (fixed)** | Full-DOF cubic checkpoint previously recorded only longitudinal $E_z$, not the full orthotropic constitutive law; changes to transverse or shear moduli could silently alter replayed stress without changing $E_z$. | **Fixed:** version 2 stores/checks complete 6×6 matrix; new negative test. Shipped smoke DOFs preserved and checkpoint metadata updated. Older v1 checkpoints deliberately rejected. |
| **P1 — physical modeling assumption** | Cubic default enforces zero face displacement **and fixed nodal rotations**. That is not automatically equivalent to classical displacement-only 3D Cauchy clamp. Alternative `--displacement-only-clamp` exists; no demonstrated equivalence across the dense study. | Disclose BC choice, independently compare before attributing beam-theory deviations solely to FEM discretization. |
| **P1 — sampled stress interpretation** | Hermite legacy VTUs use arithmetic means of eight directly evaluated subcell-corner stresses for *cell coloring*, while the new cubic replay exports eight separate one-sided point stresses. Nodal Hex27/64 maps use corner extrema. These are distinct processing pipelines and cannot be compared as equivalent stress maxima. | Documented. No actual centroids are evaluated in these paths. Never treat a mean as a peak. |
| **P1 — historical study design** | Axial refinement 5×5×9→129 does not monotonically reduce tip displacement vs beam theory; 5×5×129 tip discrepancy is **1.773562%**. Section positions change with mesh; off-midpoint 5×5×129 accuracy unmeasured. | Chart and text correctly retain both trends; do not call the largest mesh 'most accurate' across QoIs. |
| **P1 — validation portability** | PyVista absent in this audit environment; two GUI tests skipped. Windows 11 GUI, VS Code tasks, and expensive 192-source-element runs not executed end to end. | **Open validation gate.** Static existence checks are not execution of every documented command. |
| **P2 — publishing rights** | No `LICENSE` found. | Author must choose license if public reuse is intended. No license invented. |
| **P2 — dependency repeatability** | `requirements.txt` uses version ranges, not a tested Python/NumPy/SciPy/VTK lockfile or Windows environment report. | State platform and freeze known-good versions after Windows validation. No claim of Windows compatibility test. |
| **P2 — terminology** | Redwood elastic anisotropy ratios are a research model and approximate one-dimensional beam theory is **not** an exact three-dimensional orthotropic elasticity solution. | Keep tutorial/portfolio framing, prohibit structural certification or design safety inference. |

## Independent source checks and reproducibility evidence

From the uploaded archive prior to this patch:

- `python -m unittest discover -s tests -q`: **114 tests, 2 skipped**.
- `python -m tools.audit_documentation`: **37 Markdown files, 104 internal links and 274 Python command references**, no missing targets/entrypoints. This is a *static* check, not a command-execution proof.
- `python -m tests.checks.source_provenance`: **64 matching recorded active source hashes** (identity, not correctness).
- `python -m tools.mechanics_benchmarks`: nine research benchmark checks passed; status remains **NOT CERTIFIED**.
- `python -m tests.checks.certification_readiness`: explicitly reports no 3D stress error certificate.
- `python -m tools.plot_cubic_densification`: regenerated the SVG from committed numerical JSONs.
- Python compile check succeeded. VTK readable without PyVista; actual desktop renderer not tested.
- v2-specific test verifies a checkpoint with altered shear stiffness is rejected even when $E_z$ is unchanged.
- **Final revised build:** 116 tests run (114 passed, 2 GUI skips); 38 Markdown files, 108 links and 279 Python command references checked; 64 source hash checks passed; Python compilation succeeded. The new centroid-diagnostic regression passed.

### Confirmed default theoretical benchmark

| Input/output | Recomputed value |
|---|---:|
| Width × depth × length | 4 × 8 × 120 in |
| Gravity force density | 0.0174 lbf/in³ |
| End load | −12.8 lbf |
| Total distributed self-weight | 66.816 lbf |
| Total downward load | 79.616 lbf |
| Gravity-plus-end fixed-end moment magnitude | 5,544.96 lbf·in |
| Nominal outer-fiber bending stress | 129.96 psi |
| Classical clamp-centerline max shear | 3.732 psi |
| Euler–Bernoulli tip deflection | −0.116149091 in |
| Historic 5×5×129 FEM tip displacement | −0.118209067 in |

**Local-stress validity warning:** 16-element Hex27/Hex64 samples have
maxima that depend strongly on inset sampling location. The recorded Hex64
von Mises maximum of 129.693 psi happens to be near 129.96 psi nominal bending,
but **these are different stress measures**. At one stored Hex64 maximum
longitudinal-stress sample, $|\sigma_{zz}|=135.982$ psi versus a local beam
reference of 128.415 psi. Neither comparison establishes convergence.

## Publication gate — outstanding manual/expensive checks

1. Run `python -m unittest discover -s tests -v` and the [release checklist](RELEASE_CHECKLIST.md) on a supported Windows desktop with PyVista installed; record Python and dependency versions and ensure all viewer shortcuts work.
2. Run 5×5×129 **with `--save-solution`**, inspect its KKT linear residual and independently assess conditioning/forward accuracy. Do **not** substitute historical JSON for saved DOFs.
3. Repeat local inset-corner stress replay on 5×5×33, 5×5×65, and 5×5×129. Compare at **fixed physical coordinates** in the interior and study local/clamp peaks separately. Report $\sigma_{zz}$, $\tau_{yz}$ and von Mises separately.
4. Run the displacement-only clamp sensitivity and evaluate convergence to a 3D reference, not solely approximate Euler–Bernoulli theory.
5. Measure full production Hex27/Hex64 memory and runtime on the intended machine; do not infer from a 16-element example.
6. Consider independent scientific peer review and select a `LICENSE` if redistribution/reuse is intended.

### Release-oriented fast commands (project root)

```powershell
python -m unittest discover -s tests -v
python -m tools.audit_documentation
python -m tests.checks.source_provenance
python -m tools.mechanics_benchmarks
python -m tests.checks.certification_readiness
python -m tools.plot_cubic_densification
python -m unittest tests.test_cubic_replay_fields tests.test_gravity_reference_contract -v
```

**Do not interpret green fast checks as a passed structural validation gate.**
