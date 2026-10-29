## 2026-10-29 — publication-readiness re-audit (research release)

- Independently re-ran 114 existing tests and the theoretical/load, provenance,
  numerical benchmark and documentation gate checks on the supplied archive;
  recorded outstanding scientific and Windows-GUI validation in
  `docs/PUBLICATION_AUDIT_20261029.md`.
- Retired the legacy element-center stress-divergence diagnostic; residual checks now sample non-centroid 2×2×2 or 4×4×4 points with regression guards.
- Hardened cubic solved-DOF checkpoint v2: store the full 6×6 orthotropic
  matrix and reject changed constitutive data when replaying stresses.
  Retained the shipped smoke DOFs and updated only their archive provenance.
- Added a negative regression for changed transverse/shear constants, and
  prominently distinguished section-resultant agreement from local-stress
  convergence and solver accuracy in the README.
- No production-size stress result was fabricated or rerun, no E calibration,
  and no mathematical FE kernel or load formulation was changed.

## 2026-10-29 — data traceability and documentation/command audit

- Rebuilt the Hermite densification SVG from saved solved JSON. Added 7×7×17 as an unconnected transverse-refinement diamond, not a fictitious 5×5×257 axial solve.
- Added `tools.plot_cubic_densification` and `tools.plot_p_order_comparison`; replaced the stale command to `outputs/plot_p_order_mesh_comparison.py` with a working tracked module. The Hex27/Hex64 192-cell plot is based on a labeled saved summary; raw 192-cell archives are not included.
- Included new small 3×3×5 Hex27 and Hex64 solved field NPZ archives, JSON reports and source-section CSV for immediate viewer inspection. Eliminated an unused byte-for-byte duplicate image folder named `pictyres/`; retained `pics/`.
- Added `tools.audit_documentation` and regression tests for local links, CLI entrypoints and scientific plot provenance.
- Added `matplotlib` as a requirement, clarified generate-versus-shipped `outputs/` paths, removed stale README heading placement, and unified the duplicated root OOP documents behind canonical docs links.
- Expanded the nodal convergence human-readable report with sampled section shear and moment and the FE/beam slice image with absolute `tau_yz` shear stress; this is not a certified pointwise stress bound.
- See `docs/AUDIT_20261029.md` for validation evidence, limitations and follow-up gates.

## 2026-10-29 — Independent Hex27/Hex64 OOP refactor

- Created separate `_OOP` project; source ZIP and recovered-rotation Hermite
  numerical algorithms remain untouched by this nodal refactor.
- Added element/facet OOP, named faces, cached Gauss rules, lazy shared Ke,
  streamed `GlobalAssembler`, and explicit quadrature-integrated traction.
- Extracted readable section, equilibrium, and JSON reporting functions.
- Frozen legacy Hex27/Hex64 element math and solved baseline for regression.
- Added six-facet force/moment regression, mesh/assembly checks, tutorial,
  PowerShell commands, VS Code tasks, and tracked `docs/OOP_TODO.md`.
- No certified error bound; large 192-cell run was not performed.

## 2026-10-26 — Clarify shear-comparison scope and mechanics references

- Defined the reported 1.669% shear discrepancy beside the README result and benchmark table, including its discrete norm, sampling region and equilibrium limitation. Added 3×5×33/41/49 axial refinements: the sampled difference falls below 1% at 41 nodes, while tip disagreement rises; signed projections show a small, mesh-dependent aggregate over- or under-shoot.
- Labeled the 5×5×129 shear figure in the README as a single axial-midpoint section-resultant comparison and linked the local-balance and transverse-warping evidence.
- Added the signed two-view shear guide to the documentation index and clarified that the mechanics audit uses Beer's mechanics-of-materials stress conventions; these checks remain research verification, not certification.

## 2026-10-27 — Two-view signed shear stress in PyVista

- New `--shear-comparison` renders synchronized classical `VQ/(It)` and raw source-Hermite `tau_yz` with a common signed scale, cell-wireframe and linked camera.
- Separate FE and reference stresses are evaluated at matching one-sided **original source-element corner** positions; dual VTU duplicates display-cell points to retain discontinuities without false nodal averaging.
- Added `--shear-screenshot`, `--shear-export-only`, and saved-VTU standalone viewer (`python -m core.visualization.shear_comparison ...`).
- Added `docs/DUAL_SHEAR_PYVISTA.md` and 5 no-GUI regression tests. No solver/B-matrix/elasticity change. Theoretical shear is a first-order beam reference, not an exact 3D or certified result.

## 2026-10-27 — Transverse cubic shear-warping completeness and independently resolved shear

- Established root cause of misleading 3×3×9 shear benchmark: cubic `u_z=y^3` warping is not reproduced with only three transverse `y` nodal layers. Full six-strain patch test now checks the source cubic `B @ T` directly.
- Reproduced solved 3×5×25 original-source-vertex **raw** shear L2 1.669% versus the beam reference (3×3×9: 73.875%), plus sectional shear 0.604%; 3×4×25 also passes 5% nodal shear.
- Added optional `--shear-warping-audit`, fail-visible `ny=3` warning, separate kinematic-patch JSON, new regression tests and benchmark report.
- **Did not** alter cubic element stiffness, formulation, material, loads, constraints, raw FE stresses, or make a certification claim.
- Recorded that penalty traction stabilization and beam-like stress substitution were rejected as fixes for the actual FE shear field.

## 2026-10-27 — Raw cubic shear/moment consistency investigation

- Added `--shear-balance` and independent axial-shear diagnostic without modifying
  solver stiffness, strain/stress operator, geometry, loading or BCs.
- Reports pointwise sampled, cell-averaged section shear and interpolated
  derivative of FE bending moment **separately**, with conservative research-only status.
- Reproduced 3×3×9 discrepancy and compared fixed versus displacement-only
  clamps; new reports avoid confusing averaged errors with local stress accuracy.
- Added four numerical regression tests and detailed engineering handoff.

## 2026-10-26 — Independent Mechanics of Materials verification layer

- Added exact 1D signed load, bending-stress, shear-profile and energy
  references in `core/reference/mechanics_of_materials.py` (no calibration).
- Added independent six-component Cauchy stress cube, transformation,
  traction, principal stress, constant-stress cube balance and engineering
  shear energy checks.
- Added prescribed physical cubic source stress patches for tension,
  compression, pure yz shear and Poisson-consistent quadratic pure bending.
- Added opt-in `--mechanics-audit` in cubic runner, producing direct 3D FE
  stress/centerline samples, distinct 1D theory values, and source-Gauss
  Clapeyron work-energy consistency (not a posteriori certification).
- Added 10 unit tests and CLI JSON checks; documented reproduced 3×3×5
  solve and explicit certification limitations.

## 2026-10-26 — Higher-density cubic section-resultant study

- Extended the fixed-5×5 axial series through 129 nodes: midpoint shear
  discrepancy versus beam theory decreased to 0.002908%, compared with
  20.565722% on 3×3×9. Tip and moment discrepancies remain separately reported.
- Added saved 5×5×129 results, a README densification chart and a detailed
  higher-density study with solver diagnostics and reproducible commands.
- Added configurable source-section axial positions; a three-slice 5×5×65
  check measured a maximum shear-resultant discrepancy of 0.022961%.
- Validated the CLI and existing suite: 47 tests ran, one graphics check
  skipped. These are measured section-resultant trends, not certified 3D
  stress errors.

## 2026-10-26 — Research source-stress equilibrium diagnostic

- Added internal-face traction-jump audit from unmodified FE stress; random
  nodal-data checks also show near-zero jumps on the structured cubic grid.
- Added finite-difference sampled strong-form `div(sigma)+body_force` balance
  residual diagnostic and independent affine/quadratic manufactured tests.
- Separated midpoint versus 8-Gauss-point residual results, emphasizing
  strong sampling sensitivity and lack of a mathematical upper bound.
- Added optional `--equilibrium-diagnostics` and
  `--volume-residual-order {1,2}`; default remains economical midpoint
  sampling, and all diagnostics remain explicitly uncertified.

## 2026-10-26 — Source FE sectional integration and controlled cubic refinements

- Added direct quadrature section-resultant integration of source-element stress,
  separated from visualization-subcell sampling and stored with every cubic run.
- Added `--section-only` for inexpensive non-visual QoI output and
  `tools.section_refinement` for comparison of equivalent research cases.
- Measured 3×3×17 and 3×3×33 axial and 5×5×9 transverse/axial cubic
  refinements; the tip and sectional errors are non-monotonic, so no claim of
  universal refinement improvement is made.
- Recorded a displacement-only clamp comparison explicitly separately from
  the default fixed-rotation benchmark.
- Added source-section and comparator unit tests; kept certification false.
- A 5×5×17 local run exceeded the two-minute limit and is not reported as a
  completed benchmark.

# Change log

## 2026-10-26 — Cubic-refinement and a posteriori readiness documentation

- Added a README section explaining a **controlled cubic mesh study**, with
  exact PowerShell commands and explicit (uncertified) acceptance targets.
- Added `docs/CUBIC_REFINEMENT_AND_APOSTERIORI.md`, distinguishing beam-theory
  discrepancies, mesh-to-mesh QoI changes, and rigorous a posteriori bounds.
- Added `tools/cubic_refinement.py` for offline comparison of cubic runs; it
  validates common loads, E, beam reference and boundary conditions, and
  warns when the stress-sampling settings differ.
- New stress audits save `sample_subdivisions` explicitly; old archived audits
  remain unmodified, with their sampling inferred from recorded section counts.
- Expanded the fail-closed certification status and added tool unit tests.
  No expensive refined-mesh solve was run and no certification claim was added.

## 2026-10-26 — Presentation and Redwood research cleanup

- Kept `E_L=1,100,000 psi` and the published Redwood design-value table.
- One canonical root `README.md`: all material, deflection, cubic stress and
  research-limitation tables; moved the documentation navigation to `docs/INDEX.md`.
- Corrected default displacement to nonnegative **vertical magnitude `|Uy|`**;
  blue is zero/least, red maximum/greatest. Signed `Uy` remains available on `U`.
- Explicit white 10-point scalar bar and HUD, with field, units, and color meaning.
- Added yellow cube support markers on original fixed-end FE nodes, and 3D-anchored
  `FIXED` / `FREE / LOADED` labels for all three viewer formulations.
- Curated VTUs contain original source-node flags and a new positive-magnitude
  field; displayed Hermite displacement/stress numbers were not re-calibrated.
- Re-pinned source-integrity hashes after the intentional viewer and export edits.

**Not claimed:** proof of complete three-dimensional stress accuracy, structural
design approval, or a mathematically certified a posteriori error bound.

## 2026-10-26 — Distribution audit corrections

- Fixed two independently reproducible release-check failures (duplicate README
  and out-of-date source hash). Removed `docs/README.md`; keep `docs/INDEX.md`.
- Made `|Uy|` the **actual** default in both the sampler and tutorial CLI;
  formerly the changelog claimed this without the runtime default changing.
- Added fixed/free source-node flags and absolute-value fields to all precomputed
  VTUs, and to new exports from the shared sampler.
- Hardened point-label actor lifecycle so PyVista actor-map names are not assumed.
- Added numerical/VTK assertions and documented Windows-only GUI checks.
- No solver stiffness/constitutive calibration or displacement results were changed.

## 2026-10-29 — Verification class extraction and editor cleanup
- Extracted rectangular-cell polynomial and one-cell mechanics checks to `core/verification/nodal_hex.py`, exposed through `PolynomialPatchVerifier` and `ElementMechanicsVerifier`.
- Preserved the `tools.p_order_study` compatibility functions and `--mode verify` JSON output.
- Removed unused Hex27/Hex64 and 1D Lagrange imports and simplified the optional AWC JSON comparison check to avoid ambiguous nested access warnings.
- Added `tests/test_nodal_verifier_classes.py` and `docs/VERIFICATION_ARCHITECTURE.md`.
- No changes to element stiffness integration, Gaussian rules, material parameters, or global assembly.

### 2026-10-29 — Reproducible nodal-hex convergence workflow
- Added `tools.convergence_study` with mesh/solve/report/view/all stages and conservative >8000-DOF solve protection.
- Added optional raw solved-field archival (`.npz`) to `solve_study` without changing its stiffness or linear algebra.
- Added FE source-center tension/compression/von-Mises and nodal displacement maps alongside labeled approximate Euler–Bernoulli reference maps.
- Added per-mesh Markdown/JSON reports and a full step-by-step PowerShell manual in `docs/CONVERGENCE_STUDY.md`.
- Added regression tests and a three-level Hex27 example. All 93 regression tests passed (2 skipped).
