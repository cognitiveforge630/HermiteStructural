# HermiteStructural_CubicHigherOrder_OOP — tracked checklist

Created 2026-10-29 from the user-supplied
`HermiteStructural_CubicHigherOrder(2).zip`.  **Source ZIP preserved.**
This checklist is deliberately conservative: `[x]` means implemented and
checked locally; `[ ]` means incomplete, unverified, or future improvement.

## Phase 1 — Preserve and document baseline

- [x] Create separate `HermiteStructural_CubicHigherOrder_OOP` project; do not modify original ZIP/repository.
- [x] Save baseline Hex27/Hex64 3×3×5 solved and verification outputs before refactoring (`docs/verification/legacy_baseline_3x3x5.json`).
- [x] Document X/Y/Z names, tensor indices, and engineering Voigt ordering.
- [x] Document what the tests can and cannot establish; never conflate beam comparison with certification.

## Phase 2 — OOP element foundations

- [x] Implement `HexElement` common class and `Hex27Element`, `Hex64Element` subclasses.
- [x] Expose eight descriptive corner node getters, full `node_ids`, and local coordinates.
- [x] Introduce `HexFacet` instances accessible via `-Z,+Z,-Y,+Y,-X,+X`.
- [x] Expose facet normals, area, facet node IDs, and tensor-local IDs.
- [x] Cache reference 1D Gaussian points/weights, volume rules, and facet rules by `(p,q,face)`.

## Phase 3 — Mathematical operations

- [x] Move scalar shape functions, shape derivatives, and engineering strain matrix to element mathematics module.
- [x] Implement dedicated 3D body-force quadrature and 2D six-face traction quadrature.
- [x] Implement per-instance local force vector (`clear_force_vector`, `get_force_vector`, `add_force_vector`).
- [x] Implement local `Ke` construction, lazy getter, and shared-stiffness option.
- [x] Verify affine volume and surface Jacobians, six face loads, forces/moments, and reference kernel equivalence.
- [ ] Extend to skewed/curved isoparametric hex geometries (separate validated Jacobian and normals required).
- [ ] Add spatially varying surface traction callbacks, pressure normalization, and tests.

## Phase 4 — Mesh, assembly, and RAM

- [x] Introduce `StructuredNodalHexMesh` lazy per-cell element factory with shared connectivity.
- [x] Introduce `GlobalAssembler` class for sparse K and F assembly.
- [x] Retain shared `Ke` and gravity load for equivalent affine cells.
- [x] Use bounded COO element batches flushed into CSR; avoid all-mesh dense triplet lists.
- [x] Benchmark **small 16-cell assembly only**, Hex27 and Hex64, against legacy-style COO in separate Linux processes (`docs/verification/assembly_benchmark_*_linux.json`).
- [x] Do not use unreliable `__del__` destructors for memory control.
- [ ] Benchmark **192-cell actual process peak RAM**, runtime, and matrix storage against original project in the same environment (and repeat on Windows).
- [ ] Assess assembly strategies and reusable mapped matrix for non-congruent elements.

## Phase 5 — Regression and readable reports

- [x] Extract named one-cell residual, tolerance, and verification-report functions.
- [x] Extract named section-integration, equilibrium, and full-study JSON-report functions.
- [x] Keep existing JSON result keys and non-certified status; add assembly information only if documented.
- [x] Copy frozen legacy math into test-only reference and compare `Ke`, body `Fe`, and +Z `Fe`.
- [x] Test all six face resultant forces and moments for both polynomial orders.
- [x] Test batch-sparse K and F against independently reconstructed legacy COO.
- [x] Run baseline 3×3×5 solved comparisons, rigid modes, patch and symmetry tests.
- [x] Run complete Python unittest suite in this environment (104 run: 102 passed, 2 graphics skips in the 2026-10-29 revision; see docs/AUDIT_20261029.md).
- [ ] Validate the same suite on Windows 11 with the user's Python/SciPy versions.

## Phase 6 — Tutorial and release

- [x] Add short self-contained one-element/face/assembly example `python -m tools.demo_nodal_oop --order 3`.
- [x] Add `docs/OOP_TUTORIAL.md` with formulas, coordinate conventions, and human-readable examples.
- [x] Add VS Code tasks for OOP demonstration and fast OOP tests.
- [x] Regenerate saved cubic-densification plot from published solved source-section JSON, keeping the separate 7×7×17 transverse-refinement point unconnected.
- [x] Repair plot regeneration commands, broken narrative output claims, include Matplotlib dependency, and statically check documentation links and Python entrypoints.
- [x] Add separate |tau_yz| shear comparison to Hex27/Hex64 paired field-map output and section shear/moment to the human-readable convergence report.
- [x] Explain the independent C0 solver's distinction from the legacy Hermite solvers.
- [x] Keep release research-only; no new error-bound or engineering safety certification claim.
- [ ] Consider applying OOP patterns to recovered-rotation Hermite element family **only after** independent math and patch verification.
- [ ] Document measured speed/RAM changes on representative full-size workstation runs (never infer from theoretical COO estimates).

## Acceptance criteria

**Must preserve:** local `Ke/Fe`, global `K/F` (numeric values, allowing FP
summation differences), displacement QoIs, section resultants, reactions,
force/moment and work-energy identities. On tests, no reduced integration,
no silent swapping of face axes, and no alteration of material parameters.

**Current status:** Significant OOP conversion for the *conventional nodal
Hex27/Hex64 path* is implemented and regression-tested. The separate Hermite
recovered-rotation solvers are deliberately **not** refactored. Full-size and Windows peak-RAM benchmarking and non-affine/general traction support remain.

## Mandatory future-work invariant: no centroid stress

**NEVER evaluate stress at element or visualization-cell centroids** for a new
field result, figure, stress peak, or convergence claim. Use original element
corner or near-extremity inset corner samples directly from the FE strain
operator. See [NO_CENTROID_STRESS_POLICY.md](NO_CENTROID_STRESS_POLICY.md).
Any 16-element case is explicitly a smoke test with unvalidated stress values.
