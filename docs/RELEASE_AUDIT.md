# HermiteStructural public-sharing audit — 2026-10-26

**Purpose:** reduce avoidable questions about implementation accuracy,
scientific claims and reproducibility. This is an engineering-research
software review, **not** an independent numerical certification or code audit.

## Issues found in the supplied 2026-10-26 archive and corrected

| Original issue | Disposition in revised archive |
|---|---|
| `docs/README.md` duplicated the root README, causing the one-README unit test to fail. | Removed `docs/README.md`; `docs/INDEX.md` remains a documentation index. |
| The source-integrity check failed on `pyvista_viewer.py`, despite the change log claiming a successful source pin. | Recomputed hashes for the reviewed active source; source identity check now passes. |
| Release notes described `|Uy|` as the default, but `sampling.DEFAULT_FIELD` and tutorial `--field` defaulted to signed `Uy`. | Both now actually default to nonnegative vertical-displacement magnitude; signed `Uy` remains available. |
| Saved VTUs lacked the recorded source-node flags and absolute stress/displacement fields. | Updated all four provided VTUs and the sampler for new runs; VTK metadata checks verify nine fixed and nine free original nodes per source mesh. |
| Support labels previously caused PyVista actor-name lookup exceptions. | The renderer now retains returned actors instead of indexing the actor-name map for 2D labels; actual GUI check is still pending. |
| An internal code comment described the end-face load as divided equally between nodes. | Corrected: the load vector is assembled by integrated consistent Hermite face traction. |
| Stress statistics could be mistaken for certified global stress accuracy. | README states these are interior σzz *beam-theory comparisons only*. The shear discrepancies at visualization subdivisions 2 and 3 are 12.23% and 19.82%, respectively. |

## Verification performed on the revised source

- Complete fast unit-test discovery; a few graphics-specific tests require
  PyVista and a working Windows graphics stack.
- Three-formulation rigid/affine/operator checks and a random interelement
  matched-face displacement continuity test on a small structured grid.
- Source SHA-256 snapshot consistency and a duplicate-README check.
- All four saved VTUs read by VTK, with field magnitudes, original-node
  markers, displacement sign/magnitude and cell-stress identities checked.
- **Fresh cubic 3×3×9 solve** (not just loading old JSON) reproduces
  −0.115608854-in tip `Uy`, −0.116149091-in beam reference, and
  2.162516% interior σzz relative-L² beam discrepancy.
- Verified Redwood grade design values against the linked technical data
  sheet, but **not** a complete graded-timber constitutive/failure model.

## Remaining gates before calling this public-facing package validated

1. **Run the Windows/PyVista manual checklist** in [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md), especially keyboard field switching and named support labels. There was no PyVista installation available to this review environment. Syntax tests cannot substitute for this.
2. **Select a software license** if you want other developers to reuse the source; no license was assumed here.
3. Do not claim the 3D stress distribution is accurate or certified based on the current 2.16% *interior longitudinal* comparison. Shear and boundary effects remain unresolved. No mathematically guaranteed a posteriori bound is available.
4. No independent physical experiment, full h/p mesh-convergence study or expert review of the anisotropic design assumptions has been performed here. Linear Hermite is notably over-stiff on the coarse benchmark.
5. If publishing as an engineering service, separately evaluate professional licensing, material design adjustments, connections, stability, load cases, and other application-specific requirements.

**Appropriate description:** "Python implementation and visualization of
experimental six-DOF recovered-rotation Hermite hexahedra, with reproducible
cantilever benchmarks, material-aware preliminary screening, and documented
limitations." Do not describe the package as "certified structural analysis".
