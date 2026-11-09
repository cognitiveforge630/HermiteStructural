# Migration from the uploaded HermiteStructural-main(2).zip

This distribution is a **replacement clean working tree**, not an in-place Git commit. Back up the existing clone before replacing its files and review a `git diff` prior to commit/push. The original uploaded archive is unchanged.

## Migrated active code

| Before | Now |
|---|---|
| `legacy_hermite/FEMHermiteBeamRegion.py` | `core/elements/base.py` |
| `legacy_hermite/recovered_rotation_hex8.py` | `core/elements/linear.py` |
| `legacy_hermite/quadratic_rotation_hex8.py` | `core/elements/quadratic.py` |
| `legacy_hermite/cubic_rotation_hex8.py` | `core/elements/cubic.py` |
| `legacy_hermite/face_dirichlet.py` | `core/solvers/constraints.py` |
| `legacy_hermite/linear_solvers.py` | `core/solvers/linear.py` |
| Parts of `legacy_hermite/solver.py` | `core/solvers/beam.py`, `core/materials/redwood.py`, `core/solvers/loads.py`, `core/reference/awc.py` |
| `run_higher_order_hermite.py` numerical runner | `tutorials/cantilever.py`; root CLI remains a very small wrapper |
| `orthotropic_companion/` matrix calculation | `studies/orthotropic/comparison.py` and `tests/test_material.py` |

The common Hermite base-class mesh generation now calls dedicated `core/mesh/structured.py` functions, preserving the old node/connectivity ordering. The nonlocal interpolation, primitive matrices, element-transforms, stiffness assembly and face KKT constraint numerical kernels were migrated rather than rewritten. This minimizes scientific changes during a folder refactor.

## Deliberately excluded

Old first-order Cauchy-specific Tutorial001 wrappers, outdated mesh viewers expecting previous file formats, historical session-5/6/7/8 scripts, unused reference/enrichment classes, `calibrated_longitudinal_E` and other fitted-modulus pathways, old compatibility launchers, binary caches and previous results were not imported into the clean working source tree. Those items remain retrievable from the **original uploaded ZIP or Git history**, not from an ambiguously named folder in this current release.

Some old numerical certification scripts referenced absent `aposteriori/` modules and saved reports not included in the uploaded ZIP; carrying them forward as functional checks would be misleading. The new `tests/checks/certification_readiness.py` truthfully reports outstanding proof obligations, rather than silently claiming the old checks succeed.

## What was actually validated here?

Thirteen fast unit checks; source identity; mode/finite-operator checks; three real 3×3×5-node solves; and one real 4×5×17 quadratic reference solve. The 2026-10-26 stress-audit update changed the active modulus from 2m to 1.1m psi to align with the published Redwood reference grade; the earlier approximately 0.483% displacement check is now retained separately under docs/history (not rerun at the new modulus). Large cubic/linear mesh convergence and mathematically certified a posteriori error bounds were *not* rerun or proven.

## Upgrading an existing Git clone

1. Copy this folder into a **new** working directory so you have a known backup.
2. Review `docs/MIGRATION.md` and execute the tests and benchmarks.
3. Move the clean tree into a **fresh branch** of your Git clone (or replace tracked contents intentionally); old files will otherwise remain from earlier commits.
4. Inspect `git status` and `git diff --stat` before committing; push only after you have reviewed the deletions and replacements.

## 2026-10-26 revised Redwood and cubic stress scope

Current source-grade modulus is 1.1 million psi. Sample JSON and VTU arrays reflect the new modulus. Cubic 3×3×9 raw longitudinal stress is compared with a distinct beam-theory reference; sectional moment/shear defects and unscreened wood modes are documented rather than hidden. See CUBIC_STRESS_AND_REDWOOD.md.
