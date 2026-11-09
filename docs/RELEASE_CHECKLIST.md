# Public release checklist (updated 2026-10-29)

This is a **research prototype and programming portfolio**, not a calibrated
or certified structural engineering package. Run this checklist on the final
published ZIP and commit, *not* a different working tree.

## Automated (run from project root)

```powershell
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python -m unittest tests.test_mechanics_benchmarks -v
python -m tools.mechanics_benchmarks
python -m tests.checks.element_modes
python -m tests.checks.source_provenance
python -m tests.checks.certification_readiness
```

Expected: core unit checks pass; certification readiness reports **NOT
CERTIFIED**. The last command exits normally as a *status report*; to test
fail-closed behavior use `--require-certified`, which should exit nonzero.

## Hands-on PyVista acceptance (Windows desktop, not verified here)

```powershell
python -m core.visualization.pyvista_viewer samples/cubic_3x3x9.vtu
python -m core.visualization.pyvista_viewer samples/cubic_3x3x9.vtu --field "Sigma zz (psi)"
python -m core.visualization.pyvista_viewer samples/cubic_3x3x9.vtu --field "Beam reference sigma zz (psi)"
python -m core.visualization.pyvista_viewer samples/cubic_3x3x9.vtu --field "Absolute longitudinal stress (psi)"
```

- [ ] All commands open a nonblank model without a traceback.
- [ ] Default `|Uy|` is near **zero (blue)** at fixed `z=0`, greatest
  **(red)** at free `z=L`.
- [ ] Nine yellow source-node cubes mark fixed `z=0`; support and free-end
  labels remain attached to the correct end when rotating.
- [ ] `6` signed σzz has positive/red and negative/blue surfaces; `A` has
  blue=least and red=greatest **absolute** longitudinal stress.
- [ ] `0`, `B`, `D`, `M`, `C`, `E`, `+`, `-`, `S`, `R` work without disappearing
  geometry, missing actors or unreadable color legends.
- [ ] Scalar bar and HUD have readable 10-point white text against dark
  background. Deformation scale is visibly **display only**.
- [ ] `D` or `--deform-scale 1` makes deformation magnitude interpretation
  consistent with the reported ~0.116-in free-end deflection.
- [ ] Repeat with linear and quadratic precomputed VTUs.
- [ ] Fresh cubic 3×3×9 solve produces agreement with shipped JSON within
  meaningful numerical tolerance; check no PyVista/VTK warnings.
- [ ] Confirm one root README, no `.pyc` / `__pycache__`, no stale results,
  and no unsupported structural certification promises.

## Publishing and licensing decision

- [ ] Author chooses and adds a software license (`LICENSE`) appropriate to
  intended source reuse. Public visibility on GitHub is **not** a reuse license.
- [ ] Public repository points to this exact tested release, and includes
  documented dependency versions and a self-contained demonstration.
- [ ] Do not promise engineering software certification, graded-material
  approval, or licensed engineering services based solely on these examples.

## Evidence and limitations

Automated tests establish source-level, example-level and VTU numerical
consistency; **they do not guarantee cross-version GUI compatibility**.
The current environment lacked PyVista and could not perform a graphical
window smoke test. The historical 4×5×17 quadratic benchmark uses a different
2,000,000-psi modulus and is not a revalidated graded material result.

The grade-value material sheet is a Redwood *design reference*, not an
ultimate failure-strength dataset. The 4×8×120-in beam dimensions represent
an idealized solid, not necessarily dressed or grade-stamped stock.
Do not describe the 2.16% stress comparison as an a posteriori upper bound:
it compares interior **σzz** against elementary beam theory only and is
sensitive to sampling; the **shear field is not equivalently validated**.

**Publication gate:** do not call the release "GUI tested" or "fully verified"
until the Windows graphical checks above have been recorded as passed.

## Additional prepublication stress/reproducibility gates (29 October)

- [ ] Run a full 5×5×129 **new** solve and save the cubic v2 displacement checkpoint.
- [ ] Compare *local* FEM $\tau_{yz}$ with beam theory at identical physical inset-corner positions across three or more meshes, avoiding the immediate clamp region for a beam-theory comparison.
- [ ] Verify results at fixed physical positions, not simply the maximum of mesh-dependent sample sets.
- [ ] Independently check numerical conditioning and error sensitivity of the historic dense KKT solve: the recorded RHS-relative residual is about 0.274%, despite its tiny backward error.
- [ ] Verify high-density stress extrema under refinement and establish whether clamp effects are singular, mesh-dependent, or physical.
- [ ] Run the VS Code tasks and PyVista on target Windows/Python builds, including the dual-shear viewer; no such desktop GUI run was completed here.
- [ ] Select a software license if reuse by third parties is intended.

See [publication audit](PUBLICATION_AUDIT_20261029.md). A successful unit suite
**does not** mark the unchecked engineering evidence as passed.
