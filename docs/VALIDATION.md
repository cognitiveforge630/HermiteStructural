# Verification and certification status

## Current checks

1. `python -m unittest discover -s tests -v` — active element construction, dimensions, **rigid-body and affine strain patch tests** across all three formulations, structured mesh, prescribed orthotropic material, independent compliance formulas, traction and body force integrated resultants, and invalid-runner parameters. **Fast.**
2. `python -m tests.checks.element_modes` — checks finite recovered operators, positive element Jacobian, recovery support, and rank indicators for all three methods. **No global solve.**
3. `python -m tests.checks.source_provenance` — compares normalized source bytes to `docs/ACTIVE_SOURCE_SHA256.json`. **Release source identity only**, NOT a correctness or numerical accuracy proof.
4. `python -m tests.checks.certification_readiness` — machine-readable fail-closed certification-gate checklist. **Diagnostic, not a certificate.** Add `--require-certified` to deliberately return an error until all independent gates are satisfied.
5. `python -m tutorials.<order>` — optional full solve plus signed beam-theory displacement comparison. **Cost depends strongly on mesh and stencil.**

6. `python -m unittest tests.test_visualization -v` — independent no-GUI recovered-displacement and stress-sampling unit tests across all three formulations. PyVista GUI availability and graphics driver behavior must also be checked on the target desktop.

## Invariants to retain in future developments

* Eight-vertex geometry/connectivity and stable nodal ordering across variants.
* Physical nodal rotations and the nonlocal `element_transform` correctly included in stiffness and loads.
* End traction and self-weight integrated using the same element interpolation used in assembly.
* `E_z=1,100,000 psi` from the published Redwood reference grade; prescribed research orthotropic ratios are not adjusted to fit a displacement. The previous E=2m result remains historical only.
* Consistent definition of global `uy` displacement and of volume/body-force density.
* Independently checked KKT constraints and treatment of rotations.
* No false mathematical bound claim; distinguish mesh-difference QoI, beam comparison, and a posteriori estimator.

## Missing for rigorous a posteriori certification

A global, admissible displacement reconstruction with verified boundary trace; valid and equilibrated admissible stresses on every relevant surface; material and load assumptions; stable numerical integration and rounding/linear-solve enclosures; mathematically justified energy/QoI bound transfer; separate mesh convergence and singularity analysis; independent signoff where needed. No current check passes all these gates.

**Upstream source note.** `UPSTREAM_SOURCE_SHA256.json` preserves the supplied archive's original manifest only as migration evidence. Paths and imports intentionally changed when active modules were moved; the active identity snapshot is `ACTIVE_SOURCE_SHA256.json`. Comparing either hash snapshot does not certify the FEM model.

## Cubic stress checks

`python -m unittest tests.test_cubic_stress -v` verifies grade basis, load failure screening, signed beam stress, preservation of unmodified FE stress, benchmark interior bending disagreement and the non-certified designation. Cubic tutorial additionally writes raw FE stresses, separate analytical stress references and a section moment/shear report. These are research-only validation indicators; the FE shear distribution remains less accurate.

## Additional Mechanics of Materials verification (October 2026)

- `python -m unittest tests.test_mechanics_benchmarks -v` — analytical
  point/UDL bending references; sigma/tau signs; stress-cube 45-degree
  transform, traction and principal stress; energy and shear-strain factors;
  prescribed cubic source-element exact tension, compression, shear and bending.
- `python -m tools.mechanics_benchmarks` — nine fail-closed quick research
  checks, written to `outputs/mechanics_benchmarks.json`. A passing result
  deliberately reads `PASS_RESEARCH_BENCHMARKS_NOT_CERTIFIED`.
- `python -m tutorials.cubic --nx 3 --ny 3 --nz 5 --section-only --mechanics-audit`
  — optional *solved* FE displacement/stress sample and Clapeyron work-energy
  comparison. Extra source quadrature; no GUI or a posteriori claim.

The verification report documents methods, conventions and sampled values:
[Mechanics of Materials verification](MECHANICS_OF_MATERIALS_VERIFICATION.md).
