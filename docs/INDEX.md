# Documentation index

**OOP Hex27/Hex64 refactor:** [Tutorial](OOP_TUTORIAL.md) · [TODO / acceptance checklist](OOP_TODO.md) · [Build and test report](BUILD_REPORT_OOP.md)

- [29 October 2026 project audit](AUDIT_20261029.md): documented link/CLI regression checks, data-derived plot and scope of unsolved work.

## Project documentation (the sole project README is ../README.md)

- [Hex27 and Hex64 higher-node p-study](HIGHER_NODE_HEX27_HEX64_STUDY_20261027.md): latest matched 3×5×25 results, available two-mesh plot, commands, diagnostics, reproducibility limits
- [Architecture](ARCHITECTURE.md): code ownership and boundary-condition choices
- [Formulations](FORMULATIONS.md): exactly what linear/quadratic/cubic recovery means
- [Commands](COMMANDS.md): sequential Windows PowerShell commands
- [Example results](EXAMPLE_RESULTS.md): small-grid baseline and reproduced quadratic result
- [Mechanics of Materials verification](MECHANICS_OF_MATERIALS_VERIFICATION.md): independent stress cube, bending, shear and energy checks (research only)
- [Two-view PyVista shear tutorial](DUAL_SHEAR_PYVISTA.md): signed classical beam-reference and raw Hermite shear fields with sampling and interpretation limits
- [Validation](VALIDATION.md): tested checks and missing proof obligations
- [Linear research plan](LINEAR_RESEARCH_PLAN.md): systematically diagnose excessive stiffness
- [Cubic refinement and a posteriori roadmap](CUBIC_REFINEMENT_AND_APOSTERIORI.md): controlled mesh-density experiment and the missing certification steps
- [Migration](MIGRATION.md): included/excluded source and updated paths
- [Certification status](CERTIFICATION_STATUS.json): explicit machine-readable *not certified* state
- [Active source pin](ACTIVE_SOURCE_SHA256.json): distribution identity for moved source
- [Original source pins](UPSTREAM_SOURCE_SHA256.json): evidence from supplied archive, not active checks

- [Back to project README](../README.md): primary getting-started guide and all result tables

- [Release checklist](RELEASE_CHECKLIST.md): manual Windows PyVista acceptance and release gates
- [Release audit](RELEASE_AUDIT.md): documented corrections and outstanding review gates

- [New October 2026 source-element section equilibrium study](SECTION_EQUILIBRIUM_REFINEMENT_20261026.md) — direct Gauss resultants and sampled 3D balance diagnostics; not certified.
- [Higher-density cubic study](DENSIFICATION_STUDY_20261026.md) — plotted five-level 5×5 axial series plus separate 7×7×17 transverse level, exact measured JSON input and plot-regeneration command, and section-location sensitivity.

- [2026-10-27 sectional-shear oscillation and moment-balance investigation](SHEAR_BALANCE_INVESTIGATION_20261027.md): raw FE stress, refined sampling and limits of midpoint agreement.

- [2026-10-27 source-corner stress audit](CORNER_STRESS_EVALUATION_20261027.md): one-sided original nodal stresses, inset option, changed sample norms, and integration safeguards.
- [Transverse cubic shear recovery and original-nodal benchmark (27 Oct 2026)](TRANSVERSE_CUBIC_SHEAR_RECOVERY_20261027.md): cubic warping completeness, the 3×5 axial-refinement series through 49 nodes, signed beam-reference comparisons, and limitations.

**Publication readiness:** [29 October independent-style code and claims audit](PUBLICATION_AUDIT_20261029.md) (open release gates and evidence limits).
