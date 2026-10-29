# Mandatory policy: never calculate/report stresses at element centroids

**Status: required for all new runs, exports, viewers, reports and examples.**

1. **Do not evaluate finite-element stress at element or visualization-subcell centroids** for stress maps, extrema, comparisons, or reports. Do not substitute a single source-cell midpoint stress for an entire element or its fine display cells. A centroid may still be used as a geometric quantity in a moment formula, or for locating a section; it is **not** a stress-sampling site.
2. Evaluate the true solver strain operator `B(x)` and constitutive law `sigma(x) = D @ B(x) @ U` at **one-sided source-element corners** or at **inset corners** near the eight extremities of every *visualization subcell*. No theoretical beam stress may be mislabeled as FE stress. For the nodal Hex27/Hex64 workflow the standard inset fraction is **0.95 of each subcell half-width**, giving eight locations per subcell and exactly `8 * (nx-1)*(ny-1)*(nz-1) * p^3` raw stress samples.
3. Preserve distinct one-sided FE stresses at element interfaces. Do **not** smooth/average stresses between neighboring source elements. The displayed nodal Hex27/Hex64 subcell colors are the *most extreme directly sampled* corner values (maximum or most negative, depending on the field), **not** stress evaluations at cell centers, not continuous interpolation, and not guaranteed 3D maxima.
4. **Small meshes with at most 16 source elements are SMOKE TESTS ONLY.** Show this warning prominently in README sections, command examples, screenshot overlays, summary JSON and report plots. Report sampled stress extrema only as diagnostics: **NOT VALIDATED FOR STRESS; NOT FOR DESIGN.** A dense stress-convergence study is still required.
5. Reject outdated NPZ files that use single source-center stress evaluations. `tools.nodal_field_viewer` and `tools.convergence_study --stage view` require archive scheme `inset_subcell_corners_v2`. Regenerate from the saved numerical solve command before viewing stress. Old images and their stress numbers are **superseded**, not compatible with this sampling convention.
6. Numerical integration for element assembly and accurate cross-sectional resultants **must still use the appropriate Gaussian quadrature**. This rule does not ban Gaussian points or midspan cross-sectional locations; it bans substituting centroid-only values for stress results. Never trade a valid integration rule for a corner-only quadrature approximation.
7. The prior recovered-rotation Hermite viewer uses source-element subcell corners (inset optional) and keeps the raw eight one-sided stress samples for audit. Historical figures comparing subcell-center samples are explicitly archival and do not authorize their future use.

## Reproduce the corrected 16-element examples (research smoke test only)

```powershell
python -m tools.p_order_study --mode solve --order 2 --nx 3 --ny 3 --nz 5 --field-output outputs/hex27_3x3x5_fields.npz
python -m tools.p_order_study --mode solve --order 3 --nx 3 --ny 3 --nz 5 --field-output outputs/hex64_3x3x5_fields.npz
python -m tools.nodal_field_viewer outputs/hex64_3x3x5_fields.npz --field "Von Mises (psi)"
python -m unittest tests.test_nodal_field_workflow -v
```

**Warning:** The numerical solution is unchanged by the archive exporter modification; the stress sampling locations and displayed sampled maxima **have changed**. The prior ~52.6 psi color-scale maximum was from deprecated source-center sampling; it is not a usable maximum-stress benchmark. Current sample extrema likewise cannot establish a physical maximum or a certified bound.

## Restartable dense cubic Hermite shear fields

`tools/replay_cubic_fields.py` reads full DOFs created by
`python -m tutorials.cubic --save-solution ...` and exports **separate VTK
point arrays for each one-sided 95%-inset corner stress**. It never reads
centroid stress and does not take extrema from cell-averaged fields. To
compare FE shear against gravity-inclusive 3.732 psi classical shear, use
matching sample coordinates and the saved JSON report; see
[`CUBIC_129_LOCAL_STRESS_WORKFLOW.md`](CUBIC_129_LOCAL_STRESS_WORKFLOW.md).

## Strong-form residual diagnostics

A final publication audit also found a legacy *stress derivative* diagnostic
that could sample at the element centroid using a one-point Gauss rule. This
option has been retired. `--volume-residual-order` accepts `2` or `4`;
`point_residual` refuses an explicit `(0,0,0)` parent coordinate and otherwise
uses non-centroid interior points. This is a research strong-form balance
check, not the assembly/section quadrature and not a certified stress norm.
