"""Common cantilever runner; use python -m tutorials.linear|quadratic|cubic.

Units: in, lbf, psi. Use --nx/--ny/--nz to change the structured grid.
Solves may require considerable memory/time for larger nonlocal stencils.
"""
import argparse
import json
from pathlib import Path
import numpy as np

from core.solvers.beam import build_beam
from core.solvers.loads import build_end_face_traction_load_vector, gravity_load, LOAD_DOF
from core.solvers.constraints import solve_with_face_dirichlet
from core.solvers.results import free_end_uy
from core.reference.awc import DEFAULT_WEIGHT_DENSITY, awc_reference_components
from core.reference.beam_theory import P
from core.materials.design_limits import REDWOOD_BASIS, reference_load_screen

def main(default_formulation="quadratic_rotation_hermite", default_out=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--formulation", choices=(
        "linear_rotation_hermite", "quadratic_rotation_hermite", "cubic_rotation_hermite"),
        default=default_formulation)
    p.add_argument("--nx", type=int, default=4)
    p.add_argument("--ny", type=int, default=5)
    p.add_argument("--nz", type=int, default=17)
    p.add_argument("--weight-density", type=float, default=DEFAULT_WEIGHT_DENSITY)
    p.add_argument("--end-force", type=float, default=P,
                   help="Signed Y end traction resultant in lbf (must be downward and pass the Redwood reference screen)")
    p.add_argument("--section-quadrature-order", type=int, default=4,
                   help="Gaussian order for source-element section equilibrium (cubic; default 4)")
    p.add_argument("--section-axial-positions", type=float, nargs="+", default=[-0.9,0.9],
                   metavar="ZETA", help="Source-section parent coordinates strictly between -1 and 1; default -0.9 +0.9 (inset toward source-element axial ends); explicitly pass 0 only for historical midpoint reproduction")
    p.add_argument("--volume-residual-order", type=int, choices=(2,4),default=2,
                   help="Internal stress-balance sampling: 2 or 4 Gauss points per axis; no centroid diagnostic")
    p.add_argument("--mechanics-audit", action="store_true",
                   help="Cubic: write Beer's mechanics-of-materials signed stress, stress-cube and work-energy diagnostics (not certified)")
    p.add_argument("--shear-warping-audit", action="store_true",
                   help="Cubic: run the source-B/T transverse u_z=y^3 completeness patch test (no solve certification)")
    p.add_argument("--shear-balance", action="store_true",
                   help="Cubic: audit axial shear oscillations, axial averages and FE moment-gradient balance (not certified)")
    p.add_argument("--equilibrium-diagnostics", action="store_true",
                   help="Cubic research: sample div(sigma)+gravity and internal traction jumps; no certified estimate")
    p.add_argument("--save-solution", type=Path, default=None,
                   help="Save the COMPLETE solved DOF state for later stress recovery without re-solving (cubic only)")
    p.add_argument("--section-only", action="store_true",
                   help="Cubic: skip costly visual stress sampling; source-section checks still run, sampled stress screening NOT PERFORMED")
    p.add_argument("--displacement-only-clamp", action="store_true",
                   help="Research BC comparison: enforce displacement trace without fixing nodal rotations")
    p.add_argument("--linear-solver", choices=("scipy","pypardiso"),default="scipy",
                   help="Sparse backend; scipy works without optional PyPardiso")
    p.add_argument("--nodal-stress-audit", action="store_true",
                   help="Cubic: sample one-sided ORIGINAL source element vertices, not centroids; works with --section-only")
    p.add_argument("--stress-audit", action="store_true",
                   help="Write stress comparison and sectional equilibrium report (cubic only; cubic runs default to audit)")
    p.add_argument("--out", type=Path, default=default_out or Path("outputs") / (default_formulation.split("_")[0] + ".json"))
    # Visualization is optional per run; PyVista is a project dependency, but
    # only imported when a viewer or VTU export is requested.
    p.add_argument("--shear-comparison", action="store_true",
                   help="Cubic: open TWO synchronized PyVista shear meshes: classical beam reference versus actual Hermite FE tau_yz; export one-sided corner-stress VTU")
    p.add_argument("--shear-export-only", action="store_true",
                   help="Cubic: write the two-field corner-accurate shear VTU without launching PyVista")
    p.add_argument("--shear-screenshot", type=Path, default=None,
                   help="Cubic: save a two-panel signed shear PNG; implies --shear-comparison and uses offscreen rendering if requested")
    p.add_argument("--view", action="store_true", help="Open the interactive PyVista viewer")
    p.add_argument("--export-vtu", nargs="?", const="auto", metavar="PATH",
                   help="Export sampled Hermite displacement and stress to VTU (optional path)")
    p.add_argument("--screenshot", nargs="?", const="auto", metavar="PATH",
                   help="Render a PNG screenshot (optional path), implies --view")
    p.add_argument("--off-screen", action="store_true", help="Render without a GUI window (useful for screenshots)")
    p.add_argument("--full-screen", action="store_true", help="Launch the interactive viewer full-screen")
    p.add_argument("--edges", action="store_true", help="Show subcell edges")
    p.add_argument("--subdivisions", type=int, default=2,
                   help="Visualization subcells per original element axis (default: 2)")
    p.add_argument("--stress-sample-location", choices=("exact_corners","inset_corners"),
                   default="exact_corners", help="Source-element-sided stress evaluation: exact visualization-subcell corners (default) or inset corners; never centroid")
    p.add_argument("--stress-inset-fraction", type=float, default=0.95,
                   help="When using inset_corners: fraction of corner offset from subcell center (0 < f < 1)")
    p.add_argument("--field", default="Vertical displacement magnitude (in)",
                   choices=("Vertical displacement magnitude (in)", "Displacement magnitude (in)", "Ux (in)", "Uy (in)", "Uz (in)",
                            "Von Mises (psi)", "Absolute longitudinal stress (psi)", "Sigma xx (psi)", "Sigma yy (psi)",
                            "Sigma zz (psi)", "Tau xy (psi)", "Tau yz (psi)",
                            "Tau xz (psi)", "Maximum principal (psi)", "Minimum principal (psi)",
                            "Beam reference sigma zz (psi)","Cubic minus beam sigma zz (psi)",
                            "Beam reference tau yz (psi)","Cubic minus beam tau yz (psi)",
                            "Longitudinal tension utilization","Longitudinal compression utilization",
                            "Shear yz utilization"))
    p.add_argument("--deform-scale", type=float, default=None,
                   help="Displacement magnification (default: automatic 12%% beam-length scale)")
    a = p.parse_args()
    if min(a.nx, a.ny, a.nz) < 3:
        p.error("Recovery needs at least 3 nodes in every direction")
    if a.formulation == "cubic_rotation_hermite" and a.nz < 4:
        p.error("Cubic recovery needs at least 4 nodes in the beam direction")
    if not np.isfinite(a.weight_density) or a.weight_density <= 0:
        p.error("Weight density must be finite and positive")
    if not np.isfinite(a.end_force) or a.end_force > 0:
        p.error("--end-force must be finite and <=0 (downward load convention)")
    if a.formulation != "cubic_rotation_hermite" and a.field in (
            "Beam reference sigma zz (psi)","Cubic minus beam sigma zz (psi)",
            "Beam reference tau yz (psi)","Cubic minus beam tau yz (psi)",
            "Longitudinal tension utilization","Longitudinal compression utilization", "Shear yz utilization"):
        p.error("Cubic stress comparison fields are only available for the cubic formulation")
    if a.section_quadrature_order < 2 or a.section_quadrature_order > 8:
        p.error("--section-quadrature-order must be 2..8")
    if (not all(np.isfinite(z) and -1. < z < 1. for z in a.section_axial_positions)
            or len(set(a.section_axial_positions)) != len(a.section_axial_positions)):
        p.error("--section-axial-positions must be distinct finite values strictly between -1 and 1")
    if a.shear_warping_audit and a.formulation != "cubic_rotation_hermite":
        p.error("--shear-warping-audit is currently cubic-only")
    if a.shear_balance and a.formulation != "cubic_rotation_hermite":
        p.error("--shear-balance is currently cubic-only")
    if a.nodal_stress_audit and a.formulation != "cubic_rotation_hermite":
        p.error('--nodal-stress-audit is currently cubic-only')
    if a.mechanics_audit and a.formulation != "cubic_rotation_hermite":
        p.error("--mechanics-audit is currently cubic-only")
    if a.equilibrium_diagnostics and a.formulation != "cubic_rotation_hermite":
        p.error("--equilibrium-diagnostics is currently cubic-only")
    if a.save_solution is not None and a.formulation != "cubic_rotation_hermite":
        p.error("--save-solution currently supports cubic formulation only")
    if a.section_only and a.formulation != "cubic_rotation_hermite":
        p.error("--section-only is only available for cubic")
    if a.section_only and (a.view or a.screenshot is not None or a.export_vtu is not None
                            or a.off_screen or a.stress_audit):
        p.error("--section-only is incompatible with visual/stress-audit flags")
    if (a.shear_comparison or a.shear_export_only or a.shear_screenshot is not None) and a.formulation != "cubic_rotation_hermite":
        p.error("Shear comparison and export require the cubic formulation")
    if a.section_only and (a.shear_comparison or a.shear_export_only or a.shear_screenshot is not None):
        p.error("--section-only cannot render or export a shear comparison")
    if a.shear_export_only and (a.shear_comparison or a.shear_screenshot is not None):
        p.error("--shear-export-only cannot also request the interactive shear viewer or screenshot")
    if a.stress_audit and a.formulation != "cubic_rotation_hermite":
        p.error("The documented stress-distribution audit is currently cubic-only")
    if not np.isfinite(a.stress_inset_fraction) or not 0 < a.stress_inset_fraction < 1:
        p.error('--stress-inset-fraction must be finite and strictly between 0 and 1')
    if a.subdivisions < 1 or a.subdivisions > 12:
        p.error("--subdivisions must be between 1 and 12")
    if a.deform_scale is not None and (not np.isfinite(a.deform_scale) or a.deform_scale < 0):
        p.error("--deform-scale must be finite and nonnegative")
    basis_check = reference_load_screen(tip_force_lbf=a.end_force,
        density_lbf_in3=a.weight_density)
    if basis_check['status'] != 'PASS_PRELIMINARY':
        p.error('Redwood baseline preliminary load/serviceability screen FAILED; reduce the end load or change dimensions before running')
    beam = build_beam(nx=a.nx, ny=a.ny, nz=a.nz, formulation=a.formulation)
    print(f"Building {a.formulation}, node mesh {a.nx}x{a.ny}x{a.nz}", flush=True)
    transverse_y_cubic = None
    if a.formulation == "cubic_rotation_hermite":
        transverse_y_cubic = (0,3,0) in {tuple(t) for t in beam.cubic_terms}
        if not transverse_y_cubic:
            print("WARNING: cubic y^3 warping mode is not recoverable with only three "
                  "transverse y node layers; nodal shear may be inaccurate. "
                  "Try --ny 4 or more and verify convergence.", flush=True)
    end, _, _ = build_end_face_traction_load_vector(beam, total_load=a.end_force)
    weight = gravity_load(beam, a.weight_density)
    total = end + weight
    fy = float(total[beam.global_dofs[:, LOAD_DOF]].sum())
    expected = float(a.end_force - a.weight_density * beam.Lx * beam.Ly * beam.Lz)
    if not np.isclose(fy, expected, rtol=1e-8, atol=1e-6):
        raise RuntimeError(f"Force mismatch: {fy} vs {expected}")
    print(f"Applied y force {fy:.6f} lbf; assembling stiffness...", flush=True)
    U, solve_information = solve_with_face_dirichlet(
        beam, total, lambda x: (0., 0., 0.),
        None if a.displacement_only_clamp else lambda x: (0., 0., 0.),
        constrain_rotations=not a.displacement_only_clamp,
        linear_solver=a.linear_solver)
    tip = float(free_end_uy(beam, U))
    reference = awc_reference_components(a.weight_density, end_force=a.end_force)
    error = 100.0 * abs(tip - reference["total_uy_in"]) / abs(reference["total_uy_in"])
    assert beam.orthotropic_constants is not None
    result = {
        "formulation": a.formulation, "research_not_certified": True,
        "comparison_is_not_a_posteriori_bound": True,
        "recovery": "linear" if a.formulation.startswith("linear") else ("quadratic" if a.formulation.startswith("quadratic") else "cubic"),
        "node_mesh": [a.nx, a.ny, a.nz],
        "fem_tip_uy_in": tip, "awc_reference": reference,
        "difference_vs_awc_pct": error, "under_one_percent": bool(error < 1.0),
        "total_y_force_lbf": fy, "calibrated_modulus": False,
        "longitudinal_E_psi": float(beam.orthotropic_constants.Ez),
        "end_force_lbf":float(a.end_force),
        "redwood_reference_grade":REDWOOD_BASIS.grade,
        "redwood_reference_screen":basis_check,
        "stress_distribution_certified":False,
        "boundary_condition": ("strong KKT face-displacement trace; nodal rotations free"
            if a.displacement_only_clamp else
            "strong KKT face-displacement trace and fixed nodal rotations"),
        "linear_solver_diagnostics":solve_information,
    }
    if a.save_solution is not None:
        from core.solvers.solution_archive import save_cubic_solution
        save_cubic_solution(a.save_solution, beam, U, end_force=a.end_force,
                            weight_density=a.weight_density,
                            displacement_only_clamp=a.displacement_only_clamp,
                            fem_tip_uy_in=tip)
        result['solution_state_file'] = str(a.save_solution)
        print(f"Saved complete solved displacement state: {a.save_solution}", flush=True)
    if a.formulation == "cubic_rotation_hermite":
        result["cubic_transverse_y_cubic_recovery_present"] = transverse_y_cubic
        result["cubic_y_warping_warning"] = (None if transverse_y_cubic else
            "Only three node layers in y; exact cubic y^3 displacement warping not represented")
    if a.shear_warping_audit:
        from core.verification.shear_warping import audit_cubic_y_warping
        patch = audit_cubic_y_warping(beam)
        warp_path = a.out.with_name(a.out.stem + "_shear_warping_patch.json")
        warp_path.parent.mkdir(parents=True, exist_ok=True)
        warp_path.write_text(json.dumps(patch,indent=2) + "\n", encoding="utf8")
        result["cubic_shear_warping_patch_file"] = str(warp_path)
        result["cubic_shear_warping_patch_relative_strain_L2"] = patch["relative_discrete_strain_L2"]
        print(f"Cubic y^3 kinematic warping patch relative strain L2: "
              f"{patch['relative_discrete_strain_L2']:.6g} "
              f"(kinematic test only)",flush=True)
    if a.formulation == "cubic_rotation_hermite":
        from core.solvers.section_equilibrium import integrate_sections
        source_report=integrate_sections(
            beam,U,order=a.section_quadrature_order,
            axial_positions=a.section_axial_positions,
            end_force=a.end_force,density=a.weight_density)
        section_path=a.out.with_name(a.out.stem+"_source_sections.json")
        section_path.parent.mkdir(parents=True,exist_ok=True)
        section_path.write_text(json.dumps(source_report,indent=2)+"\n",encoding="utf-8")
        result["cubic_source_section_file"]=str(section_path)
        result["cubic_source_max_midspan_shear_y_difference_pct"]=source_report["max_midspan_shear_y_difference_pct"]
        result["cubic_source_max_midspan_moment_x_difference_pct"]=source_report["max_midspan_moment_x_difference_pct"]
        result["cubic_source_section_quadrature_order"]=a.section_quadrature_order
        result["sampled_cubic_stress_screen_performed"]=not a.section_only
        print(f"Source-Gauss section check: shear discrepancy "
              f"{source_report['max_midspan_shear_y_difference_pct']:.3f}%, "
              f"moment discrepancy {source_report['max_midspan_moment_x_difference_pct']:.3f}% "
              f"(uncertified beam reference)",flush=True)
    if a.shear_balance:
        from core.solvers.shear_balance import audit_axial_shear_balance
        shear_balance=audit_axial_shear_balance(
            beam,U,order=4,section_order=a.section_quadrature_order,
            end_force=a.end_force,density=a.weight_density)
        path=a.out.with_name(a.out.stem+'_shear_balance.json')
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(shear_balance,indent=2)+'\n',encoding='utf-8')
        result['cubic_axial_shear_balance_file']=str(path)
        result['cubic_max_interior_sampled_shear_difference_pct']=shear_balance['max_interior_sampled_shear_reference_difference_pct']
        result['cubic_max_interior_axially_averaged_shear_difference_pct']=shear_balance['max_interior_axially_averaged_shear_reference_difference_pct']
        result['cubic_max_interior_shear_moment_balance_difference_pct']=shear_balance['max_interior_sampled_shear_moment_balance_difference_pct']
        print(f"Axial shear audit: worst interior sampled "
              f"{result['cubic_max_interior_sampled_shear_difference_pct']:.3f}%, "
              f"cell-averaged {result['cubic_max_interior_axially_averaged_shear_difference_pct']:.3f}%, "
              f"FE dM/dz balance {result['cubic_max_interior_shear_moment_balance_difference_pct']:.3f}% "
              f"(all diagnostic, not certified)",flush=True)
    if a.nodal_stress_audit or a.stress_audit:
        from core.verification.nodal_stress import audit_source_vertex_stress
        nodal=audit_source_vertex_stress(beam,U,location=a.stress_sample_location,
            inset_fraction=a.stress_inset_fraction,
            end_force=a.end_force,density=a.weight_density)
        path=a.out.with_name(a.out.stem+'_source_vertex_stress.json')
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(nodal,indent=2)+'\n',encoding='utf8')
        result['source_vertex_stress_file']=str(path)
        result['source_vertex_sigma_zz_interior_L2_pct']=nodal['interior_sigma_zz_discrete_relative_L2_pct']
        result['source_vertex_tau_yz_interior_L2_pct']=nodal['interior_tau_yz_discrete_relative_L2_pct']
        print(f"Original-source-vertex stress: longitudinal L2 {result['source_vertex_sigma_zz_interior_L2_pct']:.3f}%, "
              f"shear L2 {result['source_vertex_tau_yz_interior_L2_pct']:.3f}% "
              f"(discrete, not certified)",flush=True)
    if a.mechanics_audit:
        from core.verification.solved_audit import solved_mechanics_report
        mechanics=solved_mechanics_report(beam,U,total,end_force=a.end_force,
                                           weight_density=a.weight_density)
        mechanics_path=a.out.with_name(a.out.stem+'_mechanics_audit.json')
        mechanics_path.parent.mkdir(parents=True,exist_ok=True)
        mechanics_path.write_text(json.dumps(mechanics,indent=2)+'\n',encoding='utf8')
        result['mechanics_audit_file']=str(mechanics_path)
        result['mechanics_clapeyron_relative_residual']=mechanics['clapeyron_relative_residual']
        print(f"Source-FE Clapeyron work/energy relative residual: "
              f"{mechanics['clapeyron_relative_residual']:.6g} (not a 3D error bound)",flush=True)
    if a.equilibrium_diagnostics:
        from core.solvers.strong_residual import sample_volume_residual
        from core.solvers.traction_jumps import source_traction_jumps
        volume_balance=sample_volume_residual(beam,U,density=a.weight_density,gauss_order=a.volume_residual_order)
        interface_balance=source_traction_jumps(beam,U,order=a.section_quadrature_order)
        equilibrium_path=a.out.with_name(a.out.stem+"_equilibrium_diagnostics.json")
        equilibrium_path.parent.mkdir(parents=True,exist_ok=True)
        equilibrium_path.write_text(json.dumps({
            'status':'SAMPLED_EQUILIBRIUM_DIAGNOSTICS_NOT_CERTIFIED',
            'certified_error_bound':False,
            'volume_residual':volume_balance,
            'interface_traction_jumps':interface_balance,
        },indent=2)+"\n",encoding="utf-8")
        result['cubic_equilibrium_diagnostics_file']=str(equilibrium_path)
        result['volume_balance_rms_lbf_in3']=volume_balance['weighted_sample_rms_div_sigma_plus_b_lbf_in3']
        print(f"Raw-stress equilibrium diagnostic: volume residual RMS "
              f"{result['volume_balance_rms_lbf_in3']:.6g} lbf/in^3; "
              f"not an error bound",flush=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Tip uy={tip:.9f} in; AWC={reference['total_uy_in']:.9f} in")
    print(f"Difference={error:.6f}% ; under 1%: {error < 1.0}")
    print(f"Wrote {a.out}")
    wants_view = a.view or (a.screenshot is not None) or (a.off_screen and not (a.shear_comparison or a.shear_export_only or a.shear_screenshot is not None))
    is_cubic=a.formulation=="cubic_rotation_hermite"
    wants_shear_comparison = a.shear_comparison or a.shear_screenshot is not None
    wants_shear_export = wants_shear_comparison or a.shear_export_only
    if wants_view or wants_shear_export or a.export_vtu is not None or (is_cubic and not a.section_only) or a.stress_audit:
        from core.visualization.sampling import sample_volume
        from core.visualization.pyvista_viewer import export_vtu, show_sample
        print(f"Sampling recovered Hermite field, {a.subdivisions} subcell(s) per axis...", flush=True)
        sampled = sample_volume(beam, U, subdivisions=a.subdivisions,
                                stress_sample_location=a.stress_sample_location,
                                stress_inset_fraction=a.stress_inset_fraction)
        stem = a.formulation.split("_")[0]
        if is_cubic:
            from core.visualization.cubic_stress import cubic_stress_audit, write_stress_audit
            audit=cubic_stress_audit(sampled,beam,end_force=a.end_force,density=a.weight_density)
            audit_path=a.out.with_name(a.out.stem+'_stress_audit.json')
            write_stress_audit(audit,audit_path)
            result['cubic_stress_audit_file']=str(audit_path)
            result['cubic_stress_screen_status']=audit['status']
            result['quarter_to_threequarter_span_sigma_zz_relative_L2_pct']=audit['quarter_to_threequarter_span_sigma_zz_relative_L2_pct']
            result['stress_sample_location']=audit['stress_sample_location']
            result['stress_sample_count']=audit['raw_corner_sample_count']
            a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
            print(f"Cubic {a.stress_sample_location} longitudinal-stress interior beam-L2 difference: "
                  f"{audit['quarter_to_threequarter_span_sigma_zz_relative_L2_pct']:.3f}%",flush=True)
            print(f"Stress screening: {audit['status']} (research only, not 3D stress certification)",flush=True)
            if audit['status'] != 'PASS_PRELIMINARY_SAMPLE_SCREEN':
                raise RuntimeError('Cubic sampled stress exceeds stated provisional Redwood design screen; outputs are research only')
        if a.export_vtu is not None:
            export_path = (Path("outputs") / f"{stem}.vtu" if a.export_vtu == "auto"
                           else Path(a.export_vtu))
            export_vtu(sampled, export_path)
            print(f"Exported VTU: {export_path}", flush=True)
        if wants_shear_export:
            from core.visualization.shear_comparison import (
                prepare_shear_comparison, export_shear_vtu, show_shear_grid, as_shear_pyvista_grid)
            # The paired shear mesh is the ORIGINAL source FE Hex8 topology,
            # sampled at all original element corners (one-sided): NOT the
            # extra visualization subcells used by the general stress viewer.
            shear_source_sample = sample_volume(
                beam, U, subdivisions=1,
                stress_sample_location=a.stress_sample_location,
                stress_inset_fraction=a.stress_inset_fraction)
            data = prepare_shear_comparison(shear_source_sample, beam,
                                            end_force=a.end_force,
                                            density=a.weight_density)
            shear_vtu = a.out.with_name(a.out.stem + '_shear_comparison.vtu')
            export_shear_vtu(data, shear_vtu)
            result['cubic_shear_comparison_vtu'] = str(shear_vtu)
            result['cubic_shear_comparison_stress_sources'] = {
                'beam_reference': 'classical approximate rectangular-section Saint-Venant tau_yz',
                'source_fe': 'one-sided source cubic Hermite D @ B @ (T @ U) stress at subcell corners',
                'certified_3d_error_bound': False,
            }
            a.out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
            print(f'Exported independent corner-valued shear fields: {shear_vtu}', flush=True)
            shot = (a.shear_screenshot if a.shear_screenshot is not None else
                    a.out.with_name(a.out.stem + '_shear_comparison.png') if a.off_screen else None)
            if wants_shear_comparison:
                show_shear_grid(as_shear_pyvista_grid(data), title=a.formulation,
                                off_screen=a.off_screen, full_screen=a.full_screen,
                                screenshot=shot, show_edges=a.edges)
                if shot is not None:
                    print(f'Exported dual-shear PNG: {shot}', flush=True)
        if wants_view:
            shot = (None if a.screenshot is None else
                    Path("outputs") / f"{stem}.png" if a.screenshot == "auto" else Path(a.screenshot))
            show_sample(sampled, formulation=a.formulation,
                        field=a.field, deform_scale=a.deform_scale,
                        show_edges=a.edges, full_screen=a.full_screen,
                        off_screen=a.off_screen, screenshot=shot)
    return result

if __name__ == "__main__":
    main()
