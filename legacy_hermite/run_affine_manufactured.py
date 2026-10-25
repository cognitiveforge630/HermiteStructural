"""Exactly solvable constant-strain axial traction benchmark (research only).

Prescribe u_z = strain_zz * z, u_x=u_y=0, clamped z=0,
and exact constant boundary traction t = (D*epsilon) @ n on every other
surface.  For homogeneous orthotropic Cauchy elasticity this is an exact
continuum solution.  Independently integrate the TRUE energy-norm stress
error; compare against the UNVALIDATED residual indicator.

Unlike an end-loaded bending cantilever this has an analytic exact 3D solution.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from aposteriori.residual import FACES, compute_residual_indicator, face_parent_point, face_geometry, voigt_to_stress
from solver import build_beam, solve_with_force_vector


def exact_constant_strain_load(beam, stress):
    """Integrate consistent exact outward traction across *all* natural faces."""
    n_dof = beam.grid_points.shape[0] * beam.ndof_per_node
    force = np.zeros(n_dof)
    face_counts = Counter(tuple(sorted(int(nodes[n]) for n in local))
                          for nodes in beam.elements for _, _, local in FACES)
    for e, nodes in enumerate(beam.elements):
        coords = beam.get_element_nodes(e)
        dofs = beam.get_element_global_dofs(e).ravel()
        for axis, side, local in FACES:
            key = tuple(sorted(int(nodes[n]) for n in local))
            if face_counts[key] != 1:
                continue
            if axis == 2 and side == -1 and np.allclose(coords[:4, 2], 0):
                continue  # essential boundary; displacement u=0 there
            for i, s in enumerate(beam.gauss_points):
                for j, t in enumerate(beam.gauss_points):
                    parent = face_parent_point(axis, side, s, t)
                    normal, area = face_geometry(beam, coords, parent, axis, side)
                    N_disp, _, _, _ = beam.get_hermite_displacement_matrices(*parent, coords)
                    force[dofs] += (N_disp.T @ (stress @ normal) * area *
                                    beam.gauss_weights[i] * beam.gauss_weights[j])
    return force


def true_energy_error(beam, U, exact_e):
    """Independently integrate ||epsilon(U)-exact epsilon||_D / ||exact epsilon||_D."""
    err_sq, exact_sq, field_sq = 0., 0., 0.
    for e in range(beam.get_num_elements()):
        coords = beam.get_element_nodes(e)
        ue = U[beam.get_element_global_dofs(e).ravel()]
        for i, xi in enumerate(beam.gauss_points):
            for j, eta in enumerate(beam.gauss_points):
                for k, zeta in enumerate(beam.gauss_points):
                    B, detJ, _, _ = beam.build_B_matrix(xi, eta, zeta, coords)
                    weight = detJ * beam.gauss_weights[i] * beam.gauss_weights[j] * beam.gauss_weights[k]
                    calculated = B @ ue
                    difference = calculated - exact_e
                    err_sq += float(difference @ beam.D @ difference) * weight
                    exact_sq += float(exact_e @ beam.D @ exact_e) * weight
                    field_sq += float(calculated @ beam.D @ calculated) * weight
    return {
        'true_energy_error_sqrt_lbf_in': float(np.sqrt(max(0., err_sq))),
        'exact_energy_norm_sqrt_lbf_in': float(np.sqrt(max(0., exact_sq))),
        'computed_energy_norm_sqrt_lbf_in': float(np.sqrt(max(0., field_sq))),
        'true_relative_energy_error_percent': (float(100 * np.sqrt(err_sq / exact_sq))
                                               if exact_sq > 0 else None)
    }


def evaluate(formulation='cauchy_hex8_reference', mesh=(1, 1, 2), strain_zz=1e-5,
             length_scale='diameter'):
    beam = build_beam(nx=mesh[0]+1, ny=mesh[1]+1, nz=mesh[2]+1,
                      formulation=formulation)
    exact_e = np.array([0., 0., strain_zz, 0., 0., 0.])
    exact_stress = voigt_to_stress(beam.D @ exact_e)
    F = exact_constant_strain_load(beam, exact_stress)
    U = solve_with_force_vector(beam, F)
    true = true_energy_error(beam, U, exact_e)
    residual = compute_residual_indicator(beam, U, boundary_stress=exact_stress,
                                           length_scale=length_scale)
    exact_U = np.zeros_like(U).reshape((-1,beam.ndof_per_node))
    exact_U[:,2] = beam.grid_points[:,2] * strain_zz
    nodal_error = float(np.max(np.abs(U-exact_U.ravel())))
    norm_true = true['true_energy_error_sqrt_lbf_in']
    norm_est = residual['eta_energy_units_sqrt_lbf_in']
    effectivity = (norm_est / norm_true if norm_true > 1e-8 * true['exact_energy_norm_sqrt_lbf_in'] else None)
    return {'formulation': formulation, 'mesh_elements_xyz': mesh,
            'strain_zz': strain_zz, 'length_scale': length_scale,
            'max_nodal_dof_error_in_or_rad': nodal_error,
            'exact_stress_psi':exact_stress.tolist(),
            'nonzero_load_dofs': int(np.count_nonzero(np.abs(F)>1e-13)),
            **true,
            'experimental_eta_percent': residual['estimated_energy_error_percent_experimental'],
            'experimental_eta_sqrt_lbf_in': norm_est,
            'experimental_effectivity_index': effectivity,
            'volume_sq_lbf_in': residual['volume_residual_squared_lbf_in'],
            'jump_sq_lbf_in': residual['interior_traction_jump_squared_lbf_in'],
            'boundary_sq_lbf_in': residual['natural_traction_mismatch_squared_lbf_in'],
            'status': 'EXACT_SOLUTION_TEST_ESTIMATOR_NOT_CERTIFIED'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', nargs=3, type=int, metavar=('EX','EY','EZ'), default=(1,1,2))
    parser.add_argument('--strain',type=float,default=1e-5)
    parser.add_argument('--formulations',nargs='*',default=['legacy_hermite','cauchy_hex8_reference','affine_enriched_hex8'])
    parser.add_argument('--out',type=Path,default=Path('outputs/manufactured_affine'))
    args=parser.parse_args()
    if min(args.mesh)<1 or args.strain==0: parser.error('Positive mesh sizes and nonzero strain required')
    rows=[]
    for formulation in args.formulations:
        for scale in ('diameter','geometric_mean','min_edge'):
            print(f'[{formulation}/{scale}] computing...',flush=True)
            result=evaluate(formulation,tuple(args.mesh),args.strain,scale)
            rows.append(result)
            print(f'  true energy error={result["true_relative_energy_error_percent"]:.6g}% '
                  f' experimental eta={result["experimental_eta_percent"]:.6g}%',flush=True)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'results.json').write_text(json.dumps({'rows':rows},indent=2)+'\n',encoding='utf-8')
    report=['# Constant-strain manufactured-solution verification','',
            'Exact displacement: u_z = strain_zz*z; all other components 0.','',
            'Exact normal traction is prescribed on all natural faces; z=0 fixed.',
            'Independent exact *energy-norm* error is shown beside the experimental h-scaled estimator.',
            '', '| Formulation | Element scale | TRUE relative energy error | Experimental indicator |',
            '|---|---|---:|---:|']
    for row in rows:
        report.append(f'| {row["formulation"]} | {row["length_scale"]} | '
                      f'{row["true_relative_energy_error_percent"]:.7g}% | '
                      f'{row["experimental_eta_percent"]:.7g}% |')
    report += ['', 'The legacy six-DOF Hermite interpolation is **not affine-complete**.',
               'The 3-DOF Cauchy control and 6-DOF **interior bubble** enrichment are',
               'separate research formulations, NOT corrected physical Hermite rotations.',
               'True error is available here only because this manufactured field is exactly known.',
               'The experimental indicator is not certified; percentages are not interchangeable.']
    (args.out/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    print(f'Wrote {args.out / "REPORT.md"}',flush=True)

if __name__=='__main__':main()
