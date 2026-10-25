"""Boundary audit and exact-solution accuracy studies; no certified estimator.

Defaults are deliberately small. Larger mesh sweeps are user-run. Each case
is saved immediately; solution archives permit later independent comparison.
The sine/shear problems have side tractions and are NOT the cantilever load.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import time

import numpy as np

from aposteriori.residual import FACES, face_geometry, face_parent_point, voigt_to_stress
from face_dirichlet import solve_with_face_dirichlet
from regional_energy import RegionalEnergy
from accuracy_models import build_accuracy_beam
from run_session5_refinement import parse_mesh, mesh_name, nodal_rigid_virtual_work
from solver import (build_beam, build_end_face_traction_load_vector,
                    solve_with_force_vector, solve_with_prescribed_displacement, REFERENCE_LONGITUDINAL_E)


def operators(beam, e, p):
    """Full support displacement and engineering-strain operators."""
    xyz = beam.get_element_nodes(e)
    if hasattr(beam, 'element_transform'):
        ids, transform = beam.element_transform(e)
        B, det, N = beam.primitive_B(*p, xyz)
        return ids, N @ transform, B @ transform, det
    ids = beam.get_element_global_dofs(e).ravel()
    B, det, N, _ = beam.build_B_matrix(*p, xyz)
    return ids, N, B, det


def exact_fields(beam, problem):
    """u_x=f(z), epsilon_xz=f'(z); orthotropic axes as in solver.py."""
    if problem not in ('shear', 'sine'):
        raise ValueError('Manufactured fields require shear or sine')
    if not np.allclose(beam.D[:5, 5], 0.):
        raise ValueError('Manufactured body load assumes uncoupled xz shear')
    amplitude = 1e-3
    wave = np.pi / (2 * beam.Lz)

    def values(x):
        z = x[2]
        if problem == 'shear':
            return amplitude*z/beam.Lz, amplitude/beam.Lz, 0.
        return (amplitude*np.sin(wave*z), amplitude*wave*np.cos(wave*z),
                -amplitude*wave**2*np.sin(wave*z))

    def displacement(x):
        return np.array([values(x)[0], 0., 0.])

    def strain(x):
        return np.array([0., 0., 0., 0., 0., values(x)[1]])

    def body(x):
        # For this orthotropic D, only xz shear stress is generated.
        return np.array([-beam.D[5, 5]*values(x)[2], 0., 0.])

    def traction(x, n):
        return voigt_to_stress(beam.D @ strain(x)) @ n

    return displacement, strain, body, traction


def manufactured_load(beam, body, traction, order=6):
    """Independent integration of body load and ALL natural faces.

    Integrate nonlocal loads in primitive coordinates, then apply T.T once per
    element. Cache primitive matrices across equal-size axis-aligned bricks.
    This preserves quadrature and loads while avoiding dense B/N @ T at every
    point. As elsewhere in this driver, geometry is a structured brick mesh.
    """
    load = np.zeros(len(beam.grid_points)*beam.ndof_per_node)
    faces = getattr(beam, 'faces', FACES)
    counts = Counter(tuple(sorted(int(nodes[n]) for n in local))
                     for nodes in beam.elements for _, _, local in faces)
    gp, gw = np.polynomial.legendre.leggauss(order)
    matrix_cache = {}
    nonlocal_element = hasattr(beam, 'element_transform')
    for e, nodes in enumerate(beam.elements):
        xyz = beam.get_element_nodes(e)
        geometry_key = tuple(np.round(np.ptp(xyz, axis=0), 12))
        if nonlocal_element:
            ids, transform = beam.element_transform(e)
            element_load = np.zeros(transform.shape[0])
        else:
            ids = beam.get_element_global_dofs(e).ravel()
            transform = None
            element_load = np.zeros(len(ids))

        def displacement_matrix(p):
            key = (geometry_key, tuple(float(v) for v in p))
            if key not in matrix_cache:
                matrix_cache[key] = (beam.primitive_matrices(*p, xyz)[0] if nonlocal_element
                                     else beam.get_hermite_displacement_matrices(*p, xyz)[0])
            return matrix_cache[key]

        for i, xi in enumerate(gp):
            for j, eta in enumerate(gp):
                for k, zeta in enumerate(gp):
                    p = (xi, eta, zeta)
                    N = displacement_matrix(p)
                    det = float(np.linalg.det(beam.get_hermite_jacobian(*p, xyz)))
                    if det <= 0 or not np.isfinite(det):
                        raise ValueError('Invalid element Jacobian during load integration')
                    x = beam.hex8_shape_functions(*p) @ xyz
                    element_load += N.T @ body(x)*det*gw[i]*gw[j]*gw[k]
        for axis, side, local in faces:
            if counts[tuple(sorted(int(nodes[n]) for n in local))] != 1:
                continue
            if axis == 2 and side == -1 and np.allclose(xyz[list(local), 2], 0.):
                continue
            for i, s in enumerate(gp):
                for j, t in enumerate(gp):
                    p = face_parent_point(axis, side, s, t)
                    n, area = face_geometry(beam, xyz, p, axis, side)
                    N = displacement_matrix(p)
                    x = beam.hex8_shape_functions(*p) @ xyz
                    element_load += N.T @ traction(x, n)*area*gw[i]*gw[j]
        load[ids] += transform.T @ element_load if transform is not None else element_load
    return load


def exact_error(beam, U, exact_strain, order=6, *, exclude_clamp_lengths=()):
    gp, gw = np.polynomial.legendre.leggauss(order)
    regions = RegionalEnergy(beam, exclude_clamp_lengths)
    matrix_cache = {}
    nonlocal_element = hasattr(beam, 'element_transform')
    for e in range(beam.get_num_elements()):
        xyz = beam.get_element_nodes(e)
        geometry_key = tuple(np.round(np.ptp(xyz, axis=0), 12))
        if nonlocal_element:
            ids, transform = beam.element_transform(e)
            coefficients = transform @ U[ids]
        else:
            coefficients = U[beam.get_element_global_dofs(e).ravel()]
        lo, hi = float(xyz[:, 2].min()), float(xyz[:, 2].max())
        breaks = regions.split(lo, hi)
        for zlo, zhi in zip(breaks[:-1], breaks[1:]):
            # Split a cut element exactly rather than masking its Gauss points.
            center = (zlo+zhi)/2
            parent_center = 2*(center-lo)/(hi-lo)-1
            parent_scale = (zhi-zlo)/(hi-lo)
            for i, xi in enumerate(gp):
                for j, eta in enumerate(gp):
                    for k, zeta in enumerate(gp):
                        p = (xi, eta, parent_center+parent_scale*zeta)
                        key = (geometry_key, tuple(float(v) for v in p))
                        if key not in matrix_cache:
                            if nonlocal_element:
                                B, det, _ = beam.primitive_B(*p, xyz)
                            else:
                                B, det, _, _ = beam.build_B_matrix(*p, xyz)
                            matrix_cache[key] = B, det
                        B, det = matrix_cache[key]
                        exact = exact_strain(beam.hex8_shape_functions(*p) @ xyz)
                        diff = B @ coefficients - exact
                        weight = det*parent_scale*gw[i]*gw[j]*gw[k]
                        regions.add(float(diff @ beam.D @ diff)*weight,
                                    float(exact @ beam.D @ exact)*weight, center)
    error_sq, norm_sq = regions.total
    return {'true_energy_error_percent': 100*np.sqrt(error_sq/norm_sq) if norm_sq > 0 else None,
            'true_error_norm': np.sqrt(error_sq), 'exact_energy_norm': np.sqrt(norm_sq),
            'regional_true_errors': regions.report()}


def run_case(problem, mesh, formulation, clamp, order=6, *, z_breaks=None,linear_solver='scipy'):
    start = time.perf_counter()
    beam = build_accuracy_beam(mesh, formulation, z_breaks=z_breaks)
    phases={'model_setup':time.perf_counter()-start}
    phase_start=time.perf_counter()
    fields = None
    if problem == 'cantilever':
        load, _, _ = build_end_face_traction_load_vector(beam)
    else:
        fields = exact_fields(beam, problem)
        load = manufactured_load(beam, fields[2], fields[3], order)
    phases['load_integration']=time.perf_counter()-phase_start
    phase_start=time.perf_counter()
    K=beam.build_global_K()
    phases['stiffness_assembly']=time.perf_counter()-phase_start
    phase_start=time.perf_counter()
    if formulation in ('cauchy_hex8_reference', 'cauchy_hex27_reference'):
        if clamp != 'displacement_only':
            raise ValueError('Cauchy has no rotational DOFs')
        U = solve_with_prescribed_displacement(beam,load,stiffness=K,linear_solver=linear_solver)
        meta = {}
    else:
        U, meta = solve_with_face_dirichlet(beam, load, lambda x: np.zeros(3),
            constrain_rotations=(clamp == 'historical_zero_spin'),stiffness=K,linear_solver=linear_solver)
    phases['constraints_and_solve']=time.perf_counter()-phase_start
    phase_start=time.perf_counter()
    energy_sq = float(U @ (K @ U))
    work = float(load @ U)
    # Uniform unit end traction gives an actual area average, unlike nodal mean.
    unit_end_load, _, _ = build_end_face_traction_load_vector(beam, total_load=1.)
    end_average = float(unit_end_load @ U)
    unit_x_load, _, _ = build_end_face_traction_load_vector(beam, total_load=1., load_dof=0)
    reaction = K @ U - load
    force, moment = nodal_rigid_virtual_work(beam, load)
    rf, rm = nodal_rigid_virtual_work(beam, reaction)
    # Check trace at points different from the constraint quadrature, including edges.
    trace = 0.
    for e in beam.get_elements_at_min_z():
        for s in (-1., -.43, .29, 1.):
            for t in (-1., -.27, .61, 1.):
                ids, N, _, _ = operators(beam, e, (s, t, -1.))
                trace = max(trace, float(np.max(np.abs(N @ U[ids]))))
    row = {'problem': problem, 'mesh': list(mesh), 'formulation': formulation,
           'z_breaks_in': list(z_breaks) if z_breaks is not None else None,
           'clamp': clamp, 'dofs': len(U), 'quadrature_order': order,
           'assembly_quadrature_order': beam.gauss_order,
           'linear_solver_backend':linear_solver,'linear_solve':getattr(beam,'last_linear_solve',None),
           'assembly_statistics':getattr(beam,'last_assembly_statistics',None),
           'estimated_error_percent': None, 'certified_bound': False,
           'true_error_available': fields is not None,
           'end_face_average_uy_in': end_average,
           'end_face_average_ux_in': float(unit_x_load @ U),
           'energy_norm': np.sqrt(energy_sq),
           'work_energy_relative_mismatch': abs(energy_sq-work)/max(abs(work), 1e-30),
           'maximum_independent_clamp_trace_in': trace,
           'force_balance_lbf': (np.array(force)+rf).tolist(),
           'moment_balance_lbf_in': (np.array(moment)+rm).tolist(),
           'load_resultant_lbf': force, 'load_moment_lbf_in': moment,
           'fixed_modulus_Ez_psi': REFERENCE_LONGITUDINAL_E,
           'face_constraint_stats': meta}
    if hasattr(beam, 'quadratic_rank'):
        row['minimum_recovery_rank'] = min(beam.quadratic_rank)
    if hasattr(beam, 'cubic_rank'):
        row['cubic_recovery_terms'] = beam.cubic_terms
        row['minimum_cubic_recovery_rank'] = min(beam.cubic_rank)
        row['full_3d_cubic_completeness'] = len(beam.cubic_terms) == 19
        row['maximum_cubic_recovery_condition_number'] = max(beam.cubic_condition)
    phases['diagnostics']=time.perf_counter()-phase_start
    phase_start=time.perf_counter()
    if fields is not None:
        row.update(exact_error(beam, U, fields[1], order))
    phases['true_error_integration']=time.perf_counter()-phase_start
    row['phase_timings_seconds']=phases
    row['seconds'] = time.perf_counter()-start
    return row, beam, U, load


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--problem', choices=['shear', 'sine', 'cantilever'], default='cantilever')
    parser.add_argument('--meshes', nargs='+', type=parse_mesh, default=[(2, 2, 2)])
    parser.add_argument('--formulations', nargs='+', choices=['quadratic_rotation_hermite', 'cubic_rotation_hermite', 'cauchy_hex8_reference', 'cauchy_hex27_reference'],
                        default=['quadratic_rotation_hermite', 'cauchy_hex8_reference'])
    parser.add_argument('--clamps', nargs='+', choices=['displacement_only', 'historical_zero_spin'],
                        default=['displacement_only', 'historical_zero_spin'])
    parser.add_argument('--quadrature-order', type=int, default=6)
    parser.add_argument('--linear-solver',choices=['scipy','pypardiso'],default='scipy')
    parser.add_argument('--z-breaks', nargs='+', type=float, help='Optional EZ+1 coordinates, for one mesh only')
    parser.add_argument('--out', type=Path, default=Path('outputs/session6_accuracy'))
    args = parser.parse_args()
    if args.quadrature_order < 4:
        parser.error('Use quadrature order >=4')
    if args.z_breaks is not None and len(args.meshes) != 1:
        parser.error('Use one mesh when specifying z-breaks')
    if any(min(m) < 2 for m in args.meshes):
        parser.error('Use >=2 elements in each direction to retain rank-9 recovery')
    if 'cubic_rotation_hermite' in args.formulations and any(m[2]<3 for m in args.meshes):
        parser.error('Cubic bending recovery requires >=2x2x3 elements')
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for mesh in args.meshes:
        for form in args.formulations:
            clamps = args.clamps if form in ('quadratic_rotation_hermite','cubic_rotation_hermite') else ['displacement_only']
            for clamp in clamps:
                name = f'{args.problem}_{mesh_name(mesh)}_{form}_{clamp}'
                if args.z_breaks is not None:
                    name += '_graded'
                print(f'Starting {name}', flush=True)
                row, beam, U, load = run_case(args.problem, mesh, form, clamp, args.quadrature_order,
                                             z_breaks=args.z_breaks,linear_solver=args.linear_solver)
                rows.append(row)
                (args.out / f'{name}.json').write_text(json.dumps(row, indent=2), encoding='utf-8')
                np.savez_compressed(args.out / f'{name}.npz', U=U, load=load,
                                    coordinates=beam.grid_points, elements=beam.elements,
                                    D=beam.D, metadata=json.dumps(row))
                (args.out / f'{args.problem}_summary.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
                print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
