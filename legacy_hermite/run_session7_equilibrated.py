"""Small complementary-energy estimates from saved solutions; no new FE solve.

Cantilever: analytically equilibrated rational stress plus traction-free modes.
Shear: constant symmetric stress reconstructed from prescribed natural loads.
Quadrature and floating-point defects are not enclosed, so no certified claim.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import time

from compare_session6_fields import read_solution
from run_session6_accuracy import exact_fields
from run_session6_accuracy import operators
from run_session7_matrix import interface_diagnostics
from solver import P, build_end_face_traction_load_vector
from aposteriori.equilibrated_cauchy import (CantileverAiryStress, stress_samples,
    constitutive_gap, constant_traction_stress, audit_cantilever_stress)


def evaluate_archive(path, order=5, transverse_degree=3, span_degree=3):
    total_started=time.perf_counter()
    beam, U, meta = read_solution(path)
    phases={'archive_and_model':time.perf_counter()-total_started}
    phase_started=time.perf_counter()
    trace=0.
    for e in beam.get_elements_at_min_z():
        for a,b in ((-1.,1.),(-.37,.61),(.43,-.29),(1.,-1.)):
            ids,N,_,_=operators(beam,e,(a,b,-1.))
            trace=max(trace,float(np.max(np.abs(N @ U[ids]))))
    kinematics=interface_diagnostics(beam,U)
    kinematics['maximum_sampled_clamp_displacement_in']=trace
    if trace>1e-8 or kinematics['maximum_sampled_displacement_jump_in']>1e-8:
        raise ValueError('Displacement field fails sampled kinematic admissibility audit')
    phases['kinematic_audit']=time.perf_counter()-phase_started
    phase_started=time.perf_counter()
    points, weights, computed = stress_samples(beam, U, order)
    phases['fe_stress_sampling']=time.perf_counter()-phase_started
    phase_started=time.perf_counter()
    compliance = np.linalg.inv(beam.D)
    norm_sq = float(np.einsum('nc,nc,n->',computed @ compliance,computed,weights))
    if norm_sq <= 0:
        raise ValueError('Computed energy norm is zero')
    if meta['problem'] == 'cantilever':
        # Only the specified fixed-load problem is supported, not arbitrary archives.
        expected_load, _, _ = build_end_face_traction_load_vector(beam)
        with np.load(path, allow_pickle=False) as archive:
            np.testing.assert_allclose(archive['load'], expected_load, rtol=1e-10, atol=1e-10)
        field = CantileverAiryStress(beam, P, transverse_degree, span_degree)
        base_sq = constitutive_gap(weights, field.base(points), computed, compliance)
        phases['load_check_and_particular_stress']=time.perf_counter()-phase_started
        phase_started=time.perf_counter()
        fit = field.fit(points, weights, computed, compliance)
        phases['stress_fit']=time.perf_counter()-phase_started
        phase_started=time.perf_counter()
        eta_sq = constitutive_gap(weights, field(points), computed, compliance)
        audit = audit_cantilever_stress(beam, field, P)
        if eta_sq > base_sq*(1+1e-10):
            raise ArithmeticError('Constrained stress optimization increased the gap')
        audit['gap_nonincrease'] = True
        if audit['maximum_natural_traction_mismatch_psi'] > 1e-7 or audit['maximum_sampled_divergence_lbf_in3'] > 1e-6:
            raise ArithmeticError(f'Stress admissibility audit failed: {audit}')
    elif meta['problem'] == 'shear':
        _, _, _, traction = exact_fields(beam, 'shear')
        recovered, mismatch = constant_traction_stress(beam, traction)
        eta_sq = constitutive_gap(weights, np.broadcast_to(recovered, computed.shape), computed, compliance)
        base_sq, fit = eta_sq, {'recovered_constant_stress_psi': recovered.tolist()}
        audit = {'maximum_natural_traction_mismatch_psi': mismatch, 'analytic_divergence_zero': True}
    else:
        raise ValueError('Supported loads: cantilever or zero-body-force shear only')
    eta, norm = np.sqrt(max(0., eta_sq)), np.sqrt(norm_sq)
    phases['final_gap_and_equilibrium_audit']=time.perf_counter()-phase_started
    phases['total']=time.perf_counter()-total_started
    exact_norm = meta.get('exact_energy_norm')
    true_error = meta.get('true_error_norm')
    interpretable = true_error is not None and exact_norm is not None and true_error > 1e-6*exact_norm
    return {'status': 'EQUILIBRATED_STRESS_RESEARCH_ESTIMATE', 'certified_bound': False,
            'solution': str(Path(path).resolve()), 'problem': meta['problem'],
            'clamp': meta['clamp'], 'formulation': meta['formulation'], 'mesh': meta['mesh'],
            'quadrature_order': order, 'eta_energy_norm': float(eta),
            'phase_timings_seconds':phases,
            'computed_energy_norm': float(norm),
            'estimated_energy_error_percent_computed_normalization': float(100*eta/norm),
            'unoptimized_gap_norm': float(np.sqrt(base_sq)),
            'effectivity_against_known_true_error': float(eta/true_error) if interpretable else None,
            'normalization': '100 * eta / computed_energy_norm',
            'estimated_conservative_exact_normalization_percent':
                float(100*eta/(norm-eta)) if norm>eta else None,
            'kinematic_audit': kinematics,
            'audit': audit, 'reconstruction': fit,
            'limitations': ['Numerical quadrature and floating point errors are not enclosed.',
                'Cantilever recovery subspace is independent of x and may give conservative estimates.',
                'No local error bound follows by discarding contributions near the clamp.',
                'Reference convergence and broad estimator validation remain required.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('solution', type=Path)
    parser.add_argument('--order', type=int, default=5)
    parser.add_argument('--transverse-degree', type=int, default=3)
    parser.add_argument('--span-degree', type=int, default=3)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.order < 4 or min(args.transverse_degree, args.span_degree) < 0:
        parser.error('Use order >=4 and nonnegative degrees')
    result = evaluate_archive(args.solution, args.order, args.transverse_degree, args.span_degree)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
