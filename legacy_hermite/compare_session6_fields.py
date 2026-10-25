"""Compare saved Session 6 strain fields on a common integration partition.

This measures a reference-field difference, NOT a certified error estimate.
Neither archive is assumed converged. Only current structured brick models
are supported. Refinement grids need not be nested; integrate their overlay.
"""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from run_session6_accuracy import operators
from solver import build_beam
from regional_energy import RegionalEnergy
from accuracy_models import build_accuracy_beam
from field_sampling import StructuredField, tensor_quadrature


def read_solution(path):
    with np.load(path, allow_pickle=False) as archive:
        meta = json.loads(str(archive['metadata']))
        mesh = meta['mesh']
        beam = build_accuracy_beam(mesh, meta['formulation'], meta['fixed_modulus_Ez_psi'],
                                   z_breaks=meta.get('z_breaks_in'))
        for key, expected in [('coordinates', beam.grid_points), ('elements', beam.elements), ('D', beam.D)]:
            np.testing.assert_allclose(archive[key], expected, rtol=1e-12, atol=1e-12,
                                       err_msg='Archive cannot be reproduced by current model')
        U = archive['U'].copy()
        if U.shape != (len(beam.grid_points)*beam.ndof_per_node,) or not np.all(np.isfinite(U)):
            raise ValueError('Invalid archived solution')
    return beam, U, meta


def compare_fields(left, right, order=5, *, exclude_clamp_lengths=(), allow_boundary_mismatch=False):
    a, ua, ma = left
    b, ub, mb = right
    for key in ('problem', 'clamp', 'fixed_modulus_Ez_psi'):
        if ma[key] != mb[key]:
            if key == 'clamp' and allow_boundary_mismatch:
                continue
            raise ValueError(f'Cannot compare different physical problems: {key}')
    np.testing.assert_allclose(a.D, b.D, rtol=1e-12, atol=1e-12)
    low_a, high_a = a.grid_points.min(axis=0), a.grid_points.max(axis=0)
    np.testing.assert_allclose(low_a, b.grid_points.min(axis=0))
    np.testing.assert_allclose(high_a, b.grid_points.max(axis=0))
    sample_a,sample_b=StructuredField(a,ua),StructuredField(b,ub)
    cuts = [np.unique(np.r_[sample_a.planes[k],sample_b.planes[k]]) for k in range(3)]
    regions = RegionalEnergy(a, exclude_clamp_lengths)
    cuts[2] = np.unique(np.r_[cuts[2], regions.cutoffs])

    parent_points,parent_weights=tensor_quadrature(order)
    for index in itertools.product(*(range(len(c)-1) for c in cuts)):
        lo = np.array([cuts[k][index[k]] for k in range(3)])
        hi = np.array([cuts[k][index[k]+1] for k in range(3)])
        center, scale = (lo+hi)/2, (hi-lo)/2
        ea, la, ha = sample_a.locate(center)
        eb, lb, hb = sample_b.locate(center)
        x=center+scale*parent_points
        reference=sample_b.strain(eb,2*(x-lb)/(hb-lb)-1)
        diff=sample_a.strain(ea,2*(x-la)/(ha-la)-1)-reference
        weights=np.prod(scale)*parent_weights
        regions.add(float(np.einsum('qc,qc,q->',diff @ a.D,diff,weights)),
                    float(np.einsum('qc,qc,q->',reference @ a.D,reference,weights)),center[2])
    difference_sq, reference_sq = regions.total
    return {'status': ('BOUNDARY_SENSITIVITY_DIFFERENCE_REFERENCE_UNQUALIFIED' if ma['clamp'] != mb['clamp']
                       else 'REFERENCE_DIFFERENCE_ONLY_REFERENCE_UNQUALIFIED'),
            'solution_clamp': ma['clamp'], 'reference_clamp': mb['clamp'],
            'certified_bound': False, 'reference_uncertainty_known': False,
            'difference_energy_norm': float(np.sqrt(difference_sq)),
            'reference_energy_norm': float(np.sqrt(reference_sq)),
            'reference_relative_difference_percent':
                float(100*np.sqrt(difference_sq/reference_sq)) if reference_sq > 0 else None,
            'quadrature_order': order,
            'regional_comparisons': regions.report(),
            'sampling_cache': {'solution':sample_a.statistics(),'reference':sample_b.statistics()},
            'common_cells': int(np.prod([len(c)-1 for c in cuts]))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('solution', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('--order', type=int, default=5)
    parser.add_argument('--exclude-clamp-lengths', nargs='+', type=float, default=[],
                        help='Distances in inches from z=0; report global AND regional differences')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.order < 4:
        parser.error('Use order >=4; default 5 integrates these polynomial strain products')
    result = compare_fields(read_solution(args.solution), read_solution(args.reference), args.order,
                            exclude_clamp_lengths=args.exclude_clamp_lengths)
    result.update(solution=str(args.solution.resolve()), reference=str(args.reference.resolve()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
