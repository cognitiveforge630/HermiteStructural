"""Small reproducible experimental a posteriori cantilever diagnostics.

Example: python -u run_aposteriori_cantilever.py --elements 1 2 3
No modulus calibration is performed: the supplied reference E stays fixed.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

from aposteriori import compute_residual_indicator
from aposteriori.vtk_export import write_indicator_vtk
from solver import (P, REFERENCE_LONGITUDINAL_E, build_beam,
                    build_end_face_traction_load_vector, free_end_uy,
                    solve_with_force_vector)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--elements', type=int, nargs=3, metavar=('EX', 'EY', 'EZ'), default=(1, 2, 3))
    parser.add_argument('--formulation', choices=('legacy_hermite', 'cauchy_hex8_reference', 'affine_enriched_hex8',
                                                  'recovered_rotation_hermite','quadratic_rotation_hermite'),
                        default='legacy_hermite')
    parser.add_argument('--length-scale', choices=('diameter','geometric_mean','min_edge'),
                        default='diameter')
    parser.add_argument('--out', type=Path, default=None)
    parser.add_argument('--body-order', type=int, default=None)
    parser.add_argument('--nodal-clamp-only',action='store_true',
                        help='Research diagnostic: skip strong face trace; NOT recommended for nonlocal Hermite')
    args = parser.parse_args()
    if any(v < 1 for v in args.elements):
        parser.error('Each element count must be positive')
    if args.out is None:
        args.out = Path('outputs/aposteriori_demo' if args.formulation == 'legacy_hermite'
                        else ('outputs/aposteriori_cauchy_reference' if args.formulation == 'cauchy_hex8_reference'
                              else ('outputs/aposteriori_affine_enriched' if args.formulation == 'affine_enriched_hex8'
                                    else f'outputs/aposteriori_{args.formulation}')))
    args.out.mkdir(parents=True, exist_ok=True)
    print(f'Building {args.formulation} cantilever with fixed reference Ez...', flush=True)
    ex, ey, ez = args.elements
    beam = build_beam(REFERENCE_LONGITUDINAL_E, nx=ex+1, ny=ey+1, nz=ez+1,
                      formulation=args.formulation)
    force, _, traction = build_end_face_traction_load_vector(beam)
    print(f'Solving {beam.get_num_elements()} elements, {len(force)} DOFs...', flush=True)
    strong_clamp=hasattr(beam,'element_transform') and not args.nodal_clamp_only
    if strong_clamp:
        from face_dirichlet import solve_with_face_dirichlet
        U, clamp_info=solve_with_face_dirichlet(beam, force, lambda x: (0.,0.,0.),
                                                lambda x: (0.,0.,0.))
        print(f'Strong face trace enforced with {clamp_info["independent_trace_constraints"]} independent constraints',flush=True)
    else:
        U=solve_with_force_vector(beam,force)
        clamp_info=None
    print('Integrating volume, face and energy diagnostics...', flush=True)
    result = compute_residual_indicator(beam, U, end_traction=(0, traction, 0),
                                        body_order=args.body_order, length_scale=args.length_scale)
    result.update({'formulation': args.formulation, 'mesh_elements_xyz': args.elements, 'end_load_lbf': P,
                   'end_traction_psi': traction, 'material_Ez_psi': REFERENCE_LONGITUDINAL_E,
                   'tip_uy_in': float(free_end_uy(beam, U)), 'calibrated': False,
                   'strong_face_dirichlet':strong_clamp,'face_constraint_stats':clamp_info})
    write_indicator_vtk(beam, U, result, args.out / 'experimental_indicators.vtk')
    (args.out / 'result.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    with (args.out / 'element_indicators.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['element_id', 'eta_squared_lbf_in', 'eta_sqrt_lbf_in', 'volume_in3',
                         'volume_residual_squared_lbf_in','jump_squared_lbf_in','natural_squared_lbf_in'])
        for e, (sq, eta, vol, vr, fj, bn) in enumerate(zip(
                result['element_indicator_squared_lbf_in'],
                result['element_indicator_sqrt_lbf_in'],
                result['element_volumes_in3'],
                result['element_volume_residual_squared_lbf_in'],
                result['element_jump_squared_lbf_in'],
                result['element_natural_squared_lbf_in'])):
            writer.writerow([e, sq, eta, vol, vr, fj, bn])
    label = 'EXPERIMENTAL — NOT VALIDATED OR CERTIFIED'
    pct = result['estimated_energy_error_percent_experimental']
    summary = (f'# Cantilever a posteriori prototype\n\n**{label}**\n\n'
               f'- Formulation: {args.formulation}\n'
               f'- Residual scaling: {args.length_scale} (not certified)\n'
               f'- Elements: {ex} x {ey} x {ez}\n'
               f'- Tip uy: {result["tip_uy_in"]:.9g} in\n'
               f'- Fixed longitudinal modulus: {REFERENCE_LONGITUDINAL_E:.9g} psi\n'
               f'- Experimental energy-normalized residual indicator: {pct:.6g}%\n'
               f'- Solution energy norm: {result["solution_energy_norm_sqrt_lbf_in"]:.9g} sqrt(lbf in)\n'
               f'- Volume residual squared: {result["volume_residual_squared_lbf_in"]:.9g} lbf in\n'
               f'- Interior traction jump squared: {result["interior_traction_jump_squared_lbf_in"]:.9g} lbf in\n'
               f'- Natural boundary mismatch squared: {result["natural_traction_mismatch_squared_lbf_in"]:.9g} lbf in\n\n'
               'This is **not** a certified upper bound or a percentage error in tip displacement or stress.\n'
               + ('The legacy Hermite bending/shear patch test fails and must be addressed.\n'
                if args.formulation == 'legacy_hermite' else
                'The estimator remains unvalidated; this is not a certified physical error percentage.\n'))
    (args.out / 'REPORT.md').write_text(summary, encoding='utf-8')
    print(summary, flush=True)
    print(f'Wrote diagnostics to {args.out.resolve()}', flush=True)


if __name__ == '__main__':
    main()
