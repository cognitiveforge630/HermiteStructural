"""Repeat small a posteriori experiments without changing/calibrating E.

Not a mesh convergence certificate. Each entry uses the same fixed material,
load and boundary model. Prints experimental indicators side by side.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

from aposteriori import compute_residual_indicator
from solver import (build_beam, build_end_face_traction_load_vector,
                    solve_with_force_vector, free_end_uy, REFERENCE_LONGITUDINAL_E)

DEFAULT_MESHES = [(1, 1, 2), (1, 2, 4), (2, 2, 6)]


def run_refinement(meshes=DEFAULT_MESHES, out=Path('outputs/aposteriori_refinement'),
                   formulation='legacy_hermite', length_scale='diameter'):
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for ex, ey, ez in meshes:
        print(f'[{ex}x{ey}x{ez}] Assembling / solving...', flush=True)
        beam = build_beam(REFERENCE_LONGITUDINAL_E, ex+1, ey+1, ez+1,
                          formulation=formulation)
        F, _, traction = build_end_face_traction_load_vector(beam)
        U = solve_with_force_vector(beam, F)
        print(f'[{ex}x{ey}x{ez}] Evaluating residual and face integrals...', flush=True)
        result = compute_residual_indicator(beam, U, end_traction=(0, traction, 0),
                                            length_scale=length_scale)
        row = {
            'formulation':formulation, 'length_scale':length_scale,
            'ex':ex, 'ey':ey, 'ez':ez, 'elements':beam.get_num_elements(),
            'dofs':U.size, 'tip_uy_in':float(free_end_uy(beam,U)),
            'eta_experimental_percent':result['estimated_energy_error_percent_experimental'],
            'eta_sqrt_lbf_in':result['eta_energy_units_sqrt_lbf_in'],
            'solution_energy_sqrt_lbf_in':result['solution_energy_norm_sqrt_lbf_in'],
            'volume_sq_lbf_in':result['volume_residual_squared_lbf_in'],
            'jump_sq_lbf_in':result['interior_traction_jump_squared_lbf_in'],
            'natural_sq_lbf_in':result['natural_traction_mismatch_squared_lbf_in'],
        }
        rows.append(row)
        print(f'[{ex}x{ey}x{ez}] Experimental eta: {row["eta_experimental_percent"]:.4f}% '
              f'tip uy: {row["tip_uy_in"]:.9g} in', flush=True)
    with (out/'refinement.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    (out/'refinement.json').write_text(json.dumps({
        'status':'EXPERIMENTAL_UNVALIDATED',
        'formulation': formulation,
        'length_scale': length_scale,
        'fixed_material_Ez_psi': REFERENCE_LONGITUDINAL_E,
        'rows':rows},indent=2)+'\n',encoding='utf-8')
    lines=['# Experimental residual refinement screening', '',
           f'Formulation: {formulation}; length-scale convention: {length_scale}.', '',
           '**Not a validated error estimate, accuracy claim, or convergence certificate.**','',
           '| Mesh elements | Tip uy (in) | Experimental eta (%) |',
           '|---|---:|---:|']
    for r in rows:
        lines.append(f'| {r["ex"]}x{r["ey"]}x{r["ez"]} | {r["tip_uy_in"]:.9g} | {r["eta_experimental_percent"]:.6g} |')
    lines += ['', 'The legacy Hermite element fails its affine constant-strain patch tests;',
              'the Cauchy and bubble-enriched references pass, but are distinct research formulations.',
              'No reference energy error was computed; effectivity is unknown.',
              'See docs/ESTIMATOR_STATUS.md for assumptions and interpretation.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--formulation', choices=('legacy_hermite','cauchy_hex8_reference','affine_enriched_hex8'),
                        default='legacy_hermite')
    parser.add_argument('--length-scale',choices=('diameter','geometric_mean','min_edge'),
                        default='diameter')
    parser.add_argument('--out',type=Path,default=None)
    parser.add_argument('--meshes',nargs='*',default=['1x1x2','1x2x4','2x2x6'])
    args=parser.parse_args()
    try:
        meshes=[tuple(map(int,s.lower().split('x'))) for s in args.meshes]
        if any(len(m)!=3 or min(m)<1 for m in meshes): raise ValueError()
    except (ValueError, TypeError):
        parser.error('Meshes must be positive EXxEYxEZ triplets, e.g. 1x2x4')
    if args.out is None:
        args.out = Path('outputs/aposteriori_refinement' if args.formulation == 'legacy_hermite'
                        else ('outputs/aposteriori_refinement_cauchy_reference'
                              if args.formulation == 'cauchy_hex8_reference'
                              else 'outputs/aposteriori_refinement_affine_enriched'))
    run_refinement(meshes,args.out,args.formulation,args.length_scale)


if __name__=='__main__':
    main()
