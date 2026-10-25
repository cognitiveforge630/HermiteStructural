"""Compare archived solutions against ONE reference; no FE solves.

Reports disjoint support/interior/free-end energy differences and mean-tip QoI
differences. A reference candidate is never automatically labelled exact.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np

from aposteriori.residual import FACES, face_parent_point, voigt_to_stress
from compare_session6_fields import read_solution, compare_fields
from run_session6_accuracy import operators
from solver import build_end_face_traction_load_vector
from field_sampling import StructuredField


def interface_diagnostics(beam, U):
    sampler=StructuredField(beam,U)
    faces = getattr(beam,'faces',FACES)
    owners = {}
    for e,nodes in enumerate(beam.elements):
        for axis,side,local in faces:
            owners.setdefault(tuple(sorted(int(nodes[i]) for i in local)),[]).append((e,axis,side))
    max_u, max_t = 0.,0.
    for values in owners.values():
        if len(values) != 2:
            continue
        (ea,axis,sa),(eb,_,sb) = values
        pairs=((- .37,.23),(.61,-.43))
        pa=np.array([face_parent_point(axis,sa,s,t) for s,t in pairs])
        pb=np.array([face_parent_point(axis,sb,s,t) for s,t in pairs])
        ua,ea_strain=sampler.displacement_and_strain(ea,pa)
        ub,eb_strain=sampler.displacement_and_strain(eb,pb)
        max_u=max(max_u,float(np.max(np.linalg.norm(ua-ub,axis=1))))
        sigma_jump=(ea_strain-eb_strain) @ beam.D.T
        traction_columns=((0,3,5),(3,1,4),(5,4,2))[axis]
        max_t=max(max_t,float(np.max(np.linalg.norm(sigma_jump[:,traction_columns],axis=1))))
    return {'maximum_sampled_displacement_jump_in':max_u,
            'maximum_sampled_interior_traction_jump_psi':max_t}


def band_report(comparison, length, support, free_end):
    regions = comparison['regional_comparisons']
    support_values, end_values = regions[0]['excluded'],regions[1]['retained']
    total_d = comparison['difference_energy_norm']**2
    total_r = comparison['reference_energy_norm']**2
    middle_d = max(0., total_d-support_values['difference_energy_norm_squared_lbf_in']-
                    end_values['difference_energy_norm_squared_lbf_in'])
    middle_r = max(0., total_r-support_values['reference_energy_norm_squared_lbf_in']-
                    end_values['reference_energy_norm_squared_lbf_in'])
    middle = {'difference_energy_norm_squared_lbf_in':middle_d,
              'reference_energy_norm_squared_lbf_in':middle_r,
              'relative_energy_difference_percent':100*np.sqrt(middle_d/middle_r) if middle_r>0 else None,
              'volume_percent':100*(length-support-free_end)/length,
              'share_of_global_squared_difference_percent':100*middle_d/total_d if total_d>0 else None}
    return [{'name':'support','z_range_in':[0.,support],**support_values},
            {'name':'interior','z_range_in':[support,length-free_end],**middle},
            {'name':'free_end','z_range_in':[length-free_end,length],**end_values}]


def build_row(left, reference, support=8., free_end=8., order=5):
    beam,U,meta = left
    ref_beam,ref_U,ref_meta = reference
    if support <= 0 or free_end <= 0 or support+free_end >= beam.Lz:
        raise ValueError('Choose positive, nonoverlapping support and free-end lengths')
    comparison = compare_fields(left,reference,order,
        exclude_clamp_lengths=[support,beam.Lz-free_end], allow_boundary_mismatch=True)
    load,_,_ = build_end_face_traction_load_vector(beam,total_load=1.)
    ref_load,_,_ = build_end_face_traction_load_vector(ref_beam,total_load=1.)
    tip,ref_tip = float(load @ U),float(ref_load @ ref_U)
    comparison.update(formulation=meta['formulation'],clamp=meta['clamp'],mesh=meta['mesh'],
        reference_formulation=ref_meta['formulation'],reference_mesh=ref_meta['mesh'],
        z_breaks_in=meta.get('z_breaks_in'), end_face_average_uy_in=tip,
        reference_end_face_average_uy_in=ref_tip,
        end_face_uy_reference_difference_percent=100*abs(tip-ref_tip)/abs(ref_tip) if ref_tip else None,
        disjoint_regions=band_report(comparison,beam.Lz,support,free_end),
        interface_diagnostics=interface_diagnostics(beam,U))
    return comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solutions', nargs='+',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--support-length',type=float,default=8.)
    parser.add_argument('--free-end-length',type=float,default=8.)
    parser.add_argument('--order',type=int,default=5)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.order<4:
        parser.error('Use integration order >=4')
    args.out.mkdir(parents=True,exist_ok=True)
    reference = read_solution(args.reference)
    rows = []
    for i,path in enumerate(args.solutions):
        print(f'Comparing {path.name}',flush=True)
        row = build_row(read_solution(path),reference,args.support_length,args.free_end_length,args.order)
        row['solution'],row['reference'] = str(path.resolve()),str(args.reference.resolve())
        rows.append(row)
        (args.out/f'comparison_{i+1}.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
        (args.out/'matrix.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
    lines = ['# Fixed-reference diagnostic matrix','',
        '**Reference accuracy is unqualified; percentages below are field/QoI differences.**','',
        '| Formulation | Clamp | Mesh | Global energy difference | Mean tip difference |',
        '|---|---|---|---:|---:|']
    for r in rows:
        lines.append(f"| {r['formulation']} | {r['clamp']} | {r['mesh']} | "
                     f"{r['reference_relative_difference_percent']:.5f}% | "
                     f"{r['end_face_uy_reference_difference_percent']:.5f}% |")
    lines += ['', 'The JSON files include disjoint support/interior/free-end contributions,',
              'their absolute squared energy norms, and sampled interface diagnostics.',
              'Boundary sensitivity comparisons explicitly record differing clamp modes.']
    (args.out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
