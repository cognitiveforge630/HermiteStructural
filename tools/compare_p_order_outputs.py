"""Print and export aligned research comparisons from Hermite, Hex27, Hex64 JSON.

Usage: python -m tools.compare_p_order_outputs --cubic outputs/cubic_dual_shear.json \\
   --hex27 outputs/hex27_3x5x25.json --hex64 outputs/hex64_3x5x25.json
"""
import argparse
import csv
import json
from pathlib import Path
import math


def load_row(label, path, is_cubic=False):
    d=json.loads(Path(path).read_text(encoding='utf-8'))
    if is_cubic:
        if d.get('formulation')!='cubic_rotation_hermite':
            raise ValueError(f'{path} is not cubic Hermite')
        r={
            'scheme':label,'node_mesh':d['node_mesh'],
            'E_z_psi':d['longitudinal_E_psi'],'end_force_lbf':d['end_force_lbf'],
            'total_weight_lbf':d['awc_reference']['total_weight_lbf'],
            'tip_uy_in':d['fem_tip_uy_in'],
            'beam_tip_discrepancy_pct':d['difference_vs_awc_pct'],
            'midspan_section_shear_discrepancy_pct':d['cubic_source_max_midspan_shear_y_difference_pct'],
            'midspan_section_moment_discrepancy_pct':d['cubic_source_max_midspan_moment_x_difference_pct'],
            'boundary_condition':d['boundary_condition'],
            'DOFs':'6 per Hermite vertex',
        }
    else:
        r={
            'scheme':label,'node_mesh':d['original_Hermite_vertex_mesh'],
            'E_z_psi':d['longitudinal_E_psi'],'end_force_lbf':d['end_force_lbf'],
            'total_weight_lbf':d['density_lbf_in3']*4*8*120,
            'tip_uy_in':d['tip_max_abs_uy_in'],
            'beam_tip_discrepancy_pct':d['tip_max_vs_classical_beam_pct'],
            'midspan_section_shear_discrepancy_pct':d['max_midspan_shear_y_difference_pct'],
            'midspan_section_moment_discrepancy_pct':d['max_midspan_moment_x_difference_pct'],
            'boundary_condition':d['clamp'],
            'DOFs':d['free_and_fixed_translation_DOFs'],
        }
    return r


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cubic',type=Path,required=True)
    p.add_argument('--hex27',type=Path,required=True)
    p.add_argument('--hex64',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path('outputs/p_order_side_by_side.csv'))
    a=p.parse_args(argv)
    rows=[load_row('cubic_Hermite_Hex8',a.cubic,True),
          load_row('conventional_Hex27',a.hex27),
          load_row('conventional_Hex64',a.hex64)]
    for key in ('node_mesh','E_z_psi','end_force_lbf','total_weight_lbf'):
        values=[r[key] for r in rows]
        if key=='node_mesh':
            matches=all(v==values[0] for v in values)
        else:
            matches=all(math.isclose(v,values[0],rel_tol=1e-10,abs_tol=1e-10) for v in values)
        if not matches:
            raise ValueError(f'Cannot compare incompatible {key}: {values}')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('w',newline='',encoding='utf-8') as file:
        w=csv.DictWriter(file,fieldnames=rows[0].keys())
        w.writeheader()
        for row in rows:
            w.writerow({**row,'node_mesh':'x'.join(map(str,row['node_mesh']))})
    print('Scheme | tip uy (in) | beam tip diff % | midspan shear diff % | moment diff %')
    for r in rows:
        print(f"{r['scheme']:24s} | {r['tip_uy_in']: .9f} | "
            f"{r['beam_tip_discrepancy_pct']:.4f} | "
            f"{r['midspan_section_shear_discrepancy_pct']:.4f} | "
            f"{r['midspan_section_moment_discrepancy_pct']:.4f}")
    print('WARNING: The existing cubic result may fix face rotations; the conventional '
          'elements fix only the face displacement. Re-run cubic with '
          '--displacement-only-clamp for a closer boundary-condition comparison.')
    print('None of these beam-theory discrepancies constitutes a certified 3D error bound.')
    print('CSV:',a.out)
    return rows

if __name__=='__main__':
    main()
