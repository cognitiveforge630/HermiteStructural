"""Compare cubic source-stress *section-only* research results without re-solving.

Unlike the visualization comparator, these section metrics use direct Gaussian
quadrature on the original FE stress. It does not compute a certified bound,
true 3D stress error, or full-field convergence norm.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def _finite(v, name):
    if type(v) is bool:
        raise ValueError(f'{name}: expected finite numerical value')
    try:
        out=float(v)
    except (ValueError,TypeError) as exc:
        raise ValueError(f'{name}: expected finite numerical value') from exc
    if not math.isfinite(out):
        raise ValueError(f'{name}: expected finite numerical value')
    return out


def _section_path(source, record):
    field=record.get('cubic_source_section_file')
    if not field:
        raise ValueError(f'{source}: missing source section report')
    reference=Path(field)
    candidates=[source.with_name(source.stem+'_source_sections.json'),
                source.parent/reference.name, ROOT/reference]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f'{source}: source section report was not found')


def read_case(path):
    path=Path(path)
    data=json.loads(path.read_text(encoding='utf8'))
    if data.get('formulation')!='cubic_rotation_hermite':
        raise ValueError(f'{path}: expected cubic_rotation_hermite')
    mesh=data.get('node_mesh')
    if not isinstance(mesh,list) or len(mesh)!=3 or any(type(x)!=int or x<3 for x in mesh):
        raise ValueError(f'{path}: invalid mesh')
    path_s=_section_path(path,data)
    report=json.loads(path_s.read_text(encoding='utf8'))
    if report.get('certified_error_bound') is not False or report.get('status')!='RESEARCH_SECTIONAL_DIAGNOSTIC_NOT_A_POSTERIORI':
        raise ValueError(f'{path_s}: invalid certification marker')
    if report.get('visualization_subdivisions_affect_results') is not False:
        raise ValueError(f'{path_s}: unexpected sampling dependency')
    quadrature_order=report.get('cross_section_quadrature_order')
    if type(quadrature_order)!=int or quadrature_order<2:
        raise ValueError(f'{path_s}: invalid quadrature order')
    sections=report.get('sections')
    positions=report.get('positions_parent_zeta')
    if not isinstance(positions,list) or not positions or not isinstance(sections,list) or len(sections)!=(mesh[2]-1)*len(positions):
        raise ValueError(f'{path_s}: missing section positions')
    if any(not -1<_finite(z,'parent zeta')<1 for z in positions):
        raise ValueError(f'{path_s}: invalid slice location')
    force=_finite(data['end_force_lbf'],'tip force')
    weight=_finite(data['awc_reference']['weight_per_length_lbf_in'],'weight per length')
    tip=_finite(data['fem_tip_uy_in'],'tip')
    ref_tip=_finite(data['awc_reference']['total_uy_in'],'reference tip')
    if ref_tip==0:
        raise ValueError(f'{path}: zero tip reference')
    return {
        'result_file':str(path),'source_section_file':str(path_s),
        'mesh':mesh,'section_order':quadrature_order,'positions':positions,
        'end_force_lbf':force,'weight_per_length_lbf_in':weight,
        'longitudinal_E_psi':_finite(data['longitudinal_E_psi'],'longitudinal E'),
        'boundary_condition':data['boundary_condition'],
        'beam_tip_uy_in':ref_tip, 'tip_uy_in':tip,
        'tip_vs_beam_pct':100*abs(tip-ref_tip)/abs(ref_tip),
        'maximum_midspan_shear_vs_beam_pct':_finite(report['max_midspan_shear_y_difference_pct'],'shear reference discrepancy'),
        'maximum_midspan_moment_vs_beam_pct':_finite(report['max_midspan_moment_x_difference_pct'],'moment reference discrepancy'),
        'sampled_stress_screen_performed':data.get('sampled_cubic_stress_screen_performed'),
        'algebraic_relative_backward_error':_finite(data['linear_solver_diagnostics']['linear_solve']['relative_backward_error'],'solver backward error'),
    }


def compare(cases):
    if not cases:
        raise ValueError('Supply at least one case')
    baseline=cases[0]
    for current in cases[1:]:
        for key in ('end_force_lbf','weight_per_length_lbf_in',
                    'longitudinal_E_psi','beam_tip_uy_in'):
            if not math.isclose(current[key],baseline[key],rel_tol=1e-9,abs_tol=1e-9):
                raise ValueError(f'Incompatible load/material/reference: {key}')
        for key in ('boundary_condition','section_order','positions'):
            if current[key]!=baseline[key]:
                raise ValueError(f'Incompatible section method or boundary conditions: {key}')
    pairs=[]
    for older,newer in zip(cases,cases[1:]):
        if older['mesh']==newer['mesh'] or any(y<x for x,y in zip(older['mesh'],newer['mesh'])):
            raise ValueError('Meshes must increase without coarsening any axis')
        ratios=[(x-1)/(y-1) for x,y in zip(older['mesh'],newer['mesh'])]
        pairs.append({
            'from_mesh':older['mesh'],'to_mesh':newer['mesh'],
            'ratios_h_fine_over_coarse_xyz':ratios,
            'uniform_refinement':(max(ratios)-min(ratios)<1e-12 and ratios[0]<1),
            'tip_change_pct_of_fine':100*abs(newer['tip_uy_in']-older['tip_uy_in'])/max(abs(newer['tip_uy_in']),1e-30),
            'tip_beam_discrepancy_decreased':newer['tip_vs_beam_pct']<older['tip_vs_beam_pct'],
            'maximum_shear_beam_discrepancy_decreased':newer['maximum_midspan_shear_vs_beam_pct']<older['maximum_midspan_shear_vs_beam_pct'],
        })
    return {
        'status':'RESEARCH_SECTION_QOI_COMPARISON_NOT_CERTIFIED',
        'certified_error_bound':False,
        'caveat':('The maximum discrepancies are at the available interior section locations, '
                  'which vary by mesh; their changes are NOT a field-error norm, bound or '
                  'proof of convergence. Beam theory is not the exact 3D reference.'),
        'cases':cases,'pairs':pairs,
    }


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('results',nargs='+',type=Path)
    p.add_argument('--json-out',type=Path,default=None)
    a=p.parse_args(argv)
    result=compare([read_case(p) for p in a.results])
    print('SOURCE FE SECTION COMPARISON — NOT A POSTERIORI CERTIFIED')
    print('Mesh        tip vs beam %   worst mid shear vs beam %   worst mid moment vs beam %')
    for c in result['cases']:
        print(f"{'x'.join(map(str,c['mesh'])):11s} {c['tip_vs_beam_pct']:13.4f} "
              f"{c['maximum_midspan_shear_vs_beam_pct']:26.4f} "
              f"{c['maximum_midspan_moment_vs_beam_pct']:28.4f}")
    for pair in result['pairs']:
        print(f"{pair['from_mesh']} -> {pair['to_mesh']}: tip QoI change "
              f"{pair['tip_change_pct_of_fine']:.4f}%, "
              f"uniform h: {pair['uniform_refinement']}, "
              f"tip vs beam improved: {pair['tip_beam_discrepancy_decreased']}")
    if a.json_out:
        a.json_out.parent.mkdir(parents=True,exist_ok=True)
        a.json_out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('Wrote',a.json_out)
    return result

if __name__=='__main__':main()
