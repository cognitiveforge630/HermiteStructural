"""Session 5: independently check exact bending and cantilever mesh change.

Bending: exact 3D orthotropic Cauchy elasticity solution, so energy/stress /
corner displacement errors are TRUE errors. The research residual indicator is
not certified; effectivity is inapplicable for numerically exact solutions.

Cantilever: fixed-modulus end-face traction; no exact 3D solution is asserted.
Compare mesh-to-mesh changes, optional independent Cauchy comparator, load
virtual work, and experimental residual diagnostics. Never label changes as
true errors or inferred stress-error percentages.

Examples (from the project root, globally installed Python):
    python -u run_session5_refinement.py --quick
    python -u run_session5_refinement.py --bending-meshes 2x2x4 2x2x6 \\
        --cantilever-meshes 1x2x4 2x2x4 --vtk
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from aposteriori.residual import compute_residual_indicator, stress_at_element
from run_pure_bending_manufactured import evaluate
from solver import (P, REFERENCE_LONGITUDINAL_E, build_beam,
                    build_end_face_traction_load_vector, free_end_uy,
                    solve_with_force_vector)

FORMULATIONS = ('legacy_hermite', 'cauchy_hex8_reference', 'affine_enriched_hex8',
                'recovered_rotation_hermite', 'quadratic_rotation_hermite')
NONLOCAL = ('recovered_rotation_hermite','quadratic_rotation_hermite')


def parse_mesh(text):
    try:
        mesh=tuple(int(x) for x in text.lower().split('x'))
    except ValueError as e:
        raise argparse.ArgumentTypeError('Use EXxEYxEZ, e.g. 2x2x4') from e
    if len(mesh)!=3 or min(mesh)<1:
        raise argparse.ArgumentTypeError('Use three positive element counts EXxEYxEZ')
    return mesh


def mesh_name(mesh):
    return 'x'.join(map(str,mesh))


def percent_change(current, previous):
    return float(100*abs(current-previous)/abs(current)) if current else None


def nodal_rigid_virtual_work(beam, load):
    """Check consistent load under 3 rigid translations and 3 rotations.

    For 6-DOF rotation-aware elements, spin=angular velocity. For conventional
    Cauchy Hex8 there is no rotational DOF. The vector describes a virtual
    rigid displacement delta_u=omega cross X, delta_phi=omega.
    """
    xyz=np.asarray(beam.grid_points)
    F=np.asarray(load)
    translations=[]
    moments=[]
    for k in range(3):
        v=np.zeros((len(xyz),beam.ndof_per_node))
        v[:,k]=1.
        translations.append(float(F@v.ravel()))
        v[:,:3]=np.cross(np.broadcast_to(np.eye(3)[k],xyz.shape),xyz)
        if beam.ndof_per_node==6 and beam.__class__.__name__ != 'AffineEnrichedHex8':
            v[:,3+k]=1.
        moments.append(float(F@v.ravel()))
    return translations,moments


def sample_interior_stress(beam, U):
    """Sample all meshes at one shared, nonsingular Cartesian point (psi).

    This is a mesh-comparison QoI, NOT a stress-error estimate. It avoids
    element faces in the benchmark mesh families but is not a peak stress.
    """
    location=np.array([0.125*beam.Lx,0.375*beam.Ly,0.55*beam.Lz])
    for e in range(beam.get_num_elements()):
        coords=beam.get_element_nodes(e)
        lo=coords.min(axis=0);hi=coords.max(axis=0)
        if np.all(location>=lo-1e-9) and np.all(location<=hi+1e-9):
            p=2*(location-lo)/(hi-lo)-1
            stress=stress_at_element(beam,e,U,p)
            xx,yy,zz=stress.diagonal()
            vm=math.sqrt(max(0.,0.5*((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)
                             +3*(stress[0,1]**2+stress[1,2]**2+stress[0,2]**2)))
            return {'sample_location_in':location.tolist(),
                    'section_sigma_zz_psi':float(zz),
                    'section_von_mises_psi':float(vm),
                    'section_stress_tensor_psi':stress.tolist()}
    raise AssertionError('Physical stress sample lies outside mesh')


def cantilever_case(formulation,mesh,indicator=True,write_vtk_path=None):
    beam=build_beam(REFERENCE_LONGITUDINAL_E,
                    nx=mesh[0]+1,ny=mesh[1]+1,nz=mesh[2]+1,
                    formulation=formulation)
    load,_,traction=build_end_face_traction_load_vector(beam)
    resultants,moments=nodal_rigid_virtual_work(beam,load)
    target_F=np.array([0.,float(P),0.])
    target_M=np.array([-beam.Lz*float(P),0.,float(beam.grid_points[:,0].mean()*P)])
    if not np.allclose(resultants,target_F,rtol=1e-9,atol=1e-8) or not np.allclose(moments,target_M,rtol=1e-9,atol=1e-8):
        raise ArithmeticError(f'Consistent load fails virtual-work check: F={resultants}, M={moments}')
    U=solve_with_force_vector(beam,load)
    tip=float(free_end_uy(beam,U))
    endnodes=beam.get_nodes_at_max_z()
    uy=U.reshape((-1,beam.ndof_per_node))[endnodes,1]
    energy_sq=float(U@(beam.build_global_K()@U))
    if energy_sq<0:
        raise ArithmeticError('Negative stiffness energy')
    result={
        'study':'cantilever_mesh_change_not_true_error',
        'formulation':formulation,'mesh_elements_xyz':list(mesh),
        'elements':int(beam.get_num_elements()),'dofs':int(U.size),
        'tip_uy_in':tip,'end_uy_spread_in':float(np.ptp(uy)),
        'discrete_strain_energy_lbf_in':0.5*energy_sq,
        'external_work_lbf_in':float(load@U),
        'work_energy_relative_mismatch':float(abs(energy_sq-load@U)/max(abs(energy_sq),1e-30)),
        'load_resultant_lbf':resultants,'load_moment_lbf_in':moments,
        'expected_resultant_lbf':target_F.tolist(),
        'expected_moment_lbf_in':target_M.tolist(),
        'traction_psi':float(traction),'fixed_modulus_Ez_psi':float(REFERENCE_LONGITUDINAL_E),
        'clamp':'strong_full_face' if formulation in NONLOCAL else 'fixed_clamp_nodes',
        'true_error_available':False,
    }
    result.update(sample_interior_stress(beam,U))
    if hasattr(beam,'quadratic_rank'):
        result['quadratic_patch_ranks']={str(v):int(sum(rank==v for rank in beam.quadratic_rank)) for v in set(beam.quadratic_rank)}
    if indicator:
        eta=compute_residual_indicator(beam,U,
                      end_traction=(0.,traction,0.),length_scale='geometric_mean')
        result.update({
            'experimental_indicator_percent':eta['estimated_energy_error_percent_experimental'],
            'experimental_indicator_norm_sqrt_lbf_in':eta['eta_energy_units_sqrt_lbf_in'],
            'volume_residual_squared_lbf_in':eta['volume_residual_squared_lbf_in'],
            'interior_jump_squared_lbf_in':eta['interior_traction_jump_squared_lbf_in'],
            'natural_mismatch_squared_lbf_in':eta['natural_traction_mismatch_squared_lbf_in'],
        })
        if write_vtk_path is not None:
            from aposteriori.vtk_export import write_indicator_vtk
            write_vtk_path.parent.mkdir(parents=True,exist_ok=True)
            write_indicator_vtk(beam,U,eta,write_vtk_path)
    else:
        result['experimental_indicator_percent']=None
    return result


def annotate_changes(rows, key, output_key):
    prev={}
    for r in rows:
        identifier=r['formulation']
        val=r.get(key)
        r[output_key]=(percent_change(val,prev[identifier])
                        if identifier in prev and val is not None else None)
        if val is not None:
            prev[identifier]=val


def csv_export(path,rows):
    if not rows:
        return
    columns=list(dict.fromkeys(k for row in rows for k,v in row.items() if not isinstance(v,(list,dict))))
    with path.open('w',newline='',encoding='utf-8') as fp:
        out=csv.DictWriter(fp,fieldnames=columns,extrasaction='ignore')
        out.writeheader()
        for row in rows:
            out.writerow({k:row.get(k) for k in columns})


def fmt(value,unit=''):
    if value is None:
        return 'N/A'
    return f'{value:.6g}{unit}'


def markdown_report(bending,cantilever):
    report=[
        '# Session 5 — Hermite refinement and reference diagnostics','',
        '**RESEARCH ONLY — residual indicators are not certified error bounds.**','',
        '## Exact analytical orthotropic bending','',
        'These are actual errors against an equilibrium-satisfying 3D manufactured solution.',
        'The exact displacement, strain and stress fields are integrated independently.',
        'Near-exact results cannot produce a stable or meaningful effectivity index.','',
        '| Formulation | Elements | True energy error | Stress L2 error | Corner uz error | Experimental η | Effectivity |',
        '|---|---|---:|---:|---:|---:|---:|'
    ]
    for r in bending:
        report.append(f'| {r["formulation"]} | {mesh_name(r["mesh_elements_xyz"])} | '
                      f'{fmt(r["true_energy_error_percent"],"%")} | '
                      f'{fmt(r["stress_L2_relative_error_percent"],"%")} | '
                      f'{fmt(r["corner_uz_absolute_error_in"]," in")} | '
                      f'{fmt(r["experimental_indicator_percent"],"% ")} | '
                      f'{fmt(r["experimental_effectivity"])} |')
    report.extend(['','## Applied-traction cantilever — *mesh change only*','',
                   'This loading has no supplied exact three-dimensional displacement',
                   'or stress solution. Mesh-change percentages below are **not true errors**.',
                   'The modulus remains fixed; no displacement calibration is performed.',
                   'The Cauchy reference formulation is an independent FE comparator,',
                   '**not an exact continuum reference**.','',
                   '| Formulation | Elements | Tip uy | Tip mesh change | σzz sample | Stress mesh change | Residual η | Energy/work mismatch |',
                   '|---|---|---:|---:|---:|---:|---:|---:|'])
    for r in cantilever:
        report.append(f'| {r["formulation"]} | {mesh_name(r["mesh_elements_xyz"])} | '
                      f'{fmt(r["tip_uy_in"]," in")} | '
                      f'{fmt(r.get("previous_mesh_tip_change_percent"),"% ")} | '
                      f'{fmt(r["section_sigma_zz_psi"]," psi")} | '
                      f'{fmt(r.get("previous_mesh_stress_change_percent"),"% ")} | '
                      f'{fmt(r["experimental_indicator_percent"],"% ")} | '
                      f'{fmt(r["work_energy_relative_mismatch"])} |')
    report.extend(['','## Interpretation and limitations','',
        '- Exactly reproduced quadratic fields test polynomial completeness, not convergence for arbitrary loading.',
        '- If error is at numerical roundoff, do not interpret experimental effectivity.',
        '- Strong full-face constraints are required for the nonlocal rotational formulation.',
        '- Linear recovery on some grids cannot exactly enforce nonaffine Dirichlet traces.',
        '- Finite element differences and Euler–Bernoulli displacement differences are not certified errors.',
        '- Residual scale (`geometric_mean`) is experimental and does not have a reliability proof.',
        '- Full 3D equilibrium and boundary-condition validation remain open engineering requirements.',
        '- Cantilever stress is sampled at (x,y,z)=(0.125 Lx, 0.375 Ly, 0.55 Lz). This is NOT a peak stress.',
        '- If a quadratic recovery patch is not rank 9, the element falls back to linear recovery.',
        '- Units: in, lbf, psi, lbf·in; energy-norm quantities sqrt(lbf·in).',''])
    return '\n'.join(report)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bending-meshes',type=parse_mesh,nargs='+',default=[(2,2,4),(2,2,6)])
    p.add_argument('--cantilever-meshes',type=parse_mesh,nargs='+',default=[(2,2,4),(2,2,6)])
    p.add_argument('--bending-formulations',nargs='+',choices=FORMULATIONS,
                   default=['quadratic_rotation_hermite','cauchy_hex8_reference'])
    p.add_argument('--cantilever-formulations',nargs='+',choices=FORMULATIONS,
                   default=['quadratic_rotation_hermite','cauchy_hex8_reference'])
    p.add_argument('--skip-indicators',action='store_true',help='Faster: no residual diagnostics/VTK')
    p.add_argument('--max-indicator-elements',type=int,default=16,
                   help='Skip expensive indicators above this element count (default: 16; 0 to disable limit)')
    p.add_argument('--vtk',action='store_true',help='Write ParaView per-element indicator fields')
    p.add_argument('--quick',action='store_true',help='Run only smallest mesh from each study')
    p.add_argument('--out',type=Path,default=Path('outputs/session5_refinement'))
    args=p.parse_args(argv)
    args.out.mkdir(parents=True,exist_ok=True)
    bm=args.bending_meshes[:1] if args.quick else args.bending_meshes
    cm=args.cantilever_meshes[:1] if args.quick else args.cantilever_meshes
    bending=[];cantilever=[]
    for mesh in bm:
        for form in args.bending_formulations:
            # The quadratic trace is only exactly representable on sufficiently
            # wide meshes; do not pretend rank-deficient patches are complete.
            strong=form=='quadratic_rotation_hermite' and min(mesh)>=2
            if form=='quadratic_rotation_hermite' and not strong:
                print('WARNING: quadratic recovery has incomplete patches on '
                      f'{mesh_name(mesh)}; using nodal clamp only.',flush=True)
            print(f'[bending] {form} | {mesh_name(mesh)} | full-face={strong}',flush=True)
            indicator=not args.skip_indicators and (args.max_indicator_elements==0 or math.prod(mesh)<=args.max_indicator_elements)
            if not indicator and not args.skip_indicators:
                print('  indicator skipped: exceeds --max-indicator-elements',flush=True)
            r=evaluate(form,mesh,strong_face_dirichlet=strong,
                       compute_indicator=indicator)
            bending.append(r)
            print('  true energy',fmt(r['true_energy_error_percent'],'%'),
                  '| stress L2',fmt(r['stress_L2_relative_error_percent'],'%'),
                  '| corner uz error',fmt(r['corner_uz_absolute_error_in'],' in'),
                  '| η',fmt(r['experimental_indicator_percent'],'%'),flush=True)
    for mesh in cm:
        for form in args.cantilever_formulations:
            print(f'[cantilever] {form} | {mesh_name(mesh)}',flush=True)
            indicator=not args.skip_indicators and (args.max_indicator_elements==0 or math.prod(mesh)<=args.max_indicator_elements)
            if not indicator and not args.skip_indicators:
                print('  indicator skipped: exceeds --max-indicator-elements',flush=True)
            vtk=(args.out/f'{form}_{mesh_name(mesh)}_indicators.vtk') if args.vtk and indicator else None
            r=cantilever_case(form,mesh,indicator=indicator,
                             write_vtk_path=vtk)
            cantilever.append(r)
            print('  tip uy',fmt(r['tip_uy_in'],' in'),
                  '| η',fmt(r['experimental_indicator_percent'],'%'),
                  '| energy/work mismatch',fmt(r['work_energy_relative_mismatch']),flush=True)
    annotate_changes(cantilever,'tip_uy_in','previous_mesh_tip_change_percent')
    annotate_changes(cantilever,'section_sigma_zz_psi','previous_mesh_stress_change_percent')
    result={'status':'EXPERIMENTAL_NOT_CERTIFIED',
            'study_units':{'length':'in','force':'lbf','stress':'psi',
                           'work':'lbf in','energy_norm':'sqrt(lbf in)'},
            'bending':bending,'cantilever':cantilever}
    (args.out/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    csv_export(args.out/'bending.csv',bending)
    csv_export(args.out/'cantilever.csv',cantilever)
    (args.out/'REPORT.md').write_text(markdown_report(bending,cantilever),encoding='utf-8')
    print(f'Wrote Session 5 results to {args.out.resolve()}',flush=True)
    return result


if __name__=='__main__':
    main()
