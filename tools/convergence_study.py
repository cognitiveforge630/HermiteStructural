"""Reproducible Hex27/Hex64: mesh, solve, report, view, compare.

Run `python -m tools.convergence_study --help` from the repository root.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from core.mesh.nodal_hex import StructuredNodalHexMesh
from core.reference.beam_theory import B,H,L,E,I,P
from core.reference.awc import DEFAULT_WEIGHT_DENSITY
from tools.p_order_study import build_plan, patch_checks, element_verification, solve_study
from tools.nodal_field_viewer import require_inset_archive


def _parse_levels(value):
    levels=[]
    for item in value.split(';'):
        vals=tuple(int(v.strip()) for v in item.split(','))
        if len(vals)!=3 or min(vals)<2:
            raise argparse.ArgumentTypeError('Levels must look like 2,2,4;3,3,6 (vertex counts >=2)')
        levels.append(vals)
    return levels


def _archive_path(directory, order, mesh):
    return directory / f"hex{(order+1)**3}_{'x'.join(map(str,mesh))}"


def beam_displacement(z, density, load):
    """Classical Euler-Bernoulli y deflection along z, including uniform weight."""
    q=density*B*H
    return load*z*z*(3*L-z)/(6*E*I) - q*z*z*(6*L*L-4*L*z+z*z)/(24*E*I)


def plot_mesh(order, levels, directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=directory/'mesh_density.png'
    fig, axs=plt.subplots(len(levels), 1, figsize=(11, max(3, 2.7*len(levels))), squeeze=False)
    for j, counts in enumerate(levels):
        div=tuple(v-1 for v in counts)
        ax=axs[j,0]
        # A longitudinal XY-profile: all element edges on the x=0 cut in physical coordinates.
        yy=np.linspace(-H/2,H/2,div[1]+1)
        zz=np.linspace(0,L,div[2]+1)
        for y in yy: ax.plot([0,L],[y,y],color='#234d7b',lw=.65)
        for z in zz: ax.plot([z,z],[-H/2,H/2],color='#234d7b',lw=.65)
        ax.set(xlabel='Beam axis z (in)',ylabel='Depth y (in)',
               title=f'Hex{(order+1)**3}: {div[0]} × {div[1]} × {div[2]} physical elements; '
                     f'{math.prod(div)} cells; {build_plan(order,*counts)["free_and_fixed_translation_DOFs"]:,} DOFs')
        ax.set_aspect('auto')
    fig.tight_layout(); fig.savefig(out,dpi=160);plt.close(fig)
    return out


def plot_fields(archive, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with np.load(archive) as d:
        require_inset_archive(d)
        pos=d['coordinates']; u=d['displacement']
        # Raw samples at inset visualization-subcell corners, NEVER centroids.
        corners=d['stress_sample_coordinates']
        sig=d['fe_stress']; ref=d['beam_reference_stress']
        vm=d['fe_von_mises']; rvm=d['beam_reference_von_mises']
        density=float(d['density']); load=float(d['end_force'])
        order=int(d['order']); div=tuple(int(v) for v in d['divisions'])
    # One consistent one-sided inset corner per visualization subcell. The
    # 3D viewer independently takes sample extrema over all eight corners.
    samples=corners[:,0,:]
    cx=np.unique(np.round(samples[:,0], 9)); chosen=cx[np.argmin(abs(cx))]
    mask=np.isclose(samples[:,0],chosen)
    c=samples[mask]; s=sig[mask,0,:]; r=ref[mask,0,:]
    fv=vm[mask,0]; rv=rvm[mask,0]
    # Nodal y-z slice near x=0, plotted with a regular grid.
    ux=np.unique(np.round(pos[:,0],9)); cut=ux[np.argmin(abs(ux))]
    nod=np.isclose(pos[:,0],cut); xyz=pos[nod]; un=u[nod]
    fields=[('FE displacement magnitude (in)',xyz[:,2],xyz[:,1],np.linalg.norm(un,axis=1)),
        ('Beam-reference |uy| (in)',xyz[:,2],xyz[:,1],abs(beam_displacement(xyz[:,2],density,load))),
        ('FE tension max(sigma_zz,0) (psi)',c[:,2],c[:,1],np.maximum(s[:,2],0)),
        ('Beam-reference tension (psi)',c[:,2],c[:,1],np.maximum(r[:,2],0)),
        ('FE compression max(-sigma_zz,0) (psi)',c[:,2],c[:,1],np.maximum(-s[:,2],0)),
        ('Beam-reference compression (psi)',c[:,2],c[:,1],np.maximum(-r[:,2],0)),
        ('FE von Mises (psi)',c[:,2],c[:,1],fv),
        ('Beam-reference von Mises (psi)',c[:,2],c[:,1],rv),
        ('FE |tau_yz| (psi)',c[:,2],c[:,1],np.abs(s[:,4])),
        ('Beam-reference |tau_yz| (psi)',c[:,2],c[:,1],np.abs(r[:,4]))]
    fig,axes=plt.subplots(5,2,figsize=(15,16),layout='constrained')
    for row in range(5):
        vmax=max(float(np.max(fields[row*2][3])),float(np.max(fields[row*2+1][3])),1e-20)
        for col in range(2):
            title,z,y,values=fields[row*2+col]; ax=axes[row,col]
            # One-sided inset samples and nodal translations on structured slices.
            zz=np.unique(z); yy=np.unique(y)
            grid=np.full((len(yy),len(zz)),np.nan)
            iz=np.searchsorted(zz,z); iy=np.searchsorted(yy,y)
            grid[iy,iz]=values
            dz = float(np.median(np.diff(zz))) if len(zz)>1 else L/max(div[2],1)
            dy = float(np.median(np.diff(yy))) if len(yy)>1 else H/max(div[1],1)
            extent=(float(min(zz)-dz/2),float(max(zz)+dz/2),
                    float(min(yy)-dy/2),float(max(yy)+dy/2))
            im=ax.imshow(grid, origin='lower',extent=extent,aspect='auto',vmin=0,vmax=vmax,cmap='viridis',interpolation='nearest')
            ax.set(title=title,xlabel='z (in)',ylabel='y (in)')
            fig.colorbar(im,ax=ax,shrink=.72)
    banner=('SMOKE TEST ONLY — STRESS NOT VALIDATED / NOT FOR DESIGN'
            if np.prod(div)<=16 else 'RESEARCH ONLY — STRESS NOT CERTIFIED')
    fig.suptitle(f'Hex{(order+1)**3}, divisions {div}: FE versus approximate classical beam reference\n'
                 + banner + '\nOne directly evaluated inset corner per pixel; NEVER centroid stress',fontsize=11)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=150);plt.close(fig)
    return output


def report_text(rows, order):
    lines=[f'# Hex{(order+1)**3} cantilever mesh convergence study','',
        '## Interpretation',
        'Each row uses the same material, geometry, gravity, end traction and translation-only z=0 clamp.',
        'The AWC/Euler-Bernoulli reference is one-dimensional beam theory, **not** an exact three-dimensional solution.',
        'Neighboring-mesh changes measure numerical sensitivity, **not** a certified discretization-error bound.',
        'Stress extrema below come from **direct inset visualization-subcell corner samples**, never centroids; not guaranteed global maximum stresses.',
        'SMOKE TEST (16 elements or fewer): reported stress magnitudes are illustrative ONLY; not validated for design or final comparison.',
        'Tension and compression refer to the longitudinal normal stress sigma_zz; von Mises is nonnegative.',
        'Clamped-end peaks and corner singularities require special care. A smaller tip discrepancy alone does not prove stress convergence.','',
        '## Results','',
        '| Vertex mesh | DOFs | Tip FE uy (in) | Beam uy (in) | Tip vs beam (%) | Change vs prior (%) | Max sampled section-shear discrepancy (%) | Max sampled section-moment discrepancy (%) | Sample max tension (psi) | Sample max compression (psi) | Sample max von Mises (psi) |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['status']!='solved':
            lines.append(f"| {r['mesh']} | {r['dofs']} | {r['status']} | — | — | — | — | — | — | — | — | — |")
            continue
        def f(k): return f"{r[k]:.7g}" if isinstance(r.get(k),(float,int)) else '—'
        lines.append(f"| {r['mesh']} | {r['dofs']:,} | {f('tip')} | {f('beam')} | {f('tip_beam_difference_pct')} | {f('change_previous_pct')} | {f('shear')} | {f('moment')} | {f('tension')} | {f('compression')} | {f('von_mises')} |")
    lines+=['','## What this establishes',
        'Polynomial and single-element mechanical checks are executed before each solve.',
        'Global equilibrium and energy diagnostics are retained in each per-mesh JSON.',
        'Inspect section CSV data for shear and bending moment deviations from beam theory. Their reported maxima sample only the middle half of the beam span.',
        'An extrapolation, certified error estimator, and physically defensible stress QoI treatment are NOT provided.',
        '## Viewing','',
        'See `mesh_density.png` before solving and `*_fields.png` after solving.',
        'Each field figure pairs an inset-subcell-corner FE map with the beam-theory reference using a shared color scale. Maps include |tau_yz| shear.',
        'The reference beam displacement is only vertical; the FE displacement magnitude includes all three components. Shear maps show |tau_yz| rather than the section-integrated shear force.',
        'Run `python -m tools.convergence_study --stage view --order '+str(order)+' --levels '+ '"'+';'.join(r['mesh'].replace('x',',') for r in rows)+'"'+' --out outputs/convergence_hex'+str((order+1)**3)+'` to regenerate images.','']
    return '\n'.join(lines)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('mesh','solve','report','view','all'),default='mesh')
    parser.add_argument('--order',type=int,choices=(2,3),default=2)
    parser.add_argument('--levels',type=_parse_levels,default=_parse_levels('2,2,4;3,3,5;3,4,7'))
    parser.add_argument('--out',type=Path,default=None)
    parser.add_argument('--accept-expensive-run',action='store_true')
    parser.add_argument('--weight-density',type=float,default=DEFAULT_WEIGHT_DENSITY)
    parser.add_argument('--end-force',type=float,default=P)
    args=parser.parse_args(argv)
    dest=args.out or Path('outputs')/f'convergence_hex{(args.order+1)**3}'
    dest.mkdir(parents=True,exist_ok=True)
    if args.stage in ('mesh','all'):
        print('Mesh density:',plot_mesh(args.order,args.levels,dest),flush=True)
    if args.stage in ('solve','all'):
        if args.weight_density<=0 or args.end_force>0: parser.error('Positive weight density and nonpositive end force required')
        if not patch_checks(args.order)['passed'] or not element_verification(args.order)['passed']:
            raise AssertionError('Element verification failed; solve prohibited')
        for counts in args.levels:
            name=_archive_path(dest,args.order,counts); plan=build_plan(args.order,*counts)
            print(f'Case {counts}: {plan["free_and_fixed_translation_DOFs"]} DOFs',flush=True)
            if plan['free_and_fixed_translation_DOFs']>8000 and not args.accept_expensive_run:
                print('  SKIPPED: >8000 DOFs requires --accept-expensive-run (memory/time)',flush=True)
                continue
            res=solve_study(args.order,*counts,args.weight_density,args.end_force,
                            field_output_path=name.with_suffix('.npz'))
            res['polynomial_patch_tests']=patch_checks(args.order)
            res['one_cell_mechanics_tests']=element_verification(args.order)
            name.with_suffix('.json').write_text(json.dumps(res,indent=2)+'\n',encoding='utf-8')
            with name.with_name(name.name+'_sections.csv').open('w',newline='',encoding='utf-8') as f:
                writer=csv.DictWriter(f,fieldnames=list(res['source_sections'][0]))
                writer.writeheader();writer.writerows(res['source_sections'])
            print('  Solved; output:',name,flush=True)
    if args.stage in ('report','all'):
        rows=[]; prev=None
        for counts in args.levels:
            name=_archive_path(dest,args.order,counts); path=name.with_suffix('.json'); plan=build_plan(args.order,*counts)
            row={'mesh':'x'.join(map(str,counts)), 'dofs':plan['free_and_fixed_translation_DOFs']}
            if path.exists():
                res=json.loads(path.read_text(encoding='utf-8'))
                tip=res['tip_max_abs_uy_in']; beam=res['reference_tip_uy_in']
                stats=res['field_summary']
                if 'inset_subcell_corners_v2' not in stats.get('stress_sampling',''):
                    raise ValueError(f'OBSOLETE CENTER-SAMPLED REPORT: {path}; rerun --stage solve to regenerate numerical stresses')
                row.update(status='solved',tip=tip,beam=beam,
                    tip_beam_difference_pct=100*abs(tip-beam)/max(abs(beam),1e-30),
                    change_previous_pct=None if prev is None else 100*abs(tip-prev)/max(abs(tip),1e-30),
                    shear=res['max_midspan_shear_y_difference_pct'],
                    moment=res['max_midspan_moment_x_difference_pct'],
                    tension=stats['fe_max_tension_zz_psi'],compression=stats['fe_max_compression_zz_psi'],
                    von_mises=stats['fe_max_von_mises_psi'])
                prev=tip
            else: row['status']='not solved'
            rows.append(row)
        (dest/'report.md').write_text(report_text(rows,args.order),encoding='utf-8')
        (dest/'summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        print('Reports:',dest/'report.md',dest/'summary.json',flush=True)
    if args.stage in ('view','all'):
        for counts in args.levels:
            name=_archive_path(dest,args.order,counts)
            if name.with_suffix('.npz').exists():
                print('Heatmaps:',plot_fields(name.with_suffix('.npz'),name.with_name(name.name+'_fields.png')),flush=True)
            else: print('No solved field archive for',counts,flush=True)
    return dest

if __name__=='__main__':
    main()
