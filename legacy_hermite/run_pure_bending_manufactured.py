"""Exact non-affine Saint-Venant bending benchmark for orthotropic Cauchy elasticity.

Exact equilibrium, zero body load: sigma_zz = stress_slope*x; other stresses
are zero. Compliance gives the full 3D quadratic compatible displacement.
The z=0 end receives nonzero exact Dirichlet translations; top z=L receives
varying traction (0,0,stress_slope*x); all side faces are traction free.
This is a verification problem, NOT the end-loaded calibrated cantilever.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from solver import build_beam, solve_with_prescribed_displacement
from aposteriori.residual import compute_residual_indicator


def exact_bending_fields(beam, stress_slope):
    S = np.linalg.inv(beam.D)
    beta, gamma, alpha = stress_slope * S[:3,2]

    def displacement(xyz):
        x,y,z = np.asarray(xyz)
        return np.array([0.5*beta*x*x - 0.5*gamma*y*y - 0.5*alpha*z*z,
                         gamma*x*y, alpha*x*z])

    def strain(xyz):
        return stress_slope * float(xyz[0]) * S[:,2]

    def traction(xyz, normal):
        return np.array([0.,0.,stress_slope*xyz[0]*normal[2]])

    return displacement,strain,traction


def evaluate(formulation='cauchy_hex8_reference', mesh=(1,1,2),stress_slope=30.,
             length_scale='geometric_mean',strong_face_dirichlet=False,compute_indicator=True):
    beam=build_beam(nx=mesh[0]+1,ny=mesh[1]+1,nz=mesh[2]+1,
                    formulation=formulation)
    exact_u,exact_e,exact_t = exact_bending_fields(beam,stress_slope)
    # Natural boundary is only the z=L face: all other natural tractions vanish.
    load=np.zeros(len(beam.grid_points)*beam.ndof_per_node)
    for e in beam.get_elements_at_max_z():
        coords=beam.get_element_nodes(e)
        is_nonlocal = hasattr(beam, "element_transform")
        dofs = (beam.element_transform(e)[0] if is_nonlocal
                else beam.get_element_global_dofs(e).ravel())
        for i,xi in enumerate(beam.gauss_points):
            for j,eta in enumerate(beam.gauss_points):
                p=(xi,eta,1.)
                J=beam.get_hermite_jacobian(*p,coords)
                area=np.linalg.norm(np.cross(J[0],J[1]))
                x=beam.hex8_shape_functions(*p) @ coords
                if is_nonlocal:
                    _,Ndisp=beam.displacement_operator(e,p)
                else:
                    Ndisp,*_=beam.get_hermite_displacement_matrices(*p,coords)
                load[dofs] += Ndisp.T @ exact_t(x,np.array([0.,0.,1.])) * area * beam.gauss_weights[i]*beam.gauss_weights[j]
    S = np.linalg.inv(beam.D)
    gamma = stress_slope * S[1,2]
    alpha = stress_slope * S[2,2]
    # The legacy beam imposes essential microrotations at the clamp. Match the
    # exact continuum curl there rather than incorrectly setting it to zero.
    def exact_spin(xyz):
        x,y,z=np.asarray(xyz)
        return np.array([0.,-alpha*z,gamma*y])
    face_constraint_stats=None
    if strong_face_dirichlet:
        from face_dirichlet import solve_with_face_dirichlet
        U,face_constraint_stats=solve_with_face_dirichlet(
            beam,load,exact_u,exact_spin)
    else:
        U=solve_with_prescribed_displacement(beam,load,prescribed_translations=exact_u,
                                             prescribed_rotations=exact_spin)
    true_sq=exact_sq=0.
    stress_difference_sq=stress_exact_sq=0.
    maximum_zz_stress_error_psi=0.
    # Independent 5x5x5 quadrature, distinct from 3/4-point FE assembly.
    gps,wts=np.polynomial.legendre.leggauss(5)
    for e in range(beam.get_num_elements()):
        coords=beam.get_element_nodes(e)
        ue=U[beam.get_element_global_dofs(e).ravel()]
        for i,xi in enumerate(gps):
            for j,eta in enumerate(gps):
                for k,zeta in enumerate(gps):
                    p=(xi,eta,zeta)
                    if is_nonlocal:
                        detJ=np.linalg.det(beam.get_hermite_jacobian(*p,coords))
                        eh=beam.strain_at(e,U,p)
                    else:
                        B,detJ,*_=beam.build_B_matrix(*p,coords)
                        eh=B@ue
                    x=beam.hex8_shape_functions(*p) @ coords
                    e_exact=exact_e(x)
                    diff=eh-e_exact
                    weight=detJ*wts[i]*wts[j]*wts[k]
                    true_sq += float(diff@beam.D@diff)*weight
                    exact_sq += float(e_exact@beam.D@e_exact)*weight
                    sigma_diff=beam.D@diff
                    sigma_exact=beam.D@e_exact
                    stress_difference_sq += float(sigma_diff@sigma_diff)*weight
                    stress_exact_sq += float(sigma_exact@sigma_exact)*weight
                    maximum_zz_stress_error_psi=max(maximum_zz_stress_error_psi,abs(float(sigma_diff[2])))
    residual=(compute_residual_indicator(beam,U,boundary_traction=exact_t,
                                         length_scale=length_scale)
              if compute_indicator else None)
    actual=float(100*np.sqrt(true_sq/exact_sq))
    # The upper x/y/z corner is a physical nodal displacement QoI. It is
    # independently computable from the manufactured quadratic field.
    corner=np.argmax(beam.grid_points[:,0]+10*beam.grid_points[:,1]+100*beam.grid_points[:,2])
    corner_xyz=beam.grid_points[corner]
    predicted_corner_uz=float(U[beam.global_dofs[corner,2]])
    exact_corner_uz=float(exact_u(corner_xyz)[2])
    eta=(float(residual['estimated_energy_error_percent_experimental']) if residual else None)
    actual_norm=float(np.sqrt(true_sq))
    indicator_norm=(float(residual['eta_energy_units_sqrt_lbf_in']) if residual else None)
    return {
        'formulation':formulation,'mesh_elements_xyz':list(mesh),
        'stress_slope_psi_per_in':stress_slope,
        'true_energy_error_percent':actual,
        'true_error_norm_sqrt_lbf_in':actual_norm,
        'exact_solution_norm_sqrt_lbf_in':float(np.sqrt(exact_sq)),
        'experimental_indicator_percent':eta,
        'experimental_indicator_norm_sqrt_lbf_in':indicator_norm,
        'experimental_effectivity':(indicator_norm/actual_norm
             if residual and actual_norm>1e-6*np.sqrt(exact_sq) else None),
        'effectivity_interpretable':bool(actual_norm>1e-6*np.sqrt(exact_sq)),
        'tip_uy_in':float(U.reshape((-1,beam.ndof_per_node))[:,1].max()),
        'corner_uz_in':predicted_corner_uz,
        'exact_corner_uz_in':exact_corner_uz,
        'corner_uz_absolute_error_in':abs(predicted_corner_uz-exact_corner_uz),
        'corner_uz_relative_error_percent':(100*abs(predicted_corner_uz-exact_corner_uz)/abs(exact_corner_uz)
                                             if exact_corner_uz != 0 else None),
        'stress_L2_relative_error_percent':(100*np.sqrt(stress_difference_sq/stress_exact_sq)
                                            if stress_exact_sq > 0 else None),
        'maximum_sampled_sigma_zz_error_psi':maximum_zz_stress_error_psi,
        'quadratic_patch_ranks':({str(v):int(sum(r==v for r in beam.quadratic_rank)) for v in set(beam.quadratic_rank)}
                                  if hasattr(beam,'quadratic_rank') else None),
        'length_scale':length_scale,
        'volume_residual_sq_lbf_in':residual['volume_residual_squared_lbf_in'] if residual else None,
        'face_jump_sq_lbf_in':residual['interior_traction_jump_squared_lbf_in'] if residual else None,
        'natural_mismatch_sq_lbf_in':residual['natural_traction_mismatch_squared_lbf_in'] if residual else None,
        'strong_face_dirichlet':strong_face_dirichlet,
        'face_constraint_stats':face_constraint_stats,
        'status':'TRUE_BENDING_ERROR_KNOWN_ESTIMATOR_UNVALIDATED',
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--meshes',nargs='*',default=['1x1x2','1x2x4'])
    parser.add_argument('--formulations',nargs='*',default=['legacy_hermite','cauchy_hex8_reference','affine_enriched_hex8'])
    parser.add_argument('--length-scale',choices=['diameter','geometric_mean','min_edge'],default='geometric_mean')
    parser.add_argument('--strong-face-dirichlet',action='store_true',help='KKT trace enforcement; nonlocal research elements only')
    parser.add_argument('--skip-indicator',action='store_true',help='Skip expensive strong-residual calculation; report true error only')
    parser.add_argument('--out',type=Path,default=Path('outputs/pure_bending_manufactured'))
    args=parser.parse_args()
    rows=[]
    for m in args.meshes:
        mesh=tuple(int(s) for s in m.lower().split('x'))
        if len(mesh)!=3 or min(mesh)<1:
            parser.error('Meshes must have three positive element counts, e.g. 1x2x4')
        for name in args.formulations:
            print(f'[{name} | {m}] solving exact bending benchmark...',flush=True)
            if args.strong_face_dirichlet and name not in ('recovered_rotation_hermite','quadratic_rotation_hermite'):
                parser.error('--strong-face-dirichlet requires a nonlocal rotation Hermite formulation')
            r=evaluate(name,mesh,length_scale=args.length_scale,
                       strong_face_dirichlet=args.strong_face_dirichlet,
                       compute_indicator=not args.skip_indicator)
            rows.append(r)
            print(f'  TRUE energy error={r["true_energy_error_percent"]:.6g}% '
                  f' experimental indicator={(format(r["experimental_indicator_percent"],".6g")+"%" if r["experimental_indicator_percent"] is not None else "skipped")} ' 
                  f' effectivity={(format(r["experimental_effectivity"],".6g") if r["experimental_effectivity"] is not None else "N/A (near exact)")}',flush=True)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'results.json').write_text(json.dumps({'rows':rows},indent=2)+'\n',encoding='utf-8')
    report=['# Exact orthotropic pure-bending benchmark','',
            '**Non-affine manufactured continuum solution, not the original cantilever.**','',
            'sigma_zz = (30 psi/in) * x; sigma_xx = sigma_yy = all shear stresses = 0.',
            'Cauchy equilibrium: div sigma = 0; traction varies across z=L face.',
            'The exact nonzero displacement is prescribed at z=0.',
            'All errors are energy errors; the experimental indicator is NOT a bound.','',
            '| Formulation | Mesh | TRUE energy error | Experimental eta | Effectivity |',
            '|---|---|---:|---:|---:|']
    for r in rows:
        eta = (f'{r["experimental_indicator_percent"]:.6g}%'
               if r['experimental_indicator_percent'] is not None else 'Not computed')
        effectivity = (f'{r["experimental_effectivity"]:.6g}'
                      if r['experimental_effectivity'] is not None else 'N/A')
        report.append(f'| {r["formulation"]} | {"x".join(map(str,r["mesh_elements_xyz"]))} | '
                      f'{r["true_energy_error_percent"]:.6g}% | {eta} | {effectivity} |')
    report+=['','**WARNING:** The enriched six-DOF comparator uses interior bubble',
             'amplitudes, not physical Hermite nodal rotations. The new nonlocal',
             'recovered-rotation branches use true nodal spin variables but a wider',
             'gradient-recovery stencil; strong face Dirichlet constraints are needed.',
             'For near-exact manufactured solutions, effectivity is numerically',
             'indeterminate and reported N/A. The original solver is preserved.',
             'These tests do not certify a posteriori estimates.']
    (args.out/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    print(f'Saved {args.out.resolve()}',flush=True)

if __name__=='__main__':main()
