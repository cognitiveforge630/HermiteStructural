"""Object-oriented tensor-product Hex27/Hex64 cantilever research comparison.

Usage from repository root:
    python -m tools.p_order_study --mode plan --order 3 --nx 3 --ny 5 --nz 25
    python -m tools.p_order_study --mode verify --order 3
    python -m tools.p_order_study --mode solve --order 2 --nx 3 --ny 3 --nz 5

Unlike cubic_rotation_hermite, these are *conventional* C0 displacement-only
Lagrange hexahedra, with 3 translation DOFs per shared node. No certification.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import spsolve

from core.elements.nodal_hex import (
    element_type,
    strain_matrix,
    tensor_shapes,
)
from core.materials.redwood import redwood_orthotropic_material
from core.mesh.nodal_hex import StructuredNodalHexMesh, local_nodes
from core.reference.awc import DEFAULT_WEIGHT_DENSITY, awc_reference_components
from core.reference.beam_theory import B, E, H, L, P
from core.solvers.nodal_assembly import GlobalAssembler
from core.verification.nodal_hex import (
    ElementMechanicsVerifier,
    PolynomialPatchVerifier,
)


def build_plan(p, nx, ny, nz, quad_order=None):
    if min(nx, ny, nz) < 2:
        raise ValueError('nx/ny/nz refer to legacy Hermite vertex counts and must be >=2')
    if p not in (2, 3):
        raise ValueError('order must be 2 or 3')
    q = p + 1 if quad_order is None else int(quad_order)
    if q < p + 1:
        raise ValueError('Use >=p+1 Gauss points per axis; reduced integration is not supported')
    element_counts = (nx-1, ny-1, nz-1)
    node_counts = tuple(p*v + 1 for v in element_counts)
    num_elements = math.prod(element_counts)
    num_nodes = math.prod(node_counts)
    ndof = 3*num_nodes
    local_dof = 3*(p+1)**3
    dense_contributions = num_elements*local_dof**2
    return {
        'family':f'Hex{(p+1)**3}_tensor_Lagrange_C0',
        'polynomial_order':p,
        'original_Hermite_vertex_mesh':[nx,ny,nz],
        'physical_element_divisions':list(element_counts),
        'continuity':'C0', 'DOFs_per_node':3,
        'interpolation_nodes_per_element':(p+1)**3,
        'nodes_per_element_face':(p+1)**2,
        'actual_global_nodes_per_direction':list(node_counts),
        'global_nodes':num_nodes, 'free_and_fixed_translation_DOFs':ndof,
        'local_element_DOFs':local_dof,
        'total_dense_COO_contributions_before_coalescing':dense_contributions,
        'estimated_COO_triplets_MiB':round(24*dense_contributions/2**20,1),
        'note_memory':'Conceptual unbatched COO estimate; OOP implementation flushes batches of 8 elements into CSR. Sparse factorization may still use much more RAM.',
        'assembly_batch_elements':8,
        'gauss_order_per_axis':q,
        'constitutive_model':'same orthotropic Redwood reference parameters as the Hermite solver',
        'clamp':'all three translational DOFs fixed at z=0; no rotational DOFs',
        'loading':'uniform end traction on z=L; consistently integrated body self-weight',
        'certified_error_bound':False,
    }


def element_kernel(p, xyz_extent, D, density, traction, q):
    """Compatibility API implemented through an element's separate integrals.

    Returns the same (Ke, gravity Fe, +Z face traction Fe) as the old
    monolithic kernel, without embedding a face integral in a volume loop.
    """
    element = element_type(p)(xyz_extent, D, quadrature_order=q)
    ke = element.get_stiffness_matrix()
    fg = element.integrate_body_force((0., -density, 0.))
    element.clear_force_vector()
    ft = element.facets['+Z'].integrate_traction((0., traction, 0.))
    return ke, fg, ft


# Compatibility functions keep the established CLI and test imports stable.
def patch_checks(p, q=None):
    """Evaluate interpolation identities for a Hex27 or Hex64 element."""
    cls = element_type(p)
    element = cls((2., 3., 4.), np.eye(6), quadrature_order=q)
    return PolynomialPatchVerifier().verify(element)


def element_verification(p):
    """Evaluate the inexpensive single-cell mechanics tests."""
    cls = element_type(p)
    element = cls((2., 4., 12.), np.eye(6))
    return ElementMechanicsVerifier().verify(element)


def integrate_source_sections(p, div, dims, dx, U, D, density,
                              end_force, axial_positions, section_order):
    """Raw source-element stress integrated on actual XY Gaussian sections.

    This is a stress/resultant diagnostic, not a posteriori certification.
    The face quadrature used for *loading* is a separate element operation.
    """
    qpts,qweights=leggauss(section_order)
    section_rows=[]
    # Stress is computed from actual source element DOFs and Galerkin B, not viewer interpolation.
    for iz in range(div[2]):
        for zeta in axial_positions:
            area=0.;fx=0.;fy=0.;fz=0.;mx=0.
            z=(iz+(zeta+1)/2)*dx[2]
            for ix in range(div[0]):
                for iy in range(div[1]):
                    _,dofs=local_nodes(p,ix,iy,iz,dims)
                    uelement=U[dofs]
                    yc=-H/2+(iy+.5)*dx[1]
                    for a,xi in enumerate(qpts):
                        for b,eta in enumerate(qpts):
                            _,ds=tensor_shapes(p,(xi,eta,zeta))
                            strain=strain_matrix(ds/(dx/2)[:,None])@uelement
                            sig=D@strain
                            da=qweights[a]*qweights[b]*dx[0]*dx[1]/4
                            area+=da
                            fx+=sig[5]*da
                            fy+=sig[4]*da
                            fz+=sig[2]*da
                            mx+=(yc+eta*dx[1]/2)*sig[2]*da
            r=L-z
            beam_shear=end_force-density*B*H*r
            beam_moment=-end_force*r+0.5*density*B*H*r*r
            section_rows.append({'z_in':float(z),'source_zeta':float(zeta),
                'section_area_in2':area,'FE_shear_x_lbf':fx,
                'FE_shear_y_lbf':fy,'FE_axial_force_lbf':fz,
                'FE_bending_moment_x_lbf_in':mx,'beam_shear_y_lbf':beam_shear,
                'beam_bending_moment_x_lbf_in':beam_moment,
                'shear_y_reference_difference_pct':100*abs(fy-beam_shear)/max(abs(beam_shear),1e-30),
                'moment_x_reference_difference_pct':100*abs(mx-beam_moment)/max(abs(beam_moment),1e-30)})
    return section_rows


def build_equilibrium_report(K, U, F, coords, fixed, free, density, end_force):
    """Explain applied/reaction force, moment, residual, and work identities."""
    reaction = np.zeros_like(F)
    residual = K @ U - F
    reaction[fixed] = residual[fixed]
    rel_free = float(np.linalg.norm(residual[free]) / max(np.linalg.norm(F[free]), 1e-30))
    applied_forces=F.reshape(-1,3)
    reacted_forces=reaction.reshape(-1,3)
    total_applied=applied_forces.sum(axis=0)
    total_reaction=reacted_forces.sum(axis=0)
    total_applied_moment=np.cross(coords,applied_forces).sum(axis=0)
    total_reaction_moment=np.cross(coords,reacted_forces).sum(axis=0)
    expected_y=end_force-density*B*H*L
    expected_mx=-end_force*L+density*B*H*L*L/2
    reaction_force_defect=float(np.max(abs(total_reaction+total_applied)))
    reaction_moment_defect=float(np.max(abs(total_reaction_moment+total_applied_moment)))
    internal = 0.5 * float(U @ (K @ U))
    half_work = 0.5 * float(F @ U)
    return {
        'applied_force_xyz_lbf': total_applied.tolist(),
        'support_reaction_force_xyz_lbf': total_reaction.tolist(),
        'applied_moment_about_origin_xyz_lbf_in': total_applied_moment.tolist(),
        'support_reaction_moment_about_origin_xyz_lbf_in': total_reaction_moment.tolist(),
        'force_defect_inf_lbf': reaction_force_defect,
        'moment_defect_inf_lbf_in': reaction_moment_defect,
        'analytical_applied_Fy_lbf': expected_y,
        'analytical_applied_Mx_lbf_in': expected_mx,
        'analytical_resultant_force_Fy_error_lbf': float(total_applied[1]-expected_y),
        'analytical_resultant_Mx_error_lbf_in': float(total_applied_moment[0]-expected_mx),
        'free_DOF_relative_residual_L2': rel_free,
        'discrete_external_work_lbf_in': 2 * half_work,
        'discrete_FE_strain_energy_lbf_in': internal,
        'Clapeyron_work_energy_relative_defect': abs(internal-half_work)/max(abs(half_work),1e-30),
    }


def build_study_report(*, plan, material, density, end_force, reference,
                       tip_max, tip_mean, mid, K, fixed, equilibrium,
                       assembly_seconds, solve_seconds, total_seconds,
                       section_rows):
    """Build the backward-compatible human/JSON verification summary."""
    result={**plan,'research_not_certified':True,'comparison_is_not_a_posteriori_bound':True,
        'density_lbf_in3':density,'end_force_lbf':end_force,
        'longitudinal_E_psi':E,'orthotropic_constants':material.as_dict(),
        'reference_tip_uy_in':reference['total_uy_in'],
        'tip_max_abs_uy_in':tip_max,'tip_face_mean_uy_in':tip_mean,
        'tip_max_vs_classical_beam_pct':100*abs(tip_max-reference['total_uy_in'])/abs(reference['total_uy_in']),
        'max_midspan_shear_y_difference_pct':float(max(r['shear_y_reference_difference_pct'] for r in mid)),
        'max_midspan_moment_x_difference_pct':float(max(r['moment_x_reference_difference_pct'] for r in mid)),
        'stiffness_nnz':int(K.nnz),'clamped_DOF_count':int(fixed.size),
        'equilibrium':equilibrium,
        'runtime_seconds':{'assembly':assembly_seconds,'linear_solve':solve_seconds,
                           'total':total_seconds},
        'source_sections':section_rows,
        'caveat':('1D beam theory is not the exact 3D elasticity solution. '
            'The stress is raw constitutive Cauchy stress, not equilibrated; '
            'the clamped translation trace is not the Hermite fixed-spin boundary condition. '
            'No certified error estimate.')}
    return result


def solve_study(p,nx,ny,nz,density=DEFAULT_WEIGHT_DENSITY,end_force=P,
                q=None,axial_positions=(-0.9,0.9), section_order=None, field_output_path=None):
    plan=build_plan(p,nx,ny,nz,q)
    q=plan['gauss_order_per_axis']
    sorder=max(p+1,4) if section_order is None else section_order
    if sorder < p+1:
        raise ValueError('section quadrature should be at least p+1')
    for t in axial_positions:
        if not (-1 < t < 1):
            raise ValueError('Section locations must be strictly interior to source elements')
    t0=time.perf_counter()
    from core.materials.design_limits import reference_load_screen
    material=redwood_orthotropic_material(E)
    D=material.elastic_matrix()
    screen=reference_load_screen(tip_force_lbf=end_force, density_lbf_in3=density)
    if screen['status']!='PASS_PRELIMINARY':
        raise ValueError('Preliminary reference grade load screen did not pass')
    div=(nx-1,ny-1,nz-1)
    mesh = StructuredNodalHexMesh(p, div, (B, H, L), origin=(-B/2, -H/2, 0.))
    dims = mesh.node_shape
    coords = mesh.coordinates
    dx = mesh.element_extent
    traction = end_force / (B * H)

    # All cells here have identical affine geometry and D: share one dense Ke.
    reference_element = element_type(p)(dx, D, quadrature_order=q)
    ke = reference_element.get_stiffness_matrix()
    skew = float(np.max(abs(ke-ke.T)))
    if skew > 1e-8 * max(1., np.max(abs(ke))):
        raise ArithmeticError('Local stiffness is not symmetric')
    gravity_load = reference_element.integrate_body_force((0., -density, 0.))

    assembler = GlobalAssembler(mesh.global_dof_count, batch_elements=8)
    for (ix, iy, iz), element in mesh.iter_elements(
            D, quadrature_order=q, shared_stiffness=ke):
        # The body integral is the same for every congruent cell; reuse it.
        element.add_force_vector(gravity_load)
        # A named facet owns its own *2D* surface integral (not a volume loop).
        if iz == div[2] - 1:
            element.facets['+Z'].integrate_traction((0., traction, 0.))
        assembler.add_element(element)
        # The short-lived element can now be collected; Ke is shared.
    K, F = assembler.finalize()
    del assembler, reference_element, element, ke, gravity_load
    assembly_seconds=time.perf_counter()-t0
    # Essential displacement Dirichlet data on entire fixed face (C0 FE).
    fixed_nodes=np.flatnonzero(np.isclose(coords[:,2],0.))
    fixed = np.asarray((fixed_nodes[:, None] * 3 + np.arange(3)).ravel(), dtype=np.intp)
    free = np.asarray(np.setdiff1d(np.arange(F.size), fixed, assume_unique=True), dtype=np.intp)
    t1=time.perf_counter()
    uf = np.asarray(spsolve(K[free][:, free].tocsc(), F[free]), dtype=float).reshape(-1)
    U = np.zeros(F.size, dtype=float)
    U[free] = uf
    solve_seconds=time.perf_counter()-t1
    if not np.isfinite(U).all():
        raise ArithmeticError('Non-finite displacement solution')
    equilibrium = build_equilibrium_report(K, U, F, coords, fixed, free, density, end_force)
    nodal_end=np.flatnonzero(np.isclose(coords[:,2],L))
    tip_values=U.reshape(-1,3)[nodal_end,1]
    tip_max=float(tip_values[np.argmax(abs(tip_values))])
    tip_mean=float(tip_values.mean())
    reference=awc_reference_components(density,end_force=end_force)
    section_rows = integrate_source_sections(
        p, div, dims, dx, U, D, density, end_force, axial_positions, sorder)
    mid=[r for r in section_rows if L/4 <= r['z_in'] <=3*L/4]
    if not mid:
        raise ValueError('No sections between quarter and three-quarter span; refine nz')
    field_summary = None
    if field_output_path is not None:
        from core.visualization.nodal_convergence import export_fields
        field_summary = export_fields(field_output_path, mesh, U, D, density, end_force)
    report = build_study_report(
        plan=plan, material=material, density=density, end_force=end_force,
        reference=reference, tip_max=tip_max, tip_mean=tip_mean, mid=mid,
        K=K, fixed=fixed, equilibrium=equilibrium,
        assembly_seconds=assembly_seconds, solve_seconds=solve_seconds,
        total_seconds=time.perf_counter()-t0, section_rows=section_rows)
    if field_summary is not None:
        report['field_summary'] = field_summary
    return report


def baseline_compare(result, reference_json):
    d=json.loads(Path(reference_json).read_text(encoding='utf-8'))
    if d.get('formulation')!='cubic_rotation_hermite':
        raise ValueError('Expected cubic_rotation_hermite comparison JSON')
    for name, a, b in [('mesh',d.get('node_mesh'),result['original_Hermite_vertex_mesh']),
                       ('material Ez',d.get('longitudinal_E_psi'),result['longitudinal_E_psi']),
                       ('load',d.get('end_force_lbf'),result['end_force_lbf'])]:
        if name=='mesh':
            if list(a)!=list(b):
                raise ValueError(f'Mismatched {name}: {a} vs {b}')
        elif not np.isclose(a,b,rtol=1e-10,atol=1e-10):
            raise ValueError(f'Mismatched {name}: {a} vs {b}')
    awc = d.get('awc_reference')
    if isinstance(awc, dict):
        reference_weight = awc.get('total_weight_lbf')
        if reference_weight is not None:
            computed_weight = result['density_lbf_in3'] * B * H * L
            if not np.isclose(float(reference_weight), computed_weight):
                raise ValueError('Mismatched gravity load')
    return {'same_mesh_material_and_loads':True,
       'important_boundary_condition_difference':d['boundary_condition']+' versus displacement-only C0 clamp',
       'cubic_Hermite':{'tip_uy_in':d.get('fem_tip_uy_in'),
          'max_midspan_shear_y_difference_pct':d.get('cubic_source_max_midspan_shear_y_difference_pct'),
          'max_midspan_moment_x_difference_pct':d.get('cubic_source_max_midspan_moment_x_difference_pct')},
       'tensor_Lagrange':{'tip_uy_in':result['tip_max_abs_uy_in'],
          'max_midspan_shear_y_difference_pct':result['max_midspan_shear_y_difference_pct'],
          'max_midspan_moment_x_difference_pct':result['max_midspan_moment_x_difference_pct']}}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('plan','verify','solve'),default='plan')
    parser.add_argument('--order',type=int,choices=(2,3),default=2)
    parser.add_argument('--nx',type=int,default=3)
    parser.add_argument('--ny',type=int,default=5)
    parser.add_argument('--nz',type=int,default=25)
    parser.add_argument('--gauss-order',type=int,default=None)
    parser.add_argument('--section-order',type=int,default=None)
    parser.add_argument('--weight-density',type=float,default=DEFAULT_WEIGHT_DENSITY)
    parser.add_argument('--end-force',type=float,default=P)
    parser.add_argument('--baseline-cubic-json',type=Path)
    parser.add_argument('--accept-expensive-run',action='store_true',
        help='Explicitly accept heavy sparse factorization for systems >8000 DOFs')
    parser.add_argument('--out',type=Path,default=None)
    parser.add_argument('--field-output',type=Path,default=None,
        help='Optional .npz archive for mesh, displacement, and stress visualization')

    args=parser.parse_args(argv)

    plan=build_plan(args.order,args.nx,args.ny,args.nz,args.gauss_order)

    if args.mode=='plan':
        print(json.dumps(plan,indent=2))
        print('Plan only: no stiffness assembly or linear solve performed.')
        return plan
    
    verify=patch_checks(args.order,args.gauss_order)

    if not verify['passed']:
        raise AssertionError('Polynomial patch verification failed')
    
    element_checks=element_verification(args.order)

    if not element_checks['passed']:
        raise AssertionError('One-cell mechanical verification failed')
    
    if args.mode=='verify':
        print(json.dumps({'polynomial':verify,'one_cell_mechanics':element_checks},indent=2))
        print('Small independent-element verification only: no structural solve performed.')
        return {'polynomial':verify,'one_cell_mechanics':element_checks}
    
    if plan['free_and_fixed_translation_DOFs']>8000 and not args.accept_expensive_run:
        parser.error(f"{plan['free_and_fixed_translation_DOFs']} DOFs: large sparse direct solve; "
                     'rerun with --accept-expensive-run if you have adequate RAM/time')
        
    if args.weight_density<=0 or args.end_force>0:
        parser.error('This benchmark requires positive density and nonpositive end force')

    print('Solving',plan['family'],'mesh',plan['physical_element_divisions'],
          'nodes',plan['global_nodes'],'DOFs',plan['free_and_fixed_translation_DOFs'],flush=True)
    
    result=solve_study(args.order,args.nx,args.ny,args.nz,args.weight_density,
         args.end_force,args.gauss_order,section_order=args.section_order,
         field_output_path=args.field_output)
    
    result['polynomial_patch_tests']=verify
    result['one_cell_mechanics_tests']=element_checks

    if args.baseline_cubic_json:
        result['comparison_vs_cubic_Hermite']=baseline_compare(result,args.baseline_cubic_json)

    out=args.out or Path('outputs')/f"hex{(args.order+1)**3}_{args.nx}x{args.ny}x{args.nz}.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    csv_path=out.with_name(out.stem+'_sections.csv')

    with csv_path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(result['source_sections'][0]))
        writer.writeheader();writer.writerows(result['source_sections'])

    print('Tip max abs uy =',result['tip_max_abs_uy_in'],'in; beam reference =',result['reference_tip_uy_in'],'in')
    print('Midspan max sectional shear and moment differences:',
          result['max_midspan_shear_y_difference_pct'],result['max_midspan_moment_x_difference_pct'],'%')
    print('Equilibrium and work/energy:',json.dumps(result['equilibrium'],indent=2))
    print('Elapsed:',result['runtime_seconds'],'seconds')
    print('Wrote',out,'and',csv_path)
    if args.field_output is not None:
        print('Field archive:',args.field_output)

    return result

if __name__=='__main__':
    main()
