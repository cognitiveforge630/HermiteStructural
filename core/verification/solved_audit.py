"""Optional independent beam/stress/energy diagnostics for a solved Hermite beam.

Only computed source-element stresses are labeled FE; elementary beam values
are stored separately. This audit does NOT build an admissible stress field,
check all boundary tractions, or establish a certified error enclosure.
"""
from __future__ import annotations
import numpy as np
from core.reference.mechanics_of_materials import (
    cantilever_fields, rectangular_stresses, euler_bernoulli_bending_energy,
)
from core.verification.stress_tensor import traction, principal_stresses


def _source_at(beam,U,point):
    """Read source FE displacement and stress at one physical point in the beam."""
    p=np.asarray(point,float)
    lo=np.array((-beam.Lx/2,-beam.Ly/2,0.))
    hi=np.array((beam.Lx/2,beam.Ly/2,beam.Lz))
    if p.shape!=(3,) or not np.isfinite(p).all() or np.any(p<lo) or np.any(p>hi):
        raise ValueError('Physical point is outside the beam')
    subdivisions=np.array((beam.nx-1,beam.ny-1,beam.nz-1))
    ratio=(p-lo)/(hi-lo)
    loc=np.minimum((ratio*subdivisions).astype(int),subdivisions-1)
    ix,iy,iz=loc
    e=int((ix*(beam.ny-1)+iy)*(beam.nz-1)+iz)
    coords=beam.get_element_nodes(e)
    extent=np.ptp(coords,axis=0)
    parent=tuple((2*(p-(coords.min(axis=0)))/extent-1).tolist())
    ids,T=beam.element_transform(e)
    prim=T@U[ids]
    B,_,disp=beam.primitive_B(*parent,coords)
    stress=beam.D@(B@prim)
    displacement=disp@prim
    return displacement,stress,e


def source_elastic_strain_energy(beam, U):
    """Independent constitutive Gauss integral, 1/2 integral eps^T D eps dV.

    Uses the element's actual integration rule and FE strain operator. Equality
    with discrete work is an algebraic consistency test, not a 3D accuracy test.
    """
    U=np.asarray(U,float)
    energy=0.
    for e in range(beam.get_num_elements()):
        coords=beam.get_element_nodes(e)
        ids,T=beam.element_transform(e)
        prim=T@U[ids]
        for i,xi in enumerate(beam.gauss_points):
            for j,eta in enumerate(beam.gauss_points):
                for k,zeta in enumerate(beam.gauss_points):
                    B,detJ,_=beam.primitive_B(float(xi),float(eta),float(zeta),coords)
                    eps=B@prim
                    energy+=0.5*float(eps@(beam.D@eps))*detJ*float(
                        beam.gauss_weights[i]*beam.gauss_weights[j]*beam.gauss_weights[k])
    return float(energy)


def solved_mechanics_report(beam, U, force, *, end_force, weight_density,
                            sample_fractions=(.25,.5,.75)):
    """Return signed deflection profiles, stress-cube traces and energy identity.

    Integrates the source FE strain energy, compares to discrete external work,
    and separately reports 1D Euler--Bernoulli bending energy. No arbitrary
    acceptance threshold is imposed on approximate 1D-vs-3D comparisons.
    """
    U=np.asarray(U,float)
    force=np.asarray(force,float)
    expected=(len(beam.grid_points)*beam.ndof_per_node,)
    if U.shape!=expected or force.shape!=expected or not np.isfinite(U).all() or not np.isfinite(force).all():
        raise ValueError('Expected full finite displacement and force vectors')
    if not np.isfinite([end_force,weight_density]).all() or weight_density<0:
        raise ValueError('Invalid applied load data')
    fracs=tuple(float(f) for f in sample_fractions)
    if not fracs or not all(np.isfinite(f) and 0<f<1 for f in fracs):
        raise ValueError('Sample fractions must lie strictly within (0,1)')
    q=-weight_density*beam.Lx*beam.Ly  # signed y-load per inch
    data=[]
    for frac in fracs:
        z=frac*beam.Lz
        exact=cantilever_fields(z,length=beam.Lz,E=beam.orthotropic_constants.Ez,
                                  I=beam.Lx*beam.Ly**3/12.,P=end_force,w=q)
        fe_disp,_,e=_source_at(beam,U,(0.,0.,z))
        samples=[]
        for y in (-beam.Ly/4.,0.,beam.Ly/4.):
            _,s,s_e=_source_at(beam,U,(0.,y,z))
            ref=rectangular_stresses(y,z,width=beam.Lx,depth=beam.Ly,
                                     length=beam.Lz,P=end_force,w=q)
            samples.append({'y_in':float(y),'source_element':s_e,
                'FE_stress_voigt_psi':list(map(float,s)),
                'FE_principal_stress_psi':list(map(float,principal_stresses(s))),
                'FE_traction_on_plus_z_face_psi':list(map(float,traction(s,(0,0,1)))),
                'beam_sigma_zz_psi':float(ref['sigma_zz_psi']),
                'beam_tau_yz_psi':float(ref['tau_yz_psi']),
                'sigma_zz_difference_psi':float(s[2]-ref['sigma_zz_psi']),
                'tau_yz_difference_psi':float(s[4]-ref['tau_yz_psi'])})
        data.append({'fraction_of_length':frac,'z_in':float(z),'source_element_for_displacement':e,
                     'FE_centerline_uy_in':float(fe_disp[1]),
                     'beam_centerline_uy_in':float(exact['uy_in']),
                     'centerline_uy_difference_in':float(fe_disp[1]-exact['uy_in']),
                     'beam_shear_y_lbf':float(exact['shear_y_lbf']),
                     'beam_moment_x_lbf_in':float(exact['moment_x_lbf_in']),
                     'local_stress_samples':samples})
    internal=source_elastic_strain_energy(beam,U)
    external=float(U@force)
    reference=euler_bernoulli_bending_energy(length=beam.Lz,E=beam.orthotropic_constants.Ez,
                                              I=beam.Lx*beam.Ly**3/12.,P=end_force,w=q)
    return {'status':'SOLVED_MECHANICS_DIAGNOSTIC_NOT_CERTIFIED',
            'material_model':'3D orthotropic linear Cauchy elasticity',
            'stress_order':'xx, yy, zz, xy, yz, xz',
            'reference_model':'1D Euler-Bernoulli plus elementary parabolic shear profile',
            'discrete_external_work_lbf_in':external,
            'source_FE_strain_energy_lbf_in':internal,
            'clapeyron_relative_residual':float(abs(external-2*internal)/max(abs(external),abs(2*internal),1e-30)),
            'beam_reference_bending_energy_lbf_in':float(reference),
            'section_samples':data,
            'certified_error_bound':False,
            'caveat':('FE work-energy consistency checks implementation, not exact 3D accuracy; '
                      'beam stresses and displacement are 1D approximations; sample points may '
                      'lie on internal element faces; stress is one-sided there. No full stress '
                      'equilibration or bound has been constructed.')}
