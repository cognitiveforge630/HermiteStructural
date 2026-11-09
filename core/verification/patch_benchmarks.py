"""Fast source-element numerical patch checks against independent exact fields.

Tests are prescribed kinematic patches, NOT solved load cases, 3D stress
accuracy validations, or mathematical proofs of an a posteriori error bound.
No elasticity calibration or mesh changes are performed.
"""
from __future__ import annotations
import numpy as np


def _prescribed_dofs(beam, gradient, spin):
    """Construct nodal values and rotational DOFs for a supplied deformation."""
    xyz=beam.grid_points
    U=np.zeros((len(xyz),6))
    U[:,:3]=gradient(xyz)
    U[:,3:]=spin(xyz)
    return U.ravel()


def test_modes(beam, *, sample_points=((-0.52,0.37,-0.23),(.27,-.41,.63))):
    """Return mode-by-mode physical stress errors using actual cubic FE B.

    A uniaxial stress state requires free lateral Poisson contractions,
    *not* imposing epsilon_xx=epsilon_yy=0. Pure bending is the polynomial
    Saint-Venant displacement field with lateral curvature terms.
    """
    if beam.orthotropic_constants is None:
        raise ValueError('This benchmark requires prescribed orthotropic constants')
    m=beam.orthotropic_constants
    nzx=m.nuxz*m.Ez/m.Ex
    nzy=m.nuyz*m.Ez/m.Ey
    E=m.Ez
    k=1.7e-6  # in^-1, bending curvature
    elongation=2.3e-5
    gamma=3.5e-5
    # (name, translations at xyz, nodal spin, expected sigma_voigt at xyz)
    cases=(
        ('tension',lambda r: np.stack((-nzx*elongation*r[:,0],-nzy*elongation*r[:,1],elongation*r[:,2]),axis=1),
         lambda r: np.zeros_like(r),lambda r: np.tile([0,0,E*elongation,0,0,0],(len(r),1))),
        ('compression',lambda r: np.stack((nzx*elongation*r[:,0],nzy*elongation*r[:,1],-elongation*r[:,2]),axis=1),
         lambda r: np.zeros_like(r),lambda r: np.tile([0,0,-E*elongation,0,0,0],(len(r),1))),
        ('shear_yz',lambda r: np.stack((np.zeros(len(r)),gamma*r[:,2]/2,gamma*r[:,1]/2),axis=1),
         lambda r: np.zeros_like(r),lambda r: np.tile([0,0,0,0,m.Gyz*gamma,0],(len(r),1))),
        ('pure_bending',lambda r: np.stack((-nzx*k*r[:,0]*r[:,1],
             -0.5*k*r[:,2]**2 +0.5*nzx*k*r[:,0]**2-0.5*nzy*k*r[:,1]**2,
             k*r[:,1]*r[:,2]),axis=1),
         lambda r: np.stack((k*r[:,2],np.zeros(len(r)),nzx*k*r[:,0]),axis=1),
         lambda r: np.stack((np.zeros(len(r)),np.zeros(len(r)),E*k*r[:,1],
                            np.zeros(len(r)),np.zeros(len(r)),np.zeros(len(r))),axis=1)),
    )
    selected=tuple(int(x) for x in (0,beam.get_num_elements()//2,beam.get_num_elements()-1))
    report={}
    for name,u_fun,rot_fun,stress_fun in cases:
        U=_prescribed_dofs(beam,u_fun,rot_fun)
        absolute=0.
        relative=0.
        for e in selected:
            xyz=beam.get_element_nodes(e)
            for p in sample_points:
                phys=beam.hex8_shape_functions(*p)@xyz
                exact=stress_fun(np.asarray(phys)[None,:])[0]
                ids,T=beam.element_transform(e)
                B,_,_=beam.primitive_B(*p,xyz)
                computed=beam.D@(B@(T@U[ids]))
                absolute=max(absolute,float(np.max(np.abs(computed-exact))))
                relative=max(relative,float(np.max(np.abs(computed-exact))/max(np.max(np.abs(exact)),1.0)))
        report[name]={'maximum_absolute_stress_difference_psi':absolute,
                      'maximum_relative_stress_difference':relative,
                      'stress_field_exactly_reproduced_at_test_points':bool(absolute<1e-7)}
    return {'status':'PRESCRIBED_SOURCE_PATCH_DIAGNOSTIC_NOT_CERTIFIED',
            'formulation':type(beam).__name__, 'sampled_source_elements':list(selected),
            'sampled_parent_points':[list(x) for x in sample_points],
            'modes':report, 'certified_error_bound':False,
            'caveat':'Prescribed polynomial displacement modes, not a global solved boundary-value or a posteriori verification'}
