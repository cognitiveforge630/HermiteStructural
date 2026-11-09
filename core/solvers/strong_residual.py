"""Sampled strong-form balance residual of raw cubic FE stress (research).

For 3D Cauchy elasticity the local balance condition is div(sigma)+b=0.
This numerical central-difference residual measures its defect at interior
sample points. It is NOT a guaranteed error estimate or a validated integral;
no external traction or solver/rounding enclosure is included.
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss


def point_residual(beam, U, e, parent=(0.5773502691896257,)*3, *, density=0.0174,
                   finite_difference_step=1e-4):
    """Return physical Cartesian div(sigma)+gravity [lbf/in^3].

    Central finite differences take derivatives of the *actual polynomial
    FE stress* with respect to the parent coordinates, then map via J^-1.
    Only evaluates interior points, avoiding two-sided derivative ambiguity.
    """
    parent=np.asarray(parent,float)
    if parent.shape!=(3,) or not np.all(np.isfinite(parent)) or np.any(np.abs(parent)>=1.):
        raise ValueError('Evaluation must be strictly inside the source element')
    if np.all(parent == 0.):
        raise ValueError('Centroid-based stress diagnostics are prohibited; choose an interior non-centroid point')
    h=float(finite_difference_step)
    if not np.isfinite(h) or h<=0 or np.max(np.abs(parent))+h>=1.:
        raise ValueError('Finite-difference step must keep points inside the element')
    if not np.isfinite(density):
        raise ValueError('Gravity density must be finite')
    xyz=beam.get_element_nodes(e)
    ids,T=beam.element_transform(e)
    primitive=T@np.asarray(U,float)[ids]
    D=beam.elastic_matrix()
    derivatives=np.empty((6,3))
    for axis in range(3):
        direction=np.zeros(3); direction[axis]=h
        Bplus,_,_=beam.primitive_B(*(parent+direction),xyz)
        Bminus,_,_=beam.primitive_B(*(parent-direction),xyz)
        derivatives[:,axis]=D@((Bplus-Bminus)@primitive)/(2*h)
    J=beam.get_hermite_jacobian(*parent,xyz)
    if np.linalg.det(J)<=0:
        raise ValueError('Nonpositive mapping Jacobian')
    # J[parent_coord, physical_coord] ==> grad_global = invJ @ grad_parent.
    gradients=derivatives@np.linalg.inv(J).T
    result=np.array((gradients[0,0]+gradients[3,1]+gradients[5,2],
                     gradients[3,0]+gradients[1,1]+gradients[4,2]-density,
                     gradients[5,0]+gradients[4,1]+gradients[2,2]))
    return result


def sample_volume_residual(beam,U,*,density=0.0174,gauss_order=2,
                           finite_difference_step=1e-4):
    """Discrete source-cell diagnostic; RMS/max over interior sampling points.

    Only even Gauss orders 2 or 4 are permitted in this diagnostic: neither
    includes a source-element centroid. These samples are not a certified norm.
    """
    if type(gauss_order)!=int or gauss_order not in (2, 4):
        raise ValueError('gauss_order must be 2 or 4; center-only diagnostics prohibited')
    u=np.asarray(U,float)
    if u.shape!=(len(beam.grid_points)*6,) or not np.all(np.isfinite(u)):
        raise ValueError('Expected full finite displacement DOFs')
    if not np.isfinite(density):
        raise ValueError('Density must be finite')
    gp,gw=leggauss(gauss_order)
    rms_sum=0.
    total_volume=0.
    maximum=(-1.,None,None)
    for e in range(beam.get_num_elements()):
        coords=beam.get_element_nodes(e)
        for i,xi in enumerate(gp):
            for j,eta in enumerate(gp):
                for k,zeta in enumerate(gp):
                    parent=(float(xi),float(eta),float(zeta))
                    r=point_residual(beam,u,e,parent,density=density,
                        finite_difference_step=finite_difference_step)
                    J=beam.get_hermite_jacobian(*parent,coords)
                    dv=float(np.linalg.det(J))*gw[i]*gw[j]*gw[k]
                    if dv<=0:
                        raise ValueError('Invalid volume Jacobian')
                    rmag=float(np.linalg.norm(r))
                    rms_sum+=dv*rmag**2
                    total_volume+=dv
                    if rmag>maximum[0]:
                        maximum=(rmag,int(e),list(parent))
    return {
        'status':'SAMPLED_STRONG_FORM_RESIDUAL_NOT_CERTIFIED',
        'units':'lbf/in^3',
        'samples_per_source_cell':gauss_order**3,
        'gauss_order':gauss_order,
        'finite_difference_step_parent':float(finite_difference_step),
        'weighted_sample_rms_div_sigma_plus_b_lbf_in3':float(np.sqrt(rms_sum/total_volume)),
        'maximum_sample_div_sigma_plus_b_lbf_in3':maximum[0],
        'maximum_element_index':maximum[1],
        'maximum_parent_coordinates':maximum[2],
        'certified_bound':False,
        'caveats':'Pointwise finite-difference derivative and sampled residual, not an integral enclosure or statically admissible stress. Boundary traction defects are not covered.',
    }
