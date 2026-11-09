"""Independent, low-cost mathematical checks for affine nodal Hex27/Hex64.

The verifiers do not perform a global solve or certify an error bound.
The polynomial checks sample three interior points, not an assembled patch.
"""
from __future__ import annotations

import itertools
import numpy as np

from core.elements.nodal_hex import element_type, tensor_shapes
from core.materials.redwood import redwood_orthotropic_material
from core.reference.awc import DEFAULT_WEIGHT_DENSITY
from core.reference.beam_theory import E, P


def _element_kernel(p, xyz_extent, D, density, traction, q):
    element = element_type(p)(xyz_extent, D, quadrature_order=q)
    ke = element.get_stiffness_matrix()
    fg = element.integrate_body_force((0., -density, 0.))
    element.clear_force_vector()
    ft = element.facets['+Z'].integrate_traction((0., traction, 0.))
    return ke, fg, ft


def _patch_checks(p, q=None):
    """Analytic partition, rigid, affine and exact p-order polynomial checks."""
    q=p+1 if q is None else q
    h=np.array([2.,3.,4.])
    origin=np.array([1.2,-1.5,7.])
    offsets=np.array(list(itertools.product(np.linspace(0,1,p+1), repeat=3)))
    xyz=origin+offsets*h
    A=np.array([[.04,.025,-.03],[-.015,.035,.02],[.012,-.011,.07]])
    t=np.array([.4,-.2,.7])
    U=xyz@A.T+t
    errors={'partition_unity':0.,'derivative_partition':0.,
            'affine_displacement':0.,'affine_gradient':0.,
            'tensor_product_polynomial':0.,'tensor_product_polynomial_gradient':0.}
    # p=2: quadratic y^2, p=3: cubic y^3; exact with p-order Lagrange.
    for parent in [(-.63,.47,-.19),(.1,-.31,.61),(.77,.24,-.44)]:
        N,ds=tensor_shapes(p,parent)
        dN=ds/(h/2.)[:,None]
        x=origin+(np.asarray(parent)+1.)*h/2.
        errors['partition_unity']=max(errors['partition_unity'],abs(N.sum()-1))
        errors['derivative_partition']=max(errors['derivative_partition'],np.max(abs(dN.sum(axis=1))))
        errors['affine_displacement']=max(errors['affine_displacement'], np.max(abs(N@U-(A@x+t))))
        errors['affine_gradient']=max(errors['affine_gradient'],np.max(abs(U.T@dN.T-A)))
        f=xyz[:,1]**p
        errors['tensor_product_polynomial']=max(errors['tensor_product_polynomial'],abs(N@f-x[1]**p))
        grad=dN@f
        errors['tensor_product_polynomial_gradient']=max(errors['tensor_product_polynomial_gradient'],
            np.max(abs(grad-np.array([0.,p*x[1]**(p-1),0.]))))
    return {'order':p,'max_absolute_errors':{k:float(v) for k,v in errors.items()},
            'passed':all(v<1e-10 for v in errors.values()),
            'certified_error_bound':False}




def _element_verification(p):
    """Independent one-cell checks of consistent loads, rigid modes and energy."""
    material=redwood_orthotropic_material(E)
    D=material.elastic_matrix()
    hx,hy,hz=2.,4.,12.
    rho=DEFAULT_WEIGHT_DENSITY
    end_force=P
    ke,fg,ft=_element_kernel(p,(hx,hy,hz),D,rho,end_force/(hx*hy),p+1)
    nodes=np.array(list(itertools.product(
        np.linspace(-hx/2,hx/2,p+1),
        np.linspace(-hy/2,hy/2,p+1),
        np.linspace(0,hz,p+1))))
    applied=(fg+ft).reshape(-1,3)
    f=applied.sum(axis=0)
    mx=np.cross(nodes,applied).sum(axis=0)[0]
    # Physical rigid-body rotations produce no infinitesimal strains.
    rigid=[]
    for a in range(3):
        v=np.zeros_like(nodes)
        v[:,a]=1
        rigid.append(v.ravel())
        axis=np.zeros(3)
        axis[a]=1
        rigid.append(np.cross(axis,nodes).ravel())
    normK=float(np.linalg.norm(ke,ord=np.inf))
    rigid_res=max(float(np.linalg.norm(ke@r,ord=np.inf)/
                        max(1.,normK*np.linalg.norm(r,ord=np.inf))) for r in rigid)
    affine=np.column_stack((.01*nodes[:,0],-.02*nodes[:,1],.03*nodes[:,2])).ravel()
    strain_energy=.5*float(affine@ke@affine)
    eps=np.array([.01,-.02,.03,0,0,0])
    exact=.5*hx*hy*hz*float(eps@D@eps)
    return build_one_cell_verification_report(
        order=p, applied_force=f, applied_moment_x=mx, end_force=end_force,
        weight_density=rho, extent=(hx, hy, hz), rigid_residual=rigid_res,
        computed_energy=strain_energy, exact_energy=exact, stiffness=ke)


def calculate_one_cell_residuals(*, applied_force, applied_moment_x,
                                 end_force, weight_density, extent,
                                 rigid_residual, computed_energy,
                                 exact_energy, stiffness):
    """Named physical diagnostics, kept separate from the FE integration."""
    hx, hy, hz = extent
    expected_force_y = -weight_density * hx * hy * hz + end_force
    expected_moment_x = -end_force * hz + weight_density * hx * hy * hz * hz / 2
    return {
        'force_y_error_lbf': float(applied_force[1] - expected_force_y),
        'total_force_xz_abs_lbf': float(max(abs(applied_force[0]), abs(applied_force[2]))),
        'applied_moment_x_error_lbf_in': float(applied_moment_x - expected_moment_x),
        'rigid_translation_rotation_relative_K_residual': rigid_residual,
        'constant_strain_relative_energy_error': abs(computed_energy-exact_energy)/exact_energy,
        'local_stiffness_symmetry_relative_error': float(
            np.max(abs(stiffness-stiffness.T)) / max(np.max(abs(stiffness)), 1.0)),
    }


def one_cell_tolerances():
    """Legacy absolute acceptance tolerances (not certified error bounds)."""
    return {
        'force_y_error_lbf': 1e-9,
        'total_force_xz_abs_lbf': 1e-9,
        'applied_moment_x_error_lbf_in': 1e-8,
        'rigid_translation_rotation_relative_K_residual': 1e-10,
        'constant_strain_relative_energy_error': 1e-10,
        'local_stiffness_symmetry_relative_error': 1e-12,
    }


def build_one_cell_verification_report(*, order, **measurements):
    residuals = calculate_one_cell_residuals(**measurements)
    tolerances = one_cell_tolerances()
    return {
        'order': order,
        'residuals': residuals,
        'tolerances': tolerances,
        'passed': all(abs(residuals[key]) < limit for key, limit in tolerances.items()),
        'certified_error_bound': False,
    }



class PolynomialPatchVerifier:
    """Verify partition, affine fields and polynomial reproduction.

    :param element: Hex27Element or Hex64Element instance.
    """
    def verify(self, element):
        return _patch_checks(element.order, element.quadrature_order)


class ElementMechanicsVerifier:
    """Verify one-cell consistent loads, rigid modes, strain energy, symmetry."""
    def verify(self, element):
        # Use the element's order, so both element families share the same tests.
        return _element_verification(element.order)

