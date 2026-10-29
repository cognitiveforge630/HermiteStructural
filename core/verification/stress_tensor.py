"""Conventional symmetric Cauchy stress-cube operations.

Voigt ordering matches core.elements.base exactly: xx, yy, zz, xy, yz, xz.
The strain vector uses ENGINEERING shear gamma=2*epsilon_ij; the stress
vector uses physical shear tau_ij (no factor of two). No Cosserat couple
stress or skew stress is represented here.
"""
from __future__ import annotations
import numpy as np


def _vector(values):
    v=np.asarray(values,dtype=float)
    if v.shape!=(6,) or not np.isfinite(v).all():
        raise ValueError('Expected six finite components in xx,yy,zz,xy,yz,xz order')
    return v


def stress_matrix(s):
    xx,yy,zz,xy,yz,xz=_vector(s)
    return np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])


def strain_matrix(e):
    xx,yy,zz,xy,yz,xz=_vector(e)
    return np.array([[xx,xy/2,xz/2],[xy/2,yy,yz/2],[xz/2,yz/2,zz]])


def stress_components(t):
    t=np.asarray(t,dtype=float)
    if t.shape!=(3,3) or not np.isfinite(t).all() or not np.allclose(t,t.T,rtol=0,atol=1e-12):
        raise ValueError('Expected a finite symmetric 3x3 Cauchy stress tensor')
    return np.array([t[0,0],t[1,1],t[2,2],t[0,1],t[1,2],t[0,2]])


def traction(stress, normal):
    """Cauchy traction sigma @ n for a *unit* oriented physical normal."""
    n=np.asarray(normal,dtype=float)
    if n.shape!=(3,) or not np.isfinite(n).all() or not np.isclose(np.linalg.norm(n),1.,rtol=0,atol=1e-12):
        raise ValueError('normal must be a finite unit vector')
    return stress_matrix(stress)@n


def rotate_stress(stress, R):
    """Active rotation R of tensor and physical coordinate frame."""
    R=np.asarray(R,dtype=float)
    if R.shape!=(3,3) or not np.isfinite(R).all() or not np.allclose(R.T@R,np.eye(3),atol=1e-12) or not np.isclose(np.linalg.det(R),1.,atol=1e-12):
        raise ValueError('R must be a proper orthogonal rotation')
    return stress_components(R@stress_matrix(stress)@R.T)


def principal_stresses(stress):
    """Ordered sigma_1 >= sigma_2 >= sigma_3, all in input stress units."""
    return np.linalg.eigvalsh(stress_matrix(stress))[::-1]


def elastic_energy_density(stress, engineering_strain):
    """1/2 sigma:epsilon = 1/2 stress_voigt dot engineering_strain.

    Units lbf/in^2 (energy per volume) when stress in psi.
    """
    return 0.5*float(_vector(stress)@_vector(engineering_strain))


def cube_force_and_moment(stress, edge):
    """Resultant force/moment on a uniform-stress infinitesimal cube.

    A conceptual stress-element identity for constant symmetric stress, not a
    verification of spatial equilibrium of a computed nonuniform FE field.
    Force units = stress*edge^2; moment units = stress*edge^3.
    """
    if not np.isfinite(edge) or edge<=0:
        raise ValueError('Cube edge must be positive finite')
    force=np.zeros(3)
    moment=np.zeros(3)
    for axis in range(3):
        for sign in (-1.,1.):
            normal=np.zeros(3);normal[axis]=sign
            f=traction(stress,normal)*edge**2
            x=normal*edge/2
            force+=f
            moment+=np.cross(x,f)
    return force,moment
