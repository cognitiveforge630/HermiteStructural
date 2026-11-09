"""Independent elementary cantilever references, in / lbf / psi units.

Coordinates: clamped end z=0, loaded free end z=L, section y measured
from the centroid.  Transverse forces P and w are signed (+y positive).
These Euler--Bernoulli and elementary rectangular-section formulas are
*approximate 1D beam theory*, not exact 3D elasticity or certified FE bounds.

Reference: Beer, Johnston, DeWolf & Mazurek, Mechanics of Materials, 6th ed.,
chapters 4--6, 9 and 11; appendix D (standard beam load cases).
"""
from __future__ import annotations

import numpy as np


def cantilever_fields(z, *, length, E, I, P=0., w=0.):
    """Return signed displacement, slope, shear and bending moment at z.

    ``w`` is the signed transverse force per unit length, *not* a density.
    Moment sign is chosen so sigma_zz = M*y/I for the downward-loaded beam.
    Returns NumPy scalars/arrays matching the input z dimensions.
    """
    values = np.asarray([length, E, I, P, w], dtype=float)
    if not np.isfinite(values).all() or min(values[:3]) <= 0:
        raise ValueError('length, E and I must be positive finite values; loads finite')
    z = np.asarray(z, dtype=float)
    if not np.isfinite(z).all() or np.any(z < 0) or np.any(z > length):
        raise ValueError('z must be finite and within [0, length]')
    rem = length - z
    delta = (P*z**2*(3*length-z)/(6*E*I)
             + w*z**2*(6*length**2-4*length*z+z**2)/(24*E*I))
    slope = (P*z*(2*length-z)/(2*E*I)
             + w*z*(3*length**2-3*length*z+z**2)/(6*E*I))
    return {'uy_in':delta, 'slope_y_per_z':slope,
            'shear_y_lbf':P+w*rem,
            'moment_x_lbf_in':-P*rem-0.5*w*rem**2}


def rectangular_stresses(y, z, *, width, depth, length, P=0., w=0.):
    """Classical bending sigma_zz and parabolic transverse tau_yz (psi).

    Interior Saint-Venant reference only: ignores localized clamp/end effects,
    Poisson and 3D warping stresses. Section origin y=0 is centroid.
    """
    width,depth = float(width),float(depth)
    if not np.isfinite([width,depth]).all() or min(width,depth)<=0:
        raise ValueError('Cross-section dimensions must be positive finite')
    y = np.asarray(y,dtype=float)
    if not np.isfinite(y).all() or np.any(abs(y)>depth/2):
        raise ValueError('y must lie within the rectangular cross section')
    fields=cantilever_fields(z,length=length,E=1.,I=1.,P=P,w=w)
    I=width*depth**3/12.
    return {'sigma_zz_psi':fields['moment_x_lbf_in']*y/I,
            'tau_yz_psi':1.5*fields['shear_y_lbf']/(width*depth)*(1.-(2*y/depth)**2)}


def euler_bernoulli_bending_energy(*, length, E, I, P=0., w=0.):
    """Exact 1D bending energy 0.5*integral(M(z)^2 / EI dz), lbf*in."""
    vals=np.asarray([length,E,I,P,w],float)
    if not np.isfinite(vals).all() or min(vals[:3])<=0:
        raise ValueError('Invalid beam properties or loads')
    return (P**2*length**3/3 + P*w*length**4/4 + w**2*length**5/20)/(2*E*I)
