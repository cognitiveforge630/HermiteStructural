"""Experimental conforming affine-complete six-DOF polynomial Hex8.

Important: extra nodal columns are INTERIOR ENRICHMENT AMPLITUDES (in radians
for dimensional scaling), NOT physical nodal rotations. This is therefore a
research comparator, NOT a repaired physical 6-DOF Hermite rotation element.

u_h(X) = sum_a N_a(X) u_a + h_ref * b(xi,eta,zeta) * sum_a N_a(X) a_a
b = (1-xi**2)(1-eta**2)(1-zeta**2)

The bubble vanishes identically on all six faces, so the field is H1
conforming. Translation interpolation is affine exact; no hourglass modes if
integrated with 4-point Gauss. With an extra bubble column for each displacement
component, the 48-DOF element has only six rigid-body null modes on a cube.
"""
from __future__ import annotations
import numpy as np
from FEMHermiteBeamRegion import FEMHermiteBeamRegion

class AffineEnrichedHex8(FEMHermiteBeamRegion):
    def __init__(self,*args,**kwargs):
        if kwargs.get('gauss_order',4)<4:
            raise ValueError('Enrichment requires >=4 Gauss points along each axis')
        kwargs.setdefault('gauss_order',4)
        super().__init__(*args,**kwargs)

    def get_hermite_displacement_matrices(self,xi,eta,zeta,coords):
        N=self.hex8_shape_functions(xi,eta,zeta)
        dN=self.hex8_shape_derivatives(xi,eta,zeta)
        p=np.array([xi,eta,zeta],dtype=float)
        b=float(np.prod(1-p*p))
        db=np.array([-2*p[a]*np.prod([(1-p[t]*p[t]) for t in range(3) if t!=a]) for a in range(3)])
        # The scale converts dimensionless auxiliary amplitudes to displacement.
        # A geometric mean avoids extreme axial scaling for slender Hex8 cells.
        h_ref=float(np.prod(np.ptp(coords,axis=0))**(1/3))
        disp=np.zeros((3,48)); deriv=[np.zeros((3,48)) for _ in range(3)]
        for a in range(8):
            for c in range(3):
                disp[c,6*a+c]=N[a]
                disp[c,6*a+3+c]=h_ref*b*N[a]
                for q in range(3):
                    deriv[q][c,6*a+c]=dN[q,a]
                    deriv[q][c,6*a+3+c]=h_ref*(db[q]*N[a]+b*dN[q,a])
        return disp,*deriv
