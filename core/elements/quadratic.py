"""Research quadratic-complete nonlocal Hermite with true nodal spin DOFs.

A two-ring quadratic least-squares recovery supplies symmetric vertex
translation slopes and the mixed second derivatives. Three physical rotations
supply the skew part of the vertex gradient. This is C0 conforming on the
structured beam grids; it uses a nonlocal stiffness stencil and is not a
standard local 48x48 element, nor is it a production-certified element.

For undersampled one- and two-layer coordinate directions where a complete
quadratic fit is impossible, the recovery falls back to linear fitting;
quadratic completeness is therefore *conditional* on rank-9 recovery patches.
"""
from __future__ import annotations
from collections import defaultdict
import numpy as np
from .linear import RecoveredRotationHex8

_PAIRS = ((0,1),(1,2),(0,2))


class QuadraticRotationHex8(RecoveredRotationHex8):
    def _gradient_recovery_weights(self):
        adjacent = defaultdict(set)
        for element in self.elements:
            for n in element:
                adjacent[int(n)].update(int(v) for v in element if v != n)
        first = []
        self.hessian_weights = []
        self.quadratic_rank = []
        scale = np.array((self.Lx/(self.nx-1),self.Ly/(self.ny-1),
                          self.Lz/(self.nz-1)),float)
        for n, point in enumerate(self.grid_points):
            patch = set(adjacent[n])
            for nb in list(patch):
                patch.update(adjacent[nb])
            patch.discard(n)
            neighbors = sorted(patch)
            dx = (self.grid_points[neighbors]-point)/scale
            A = np.column_stack((dx, .5*dx[:,0]**2,.5*dx[:,1]**2,
                .5*dx[:,2]**2,dx[:,0]*dx[:,1],dx[:,1]*dx[:,2],dx[:,0]*dx[:,2]))
            rank=int(np.linalg.matrix_rank(A,tol=1e-11))
            self.quadratic_rank.append(rank)
            if rank == 9:
                pinv=np.linalg.pinv(A,rcond=1e-12)
                grad=pinv[:3]/scale[:,None]
                hessian=np.stack((pinv[6]/(scale[0]*scale[1]),
                                  pinv[7]/(scale[1]*scale[2]),
                                  pinv[8]/(scale[0]*scale[2])))
            else:
                if np.linalg.matrix_rank(dx)<3:
                    raise ValueError('Recovery neighborhood does not span 3D')
                grad=np.linalg.pinv(dx,rcond=1e-12)/scale[:,None]
                hessian=np.zeros((3,len(neighbors)))
            grad_w={v:grad[:,j] for j,v in enumerate(neighbors)}
            hess_w={v:hessian[:,j] for j,v in enumerate(neighbors)}
            grad_w[n]=-np.sum(grad,axis=1)
            hess_w[n]=-np.sum(hessian,axis=1)
            first.append(grad_w);self.hessian_weights.append(hess_w)
        return first

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.primitive_dofs=168
        if self.gauss_order<4:
            raise ValueError('Quadratic rotation research requires >=4 integration points')

    def element_transform(self,e):
        if e in self._element_transform_cache:
            return self._element_transform_cache[e]
        nodes=[int(n) for n in self.elements[e]]
        support=sorted(set(nodes).union(*(self.gradient_weights[n].keys() for n in nodes)))
        pos={n:i for i,n in enumerate(support)}
        t=np.zeros((self.primitive_dofs,len(support)*6))
        for a,node in enumerate(nodes):
            base=21*a
            t[base:base+3,6*pos[node]:6*pos[node]+3]=np.eye(3)
            for c in range(3):
                for k in range(3):
                    row=base+3+3*c+k
                    for neighbor,w in self.gradient_weights[node].items():
                        t[row,6*pos[neighbor]+c]+=0.5*w[k]
                        t[row,6*pos[neighbor]+k]+=0.5*w[c]
                for j in range(3):
                    row=base+12+3*c+j
                    for neighbor,w in self.hessian_weights[node].items():
                        t[row,6*pos[neighbor]+c]+=w[j]
            for c,k,theta,sgn in ((0,1,2,-1),(0,2,1,1),
                                   (1,0,2,1),(1,2,0,-1),
                                   (2,0,1,-1),(2,1,0,1)):
                t[base+3+3*c+k,6*pos[node]+3+theta]+=sgn
        dofs=self.global_dofs[support].ravel()
        self._element_transform_cache[e]=(dofs,t)
        return dofs,t

    def primitive_matrices(self,xi,eta,zeta,coords):
        # Hermite product with value (V) and first-derivative (R) bases.
        node_signs = ((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
                      (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))
        p=(xi,eta,zeta)
        scales=.5*np.ptp(coords,axis=0)
        disp=np.zeros((3,self.primitive_dofs))
        deriv=[np.zeros((3,self.primitive_dofs)) for _ in range(3)]
        # value, three first slopes, and xy/yz/xz mixed slopes. Evaluating
        # these tiny products in bulk avoids hundreds of np.prod reductions
        # per quadrature point without changing any interpolation polynomial.
        use_r=np.array(((0,0,0),(1,0,0),(0,1,0),(0,0,1),
                        (1,1,0),(0,1,1),(1,0,1)),dtype=bool)
        for a,sgns in enumerate(node_signs):
            v=[];r=[];dv=[];dr=[]
            for k in range(3):
                vk,rk,dvk,drk=self._node_axis_shape(p[k],sgns[k])
                v.append(vk);r.append(rk*scales[k])
                dv.append(dvk);dr.append(drk*scales[k])
            factors=np.where(use_r,r,v)
            slopes=np.where(use_r,dr,dv)
            values=factors[:,0]*factors[:,1]*factors[:,2]
            gradients=(slopes[:,0]*factors[:,1]*factors[:,2],
                       factors[:,0]*slopes[:,1]*factors[:,2],
                       factors[:,0]*factors[:,1]*slopes[:,2])
            for c in range(3):
                cols=np.array([21*a+c,21*a+3+3*c,21*a+4+3*c,21*a+5+3*c,
                               21*a+12+3*c,21*a+13+3*c,21*a+14+3*c])
                disp[c,cols]=values
                for q in range(3):
                    deriv[q][c,cols]=gradients[q]
        return disp,deriv
