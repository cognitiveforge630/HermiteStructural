"""Optional cubic-recovery Hermite research variant, six DOFs per vertex.

Three-ring shared recovery fits total-degree <=3 Taylor terms and supplies
the missing xyz mixed derivative in tensor Hermite interpolation. Pure cubic
terms in a direction need >=4 distinct vertex layers. Unsupported terms are
explicitly omitted and recorded; full 3D cubic completeness is conditional.
No historical formulation is changed.
"""
from collections import defaultdict
import itertools
import math
import numpy as np
from quadratic_rotation_hex8 import QuadraticRotationHex8


class CubicRotationHex8(QuadraticRotationHex8):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.primitive_dofs=192

    def _gradient_recovery_weights(self):
        adjacency=defaultdict(set)
        for element in self.elements:
            for n in element:
                adjacency[int(n)].update(int(v) for v in element if v!=n)
        scale=np.array([self.Lx/(self.nx-1),self.Ly/(self.ny-1),self.Lz/(self.nz-1)])
        levels=[len(np.unique(self.grid_points[:,k])) for k in range(3)]
        if min(levels)<3 or levels[2]<4:
            raise ValueError('Cubic bending recovery requires >=2x2x3 elements')
        terms=[powers for degree in (1,2,3) for powers in itertools.product(range(4),repeat=3)
               if sum(powers)==degree and all(powers[k]<levels[k] for k in range(3))]
        self.cubic_terms=[list(p) for p in terms]
        self.cubic_rank=[]
        self.cubic_condition=[]
        self.hessian_weights=[]
        self.third_weights=[]
        first=[]
        def derivative(pinv,powers):
            return pinv[terms.index(powers)]/np.prod(scale**np.array(powers))
        for node,point in enumerate(self.grid_points):
            patch={node}
            for _ in range(3):
                patch.update(set().union(*(adjacency[n] for n in list(patch))))
            patch.discard(node)
            neighbors=sorted(patch)
            dx=(self.grid_points[neighbors]-point)/scale
            matrix=np.column_stack([np.prod(dx**np.array(p),axis=1)/
                                    math.prod(math.factorial(v) for v in p) for p in terms])
            rank=int(np.linalg.matrix_rank(matrix,tol=1e-11))
            self.cubic_rank.append(rank)
            self.cubic_condition.append(float(np.linalg.cond(matrix)))
            if rank!=len(terms):
                raise ValueError(f'Cubic recovery patch {node} rank {rank} < {len(terms)}')
            pinv=np.linalg.pinv(matrix,rcond=1e-12)
            grad=np.stack([derivative(pinv,p) for p in ((1,0,0),(0,1,0),(0,0,1))])
            hess=np.stack([derivative(pinv,p) for p in ((1,1,0),(0,1,1),(1,0,1))])
            third=derivative(pinv,(1,1,1))
            gw={n:grad[:,j] for j,n in enumerate(neighbors)}
            hw={n:hess[:,j] for j,n in enumerate(neighbors)}
            tw={n:third[j] for j,n in enumerate(neighbors)}
            gw[node],hw[node],tw[node]=-grad.sum(axis=1),-hess.sum(axis=1),-third.sum()
            first.append(gw)
            self.hessian_weights.append(hw)
            self.third_weights.append(tw)
        # Compatibility rank indicator used by older accuracy reports.
        self.quadratic_rank=[9]*len(self.grid_points)
        return first

    def element_transform(self,e):
        if e in self._element_transform_cache:
            return self._element_transform_cache[e]
        # Parent fills its 21-row/node layout into a sufficiently large matrix.
        dofs,old=super().element_transform(e)
        transform=np.zeros_like(old)
        nodes=[int(n) for n in self.elements[e]]
        support=(dofs.reshape(-1,6)[:,0]//6).tolist()
        position={n:i for i,n in enumerate(support)}
        for a,node in enumerate(nodes):
            transform[24*a:24*a+21]=old[21*a:21*a+21]
            for c in range(3):
                for neighbor,w in self.third_weights[node].items():
                    transform[24*a+21+c,6*position[neighbor]+c]+=w
        self._element_transform_cache[e]=(dofs,transform)
        return dofs,transform

    def primitive_matrices(self,xi,eta,zeta,coords):
        old,dp=super().primitive_matrices(xi,eta,zeta,coords)
        disp=np.zeros((3,192))
        derivatives=[np.zeros((3,192)) for _ in range(3)]
        signs=((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
               (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))
        scales=.5*np.ptp(coords,axis=0)
        for a,sgn in enumerate(signs):
            disp[:,24*a:24*a+21]=old[:,21*a:21*a+21]
            for k in range(3):
                derivatives[k][:,24*a:24*a+21]=dp[k][:,21*a:21*a+21]
            r,dr=[],[]
            for k,p in enumerate((xi,eta,zeta)):
                _,rk,_,drk=self._node_axis_shape(p,sgn[k])
                r.append(rk*scales[k]);dr.append(drk*scales[k])
            value=r[0]*r[1]*r[2]
            slope=(dr[0]*r[1]*r[2],r[0]*dr[1]*r[2],r[0]*r[1]*dr[2])
            for c in range(3):
                disp[c,24*a+21+c]=value
                for k in range(3):
                    derivatives[k][c,24*a+21+c]=slope[k]
        return disp,derivatives
