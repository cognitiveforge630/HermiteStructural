"""Research six-DOF C1-at-vertices, C0-conforming Hermite hexahedron.

The three rotational nodal unknowns ARE physical (one half of curl u at a
vertex). Six symmetric components of the vertex displacement gradient are
recovered *linearly from globally shared nodal translations*, not independently
prescribed. Hence stiffness has a one-ring/nonlocal stencil and **cannot** be
assembled using the original element-only 48-column B matrix.

u_h = sum_a H_a u_a + sum_(a,c,k) R_(a,k) e_c G_(a,c,k),
G_a = sym(recovered_gradient_of_translations_a) + skew(theta_a).

The recovered symmetric gradient is a research choice; it is not an exact
constitutive rotational theory, no proven optimal convergence, and no a
posteriori certification is implied. Straight axis-aligned structured meshes
only. Ordinary Cauchy stress/elasticity is used.
"""
from __future__ import annotations

from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix

from .base import FEMHermiteBeamRegion


class RecoveredRotationHex8(FEMHermiteBeamRegion):
    """Six nodal DOFs: u_x,u_y,u_z, theta_x,theta_y,theta_z."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('gauss_order', 4)
        if kwargs['gauss_order'] < 3:
            raise ValueError('Recovered-rotation element requires >=3 Gauss points')
        super().__init__(*args, **kwargs)
        self.gradient_weights = self._gradient_recovery_weights()
        self.primitive_dofs = 96
        self._element_transform_cache = {}

    def _gradient_recovery_weights(self):
        connectivity = defaultdict(set)
        for element in self.elements:
            for node in element:
                connectivity[int(node)].update(int(n) for n in element if n != node)
        weights = []
        for node, x in enumerate(self.grid_points):
            neighbors = sorted(connectivity[node])
            dx = self.grid_points[neighbors] - x
            if np.linalg.matrix_rank(dx) != 3:
                raise ValueError('Nodal least-squares patch does not span 3D')
            # Grad u_c / dX_k from nodal translation differences.
            recovery = np.linalg.pinv(dx, rcond=1e-12)
            w = {n: recovery[:, j] for j, n in enumerate(neighbors)}
            w[node] = -np.sum(recovery, axis=1)
            weights.append(w)
        return weights

    def element_transform(self, e):
        """Return global support DOFs and T mapping them to 8*(u[3], grad[9])."""
        if e in self._element_transform_cache:
            return self._element_transform_cache[e]
        nodes = [int(n) for n in self.elements[e]]
        support = sorted(set(nodes).union(*(self.gradient_weights[n].keys() for n in nodes)))
        pos = {n: i for i, n in enumerate(support)}
        t = np.zeros((96, len(support)*6))
        for a, node in enumerate(nodes):
            t[12*a:12*a+3, 6*pos[node]:6*pos[node]+3] = np.eye(3)
            for c in range(3):
                for k in range(3):
                    row = 12*a + 3 + 3*c + k
                    for neighbor, w in self.gradient_weights[node].items():
                        t[row, 6*pos[neighbor]+c] += 0.5*w[k]
                        t[row, 6*pos[neighbor]+k] += 0.5*w[c]
            # skew(theta)_{i,j} = -epsilon_{i,j,k} theta_k
            for c,k,theta,sign in ((0,1,2,-1),(0,2,1,1),
                                    (1,0,2,1),(1,2,0,-1),
                                    (2,0,1,-1),(2,1,0,1)):
                t[12*a+3+3*c+k, 6*pos[node]+3+theta] += sign
        dofs = self.global_dofs[support].ravel()
        self._element_transform_cache[e] = (dofs, t)
        return dofs, t

    def primitive_matrices(self, xi, eta, zeta, coords):
        """Displacement and parent-derivative matrices with 96 primitive DOFs."""
        nh, rh, dnh, drh = self.get_hermite_shape_functions_and_derivatives(xi, eta, zeta)
        length_scales = 0.5*np.ptp(coords, axis=0)
        rh = rh * length_scales[:,None]
        drh = drh * length_scales[None,:,None]
        d = np.zeros((3, 96))
        dp = [np.zeros((3, 96)) for _ in range(3)]
        for a in range(8):
            for c in range(3):
                d[c,12*a+c] = nh[a]
                for q in range(3):
                    dp[q][c,12*a+c] = dnh[q,a]
                for k in range(3):
                    col = 12*a+3+3*c+k
                    d[c,col] = rh[k,a]
                    for q in range(3):
                        dp[q][c,col] = drh[q,k,a]
        return d,dp

    def primitive_B(self, xi, eta, zeta, coords):
        disp, derivatives = self.primitive_matrices(xi, eta, zeta, coords)
        J = self.get_hermite_jacobian(xi, eta, zeta, coords)
        detJ = float(np.linalg.det(J))
        if not np.isfinite(detJ) or detJ <= 0:
            raise ValueError('Invalid element Jacobian')
        invJ = np.linalg.inv(J)
        gradients = [invJ @ np.stack([d[c] for d in derivatives]) for c in range(3)]
        B = np.zeros((6, disp.shape[1]))
        B[0] = gradients[0][0]
        B[1] = gradients[1][1]
        B[2] = gradients[2][2]
        B[3] = gradients[0][1]+gradients[1][0]
        B[4] = gradients[1][2]+gradients[2][1]
        B[5] = gradients[0][2]+gradients[2][0]
        return B,detJ,disp

    def strain_at(self, e, U, parent):
        coords = self.get_element_nodes(e)
        B,_,_ = self.primitive_B(*parent, coords)
        dofs,T = self.element_transform(e)
        return B @ (T @ np.asarray(U)[dofs])

    def displacement_at(self, e, U, parent):
        coords=self.get_element_nodes(e)
        disp,_=self.primitive_matrices(*parent,coords)
        dofs,T=self.element_transform(e)
        return disp @ (T@np.asarray(U)[dofs])

    def displacement_operator(self, e, parent):
        """Return (support DOFs, 3 x n_support matrix), for consistent loads."""
        coords = self.get_element_nodes(e)
        disp,_ = self.primitive_matrices(*parent,coords)
        dofs,T = self.element_transform(e)
        return dofs, disp @ T

    def element_primitive_K(self, coords):
        K=np.zeros((self.primitive_dofs,self.primitive_dofs))
        for i,xi in enumerate(self.gauss_points):
            for j,eta in enumerate(self.gauss_points):
                for k,zeta in enumerate(self.gauss_points):
                    B,detJ,_ = self.primitive_B(xi,eta,zeta,coords)
                    K+= B.T @ self.D @ B * detJ * self.gauss_weights[i]*self.gauss_weights[j]*self.gauss_weights[k]
        return K

    def build_global_K(self):
        rows=[];cols=[];values=[]
        cache={}
        full_entries=assembled_entries=0
        maximum_support=maximum_active=0
        for e in range(self.get_num_elements()):
            coords=self.get_element_nodes(e)
            extents=tuple(np.round(np.ptp(coords,axis=0),12))
            if extents not in cache:
                cache[extents]=self.element_primitive_K(coords)
            dofs,T=self.element_transform(e)
            # Neighbor rotations have zero columns. Do not store their exact
            # zero contributions in COO. Retain the original dense-product
            # layout: compacting columns proved slower with the default BLAS
            # thread count on small cases. No numerical threshold is applied.
            full_entries+=len(dofs)**2
            maximum_support=max(maximum_support,len(dofs))
            active=np.flatnonzero(np.any(T!=0.,axis=0))
            maximum_active=max(maximum_active,len(active))
            Ke=T.T @ cache[extents] @ T
            rr,cc=np.nonzero(Ke)
            rows.append(dofs[rr]);cols.append(dofs[cc]);values.append(Ke[rr,cc])
            assembled_entries+=len(rr)
        ndofs = self.grid_points.shape[0]*6
        self.last_assembly_statistics={'unpruned_dense_coo_entries':full_entries,
            'contributed_nonzero_coo_entries':assembled_entries,
            'maximum_full_support_columns':maximum_support,'maximum_active_columns':maximum_active}
        return coo_matrix((np.concatenate(values),
                           (np.concatenate(rows),np.concatenate(cols))),
                          shape=(ndofs,ndofs)).tocsr()

    def build_B_matrix(self, xi, eta, zeta, coords):
        raise NotImplementedError(
            'This nonlocal element requires strain_at(e, U, parent); '
            'the old local 48-column B matrix cannot represent its stiffness')

    def build_force_vector(self):
        raise NotImplementedError('Apply loads with the global-support displacement matrix')
