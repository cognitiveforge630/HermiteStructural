"""Independent tensor-product Q2 Cauchy solid on axis-aligned brick meshes.

nx/ny/nz count coarse grid vertices, as for build_beam; each element adds
mid-edge/mid-face/interior displacement nodes. This is a reference CANDIDATE,
whose cantilever accuracy must still be qualified by refinement.
"""
import itertools
import numpy as np
from scipy.sparse import coo_matrix
from reference_hex8 import CauchyHex8Reference


class CauchyHex27Reference(CauchyHex8Reference):
    local_indices = np.array(list(itertools.product(range(3), repeat=3)))
    faces = tuple((axis, side, tuple(i for i, p in enumerate(
        itertools.product(range(3), repeat=3)) if p[axis] == (0 if side == -1 else 2)))
        for axis in range(3) for side in (-1, 1))

    def _generate_grid(self):
        self.node_counts = tuple(2*(n-1)+1 for n in (self.nx, self.ny, self.nz))
        axes = [np.linspace(-self.Lx/2, self.Lx/2, self.node_counts[0]),
                np.linspace(-self.Ly/2, self.Ly/2, self.node_counts[1]),
                np.linspace(0., self.Lz, self.node_counts[2])]
        self.grid_points = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)

    def _generate_elements(self):
        _, ny, nz = self.node_counts
        elements = []
        for i, j, k in itertools.product(range(self.nx-1), range(self.ny-1), range(self.nz-1)):
            p = self.local_indices + 2*np.array([i, j, k])
            elements.append(p[:, 0]*ny*nz+p[:, 1]*nz+p[:, 2])
        self.elements = np.asarray(elements, dtype=int)

    @staticmethod
    def axis_shapes(s):
        return np.array([.5*s*(s-1), 1-s*s, .5*s*(s+1)]), np.array([s-.5, -2*s, s+.5])

    def shapes(self, xi, eta, zeta):
        values, deriv = zip(*(self.axis_shapes(s) for s in (xi, eta, zeta)))
        N = np.einsum('i,j,k->ijk', *values).ravel()
        dN = np.array([np.einsum('i,j,k->ijk',
            *(deriv[k] if axis == k else values[k] for k in range(3))).ravel()
            for axis in range(3)])
        return N, dN

    def hex8_shape_functions(self, xi, eta, zeta):
        # Compatibility geometry API: all 27 nodal positions are mapped.
        return self.shapes(xi, eta, zeta)[0]

    def hex8_shape_derivatives(self, xi, eta, zeta):
        return self.shapes(xi, eta, zeta)[1]

    def get_hermite_displacement_matrices(self, xi, eta, zeta, coords):
        N, dN = self.shapes(xi, eta, zeta)
        disp = np.zeros((3, 81))
        dp = [np.zeros((3, 81)) for _ in range(3)]
        for c in range(3):
            disp[c, c::3] = N
            for k in range(3):
                dp[k][c, c::3] = dN[k]
        return disp, *dp

    def build_B_matrix(self, xi, eta, zeta, coords):
        N, dN = self.shapes(xi, eta, zeta)
        J = dN @ coords
        det = float(np.linalg.det(J))
        if det <= 0 or not np.isfinite(det):
            raise ValueError('Invalid Q2 Jacobian')
        grad = np.linalg.solve(J, dN)
        B = np.zeros((6, 81))
        for c in range(3):
            B[c, c::3] = grad[c]
        for row, (a, b) in enumerate(((0,1), (1,2), (0,2)), start=3):
            B[row, a::3], B[row, b::3] = grad[b], grad[a]
        disp = np.zeros((3, 81))
        for c in range(3):
            disp[c, c::3] = N
        return B, det, disp, None

    def build_element_K(self, coords):
        K = np.zeros((81, 81))
        for i, j, k in itertools.product(range(self.gauss_order), repeat=3):
            p = self.gauss_points[[i, j, k]]
            B, det, _, _ = self.build_B_matrix(*p, coords)
            K += B.T @ self.D @ B * det*np.prod(self.gauss_weights[[i,j,k]])
        return K

    def build_global_K(self):
        rows, cols, values, cache = [], [], [], {}
        for e in range(self.get_num_elements()):
            xyz = self.get_element_nodes(e)
            key = tuple(np.round(np.ptp(xyz, axis=0), 12))
            if key not in cache:
                cache[key] = self.build_element_K(xyz)
            ids = self.get_element_global_dofs(e).ravel()
            rows.append(np.repeat(ids, len(ids)))
            cols.append(np.tile(ids, len(ids)))
            values.append(cache[key].ravel())
        size = len(self.grid_points)*3
        return coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(size, size)).tocsr()

    def build_force_vector(self):
        raise NotImplementedError('Use explicit consistent loads from the accuracy drivers')
