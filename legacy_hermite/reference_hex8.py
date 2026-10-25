"""Independent affine-complete, 3-DOF Cauchy Hex8 reference element.

This is NOT a repaired six-DOF Hermite element.  It deliberately uses the
standard trilinear displacement interpolation, providing a reliable control
for isolating legacy Hermite interpolation and estimator issues.

Units and orthotropic D are those of FEMHermiteBeamRegion.  The public mesh,
geometry and B-matrix API is intentionally compatible with aposteriori/residual.
"""

import numpy as np
from scipy.sparse import coo_matrix

from FEMHermiteBeamRegion import FEMHermiteBeamRegion


class CauchyHex8Reference(FEMHermiteBeamRegion):
    """Trilinear C0 hexahedral solid, three translations per node."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ndof_per_node = 3
        self._generate_global_dofs()

    def get_hermite_displacement_matrices(self, xi, eta, zeta, coords):
        """Compatibility interface; these are *Lagrange*, not Hermite shapes."""
        N = self.hex8_shape_functions(xi, eta, zeta)
        dN = self.hex8_shape_derivatives(xi, eta, zeta)
        N_disp = np.zeros((3, 24))
        derivatives = [np.zeros((3, 24)) for _ in range(3)]
        for a in range(8):
            cols = slice(3 * a, 3 * a + 3)
            N_disp[:, cols] = np.eye(3) * N[a]
            for axis in range(3):
                derivatives[axis][:, cols] = np.eye(3) * dN[axis, a]
        return (N_disp, *derivatives)

    def build_B_matrix(self, xi, eta, zeta, coords):
        N_disp, *dN_parent = self.get_hermite_displacement_matrices(xi, eta, zeta, coords)
        J = self.get_hermite_jacobian(xi, eta, zeta, coords)
        detJ = float(np.linalg.det(J))
        if not np.isfinite(detJ) or detJ <= 0:
            raise ValueError('Invalid hexahedron Jacobian')
        # J[parent_axis, physical_axis]; gradient_X = inverse(J) @ gradient_parent.
        invJ = np.linalg.inv(J)
        gradients = [invJ @ np.stack([d[row] for d in dN_parent]) for row in range(3)]
        B = np.empty((6, 24))
        B[0] = gradients[0][0]
        B[1] = gradients[1][1]
        B[2] = gradients[2][2]
        B[3] = gradients[0][1] + gradients[1][0]
        B[4] = gradients[1][2] + gradients[2][1]
        B[5] = gradients[0][2] + gradients[2][0]
        return B, detJ, N_disp, None

    def build_element_K(self, coords):
        Ke = np.zeros((24, 24))
        for i, xi in enumerate(self.gauss_points):
            for j, eta in enumerate(self.gauss_points):
                for k, zeta in enumerate(self.gauss_points):
                    B, detJ, _, _ = self.build_B_matrix(xi, eta, zeta, coords)
                    Ke += (B.T @ self.D @ B * detJ *
                           self.gauss_weights[i] * self.gauss_weights[j] * self.gauss_weights[k])
        return Ke

    def build_global_K(self):
        rows, cols, vals = [], [], []
        stiffness_cache = {}
        for e in range(self.get_num_elements()):
            coords = self.get_element_nodes(e)
            dofs = self.get_element_global_dofs(e).ravel()
            extents = tuple(np.round(np.ptp(coords, axis=0), 12))
            if extents not in stiffness_cache:
                stiffness_cache[extents] = self.build_element_K(coords)
            Ke = stiffness_cache[extents]
            rows.append(np.repeat(dofs, 24))
            cols.append(np.tile(dofs, 24))
            vals.append(Ke.ravel())
        n = self.grid_points.shape[0] * self.ndof_per_node
        return coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()

    def build_force_vector(self):
        """Cauchy-consistent volumetric load and the existing z=L traction.

        Keeps the reference beam's legacy gamma (z body-force density) and w
        (z end-face traction) meaning without invoking six-DOF Hermite loads.
        """
        F = np.zeros(self.grid_points.shape[0] * self.ndof_per_node)
        body = np.array((0., 0., self.gamma), dtype=float)
        traction = np.array((0., 0., self.w), dtype=float)
        for e in range(self.get_num_elements()):
            coords = self.get_element_nodes(e)
            dofs = self.get_element_global_dofs(e).ravel()
            Fe = np.zeros(24)
            for i, xi in enumerate(self.gauss_points):
                for j, eta in enumerate(self.gauss_points):
                    for k, zeta in enumerate(self.gauss_points):
                        N_disp, *_ = self.get_hermite_displacement_matrices(xi, eta, zeta, coords)
                        detJ = float(np.linalg.det(self.get_hermite_jacobian(xi, eta, zeta, coords)))
                        Fe += (N_disp.T @ body * detJ *
                               self.gauss_weights[i] * self.gauss_weights[j] * self.gauss_weights[k])
            if e in self._elems_max_z:
                for i, xi in enumerate(self.gauss_points):
                    for j, eta in enumerate(self.gauss_points):
                        N_disp, *_ = self.get_hermite_displacement_matrices(xi, eta, 1., coords)
                        J = self.get_hermite_jacobian(xi, eta, 1., coords)
                        area = float(np.linalg.norm(np.cross(J[0],J[1])))
                        Fe += N_disp.T @ traction * area * self.gauss_weights[i] * self.gauss_weights[j]
            F[dofs] += Fe
        return F
