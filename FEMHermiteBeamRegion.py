import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interpn
from numpy.polynomial.legendre import leggauss

DOF_UX = 0
DOF_UY = 1
DOF_UZ = 2
DOF_MX = 3
DOF_MY = 4
DOF_MZ = 5

class FEMHermiteBeamRegion:
    """8-node hexahedral solid with Hermitian translation/rotation shapes.

    Each node has translations ``u, v, w`` and rotations ``theta_u, theta_v,
    theta_w``. The active stiffness follows the Hermitian Hexa8 formulation:
    translation DOFs use ``NH`` shape functions and rotation DOFs use ``RH``
    shape functions in a standard 6-component small-strain matrix.
    """
    def __init__(self, Lx, Ly, Lz, nx, ny, nz, E, nu, gamma, w, tol=1e-8, l_c=0.001, H_scale=1.0):
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.ndof_per_node = 6
        self.E = E
        self.nu = nu
        self.gamma = gamma
        self.w = w
        self.tol = tol
        self.l_c = l_c
        self.lambda_ = self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
        self.mu = self.E / (2 * (1 + self.nu))
        self.mu_c = 10 * self.mu
        self.gamma_c = 2 * self.mu * self.l_c ** 2
        self.H_scale = H_scale
        self._generate_grid()
        self._generate_elements()
        self._generate_global_dofs()
        self._elems_min_z = self.get_elements_at_min_z()
        self._elems_max_z = self.get_elements_at_max_z()
        self._nodes_min_z = self.get_nodes_at_min_z()
        self._nodes_max_z = self.get_nodes_at_max_z()
        self.D = self.elastic_matrix()
        # The Hermitian Hexa8 bending formulation locks under full 3-point
        # integration on slender beams. Two-point integration is stable here
        # and matches the cantilever benchmark without the 1-point hourglass.
        self.gauss_points, self.gauss_weights = leggauss(2)

    def _generate_grid(self):
        x = np.linspace(-self.Lx/2, self.Lx/2, self.nx)
        y = np.linspace(-self.Ly/2, self.Ly/2, self.ny)
        z = np.linspace(0, self.Lz, self.nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        self.grid_points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))

    def _generate_elements(self):
        elems = []
        for i in range(self.nx - 1):
            for j in range(self.ny - 1):
                for k in range(self.nz - 1):
                    n0 = i*(self.ny*self.nz) + j*self.nz + k
                    n1 = (i+1)*(self.ny*self.nz) + j*self.nz + k
                    n2 = (i+1)*(self.ny*self.nz) + (j+1)*self.nz + k
                    n3 = i*(self.ny*self.nz) + (j+1)*self.nz + k
                    n4 = i*(self.ny*self.nz) + j*self.nz + (k+1)
                    n5 = (i+1)*(self.ny*self.nz) + j*self.nz + (k+1)
                    n6 = (i+1)*(self.ny*self.nz) + (j+1)*self.nz + (k+1)
                    n7 = i*(self.ny*self.nz) + (j+1)*self.nz + (k+1)
                    elems.append([n0, n1, n2, n3, n4, n5, n6, n7])
        self.elements = np.array(elems, dtype=int)

    def _generate_global_dofs(self):
        n_nodes = self.grid_points.shape[0]
        all_dofs = np.arange(n_nodes * self.ndof_per_node)
        self.global_dofs = all_dofs.reshape(n_nodes, self.ndof_per_node)

    def get_num_elements(self):
        return self.elements.shape[0]

    def get_element_nodes(self, elem_idx):
        idx = self.elements[elem_idx]
        return self.grid_points[idx]

    def get_element_global_dofs(self, elem_idx):
        node_indices = self.elements[elem_idx]
        return self.global_dofs[node_indices]

    def get_elements_at_min_z(self):
        min_z = 0
        elems = [eid for eid in range(self.get_num_elements()) if np.any(np.isclose(self.get_element_nodes(eid)[:,2], min_z, atol=self.tol))]
        return elems

    def get_elements_at_max_z(self):
        max_z = self.Lz
        elems = [eid for eid in range(self.get_num_elements()) if np.any(np.isclose(self.get_element_nodes(eid)[:,2], max_z, atol=self.tol))]
        return elems

    def get_nodes_at_min_z(self):
        min_z = 0
        nodes = [n for n in range(self.grid_points.shape[0]) if np.isclose(self.grid_points[n,2], min_z, atol=self.tol)]
        return nodes

    def get_nodes_at_max_z(self):
        max_z = self.Lz
        nodes = [n for n in range(self.grid_points.shape[0]) if np.isclose(self.grid_points[n,2], max_z, atol=self.tol)]
        return nodes

    def hex8_shape_functions(self, xi, eta, zeta):
        N = 0.125 * np.array([
            (1 - xi) * (1 - eta) * (1 - zeta),
            (1 + xi) * (1 - eta) * (1 - zeta),
            (1 + xi) * (1 + eta) * (1 - zeta),
            (1 - xi) * (1 + eta) * (1 - zeta),
            (1 - xi) * (1 - eta) * (1 + zeta),
            (1 + xi) * (1 - eta) * (1 + zeta),
            (1 + xi) * (1 + eta) * (1 + zeta),
            (1 - xi) * (1 + eta) * (1 + zeta)
        ])
        return N

    def hex8_shape_derivatives(self, xi, eta, zeta):
        dN_dn = 0.125 * np.array([
            [ -(1 - eta) * (1 - zeta), -(1 - xi) * (1 - zeta), -(1 - xi) * (1 - eta) ],
            [ (1 - eta) * (1 - zeta), -(1 + xi) * (1 - zeta), -(1 + xi) * (1 - eta) ],
            [ (1 + eta) * (1 - zeta), (1 + xi) * (1 - zeta), -(1 + xi) * (1 + eta) ],
            [ -(1 + eta) * (1 - zeta), (1 - xi) * (1 - zeta), -(1 - xi) * (1 + eta) ],
            [ -(1 - eta) * (1 + zeta), -(1 - xi) * (1 + zeta), (1 - xi) * (1 - eta) ],
            [ (1 - eta) * (1 + zeta), -(1 + xi) * (1 + zeta), (1 + xi) * (1 - eta) ],
            [ (1 + eta) * (1 + zeta), (1 + xi) * (1 + zeta), (1 + xi) * (1 + eta) ],
            [ -(1 + eta) * (1 + zeta), (1 - xi) * (1 + zeta), (1 - xi) * (1 + eta) ]
        ]).T
        return dN_dn

    def hermite_H_functions(self, s):
        h1 = 0.25 * (2 - 3*s + s**3)
        h2 = 0.25 * (1 - s - s**2 + s**3)
        h3 = 0.25 * (2 + 3*s - s**3)
        h4 = 0.25 * (-1 - s + s**2 + s**3)
        return self.H_scale * np.array([h1, h2, h3, h4])

    def dhermite_H_functions(self, s):
        dh1 = 0.25 * (-3 + 3*s**2)
        dh2 = 0.25 * (-1 -2*s + 3*s**2)
        dh3 = 0.25 * (3 - 3*s**2)
        dh4 = 0.25 * (-1 + 2*s + 3*s**2)
        return self.H_scale * np.array([dh1, dh2, dh3, dh4])

    def get_lagrange_shapes_and_derivs(self, xi, eta, zeta):
        N_u = self.hex8_shape_functions(xi, eta, zeta)
        N_phi = N_u
        dN_dn = self.hex8_shape_derivatives(xi, eta, zeta)
        dNu_dxi = dN_dn[0, :]
        dNu_deta = dN_dn[1, :]
        dNu_dzeta = dN_dn[2, :]
        dNp_dxi = dN_dn[0, :]
        dNp_deta = dN_dn[1, :]
        dNp_dzeta = dN_dn[2, :]
        return N_u, N_phi, dNu_dxi, dNu_deta, dNu_dzeta, dNp_dxi, dNp_deta, dNp_dzeta

    def _node_axis_shape(self, s, sign):
        H = self.hermite_H_functions(s)
        dH = self.dhermite_H_functions(s)
        if sign < 0:
            return H[0], H[1], dH[0], dH[1]
        return H[2], H[3], dH[2], dH[3]

    def get_hermite_shape_functions_and_derivatives(self, xi, eta, zeta):
        node_signs = np.array([
            [-1, -1, -1],
            [1, -1, -1],
            [1, 1, -1],
            [-1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [1, 1, 1],
            [-1, 1, 1],
        ])

        NH = np.zeros(8)
        RH = np.zeros((3, 8))
        dNH = np.zeros((3, 8))
        dRH = np.zeros((3, 3, 8))

        for n, (sx, sy, sz) in enumerate(node_signs):
            vx, rx, dvx, drx = self._node_axis_shape(xi, sx)
            vy, ry, dvy, dry = self._node_axis_shape(eta, sy)
            vz, rz, dvz, drz = self._node_axis_shape(zeta, sz)

            NH[n] = vx * vy * vz
            dNH[:, n] = [
                dvx * vy * vz,
                vx * dvy * vz,
                vx * vy * dvz,
            ]

            RH[:, n] = [
                rx * vy * vz,
                vx * ry * vz,
                vx * vy * rz,
            ]
            dRH[:, 0, n] = [drx * vy * vz, rx * dvy * vz, rx * vy * dvz]
            dRH[:, 1, n] = [dvx * ry * vz, vx * dry * vz, vx * ry * dvz]
            dRH[:, 2, n] = [dvx * vy * rz, vx * dvy * rz, vx * vy * drz]

        return NH, RH, dNH, dRH

    def get_hermite_jacobian(self, xi, eta, zeta, coords):
        return self.hex8_shape_derivatives(xi, eta, zeta) @ coords

    def get_hermite_displacement_matrices(self, xi, eta, zeta, coords):
        NH, RH, dNH, dRH = self.get_hermite_shape_functions_and_derivatives(xi, eta, zeta)
        length_scales = np.array([
            0.5 * (coords[:, 0].max() - coords[:, 0].min()),
            0.5 * (coords[:, 1].max() - coords[:, 1].min()),
            0.5 * (coords[:, 2].max() - coords[:, 2].min()),
        ])
        RH = RH * length_scales[:, None]
        dRH = dRH * length_scales[None, :, None]
        N_disp = np.zeros((3, 48))
        dN_dxi = np.zeros((3, 48))
        dN_deta = np.zeros((3, 48))
        dN_dzeta = np.zeros((3, 48))

        for n in range(8):
            col = 6 * n
            N_disp[0, col + DOF_UX] = NH[n]
            N_disp[1, col + DOF_UY] = NH[n]
            N_disp[2, col + DOF_UZ] = NH[n]
            for row, dof in enumerate((DOF_UX, DOF_UY, DOF_UZ)):
                dN_dxi[row, col + dof] = dNH[0, n]
                dN_deta[row, col + dof] = dNH[1, n]
                dN_dzeta[row, col + dof] = dNH[2, n]

            # Rotations enter displacement as small rigid-rotation slopes:
            # u = theta x r. Each RH axis is scaled by the physical half-length
            # so nodal rotations remain dimensionless.
            rotation_couplings = (
                (0, DOF_MY, 2, 1.0),
                (0, DOF_MZ, 1, -1.0),
                (1, DOF_MZ, 0, 1.0),
                (1, DOF_MX, 2, -1.0),
                (2, DOF_MX, 1, 1.0),
                (2, DOF_MY, 0, -1.0),
            )
            for row, dof, axis, sign in rotation_couplings:
                N_disp[row, col + dof] += sign * RH[axis, n]
                dN_dxi[row, col + dof] += sign * dRH[0, axis, n]
                dN_deta[row, col + dof] += sign * dRH[1, axis, n]
                dN_dzeta[row, col + dof] += sign * dRH[2, axis, n]

        return N_disp, dN_dxi, dN_deta, dN_dzeta

    def get_hermite_shapes_and_derivs(self, xi, eta, zeta):
        NH, RH, dNH, dRH = self.get_hermite_shape_functions_and_derivatives(xi, eta, zeta)
        return NH, RH, dNH[0], dNH[1], dNH[2], dRH[0], dRH[1], dRH[2]

    def build_B_matrix(self, xi, eta, zeta, coords):
        N_disp, dN_dxi, dN_deta, dN_dzeta = self.get_hermite_displacement_matrices(xi, eta, zeta, coords)
        J = self.get_hermite_jacobian(xi, eta, zeta, coords)
        detJ = np.linalg.det(J)
        J_inv = np.linalg.inv(J)
        gradients = [
            J_inv.T @ np.vstack((dN_dxi[row], dN_deta[row], dN_dzeta[row]))
            for row in range(3)
        ]

        B = np.zeros((6, 48))
        B[0, :] = gradients[0][0, :]
        B[1, :] = gradients[1][1, :]
        B[2, :] = gradients[2][2, :]
        B[3, :] = gradients[0][1, :] + gradients[1][0, :]
        B[4, :] = gradients[1][2, :] + gradients[2][1, :]
        B[5, :] = gradients[0][2, :] + gradients[2][0, :]

        return B, detJ, N_disp, None

    def build_global_K(self):
        num_dofs = self.grid_points.shape[0] * self.ndof_per_node
        K_data = []
        K_row = []
        K_col = []
        for elem_idx in range(self.get_num_elements()):
            coords = self.get_element_nodes(elem_idx)
            dofs = self.get_element_global_dofs(elem_idx).ravel()
            Ke = np.zeros((48, 48))
            for i, xi in enumerate(self.gauss_points):
                for j, eta in enumerate(self.gauss_points):
                    for k, zeta in enumerate(self.gauss_points):
                        w = self.gauss_weights[i] * self.gauss_weights[j] * self.gauss_weights[k]
                        B, detJ, _, _ = self.build_B_matrix(xi, eta, zeta, coords)
                        w_det = w * detJ
                        Ke += B.T @ self.D @ B * w_det
            for a in range(48):
                for b in range(48):
                    K_row.append(dofs[a])
                    K_col.append(dofs[b])
                    K_data.append(Ke[a, b])
        K = coo_matrix((K_data, (K_row, K_col)), shape=(num_dofs, num_dofs))
        return K.tocsr()

    def build_force_vector(self):
        num_dofs = self.grid_points.shape[0] * self.ndof_per_node
        F = np.zeros(num_dofs)
        body_force = np.array([0, 0, self.gamma])
        for elem_idx in range(self.get_num_elements()):
            coords = self.get_element_nodes(elem_idx)
            dofs = self.get_element_global_dofs(elem_idx).ravel()
            Fe = np.zeros(48)
            for i, xi in enumerate(self.gauss_points):
                for j, eta in enumerate(self.gauss_points):
                    for k, zeta in enumerate(self.gauss_points):
                        w = self.gauss_weights[i] * self.gauss_weights[j] * self.gauss_weights[k]
                        J = self.get_hermite_jacobian(xi, eta, zeta, coords)
                        detJ = np.linalg.det(J)
                        NH, _, _, _ = self.get_hermite_shape_functions_and_derivatives(xi, eta, zeta)
                        N_disp = np.zeros((3, 48))
                        for n in range(8):
                            col = 6 * n
                            N_disp[:, col:col+3] = NH[n] * np.eye(3)
                        Fe += N_disp.T @ body_force * w * detJ
            F[dofs] += Fe
        surface_force = np.array([0, 0, self.w])
        for elem_idx in self._elems_max_z:
            coords = self.get_element_nodes(elem_idx)
            dofs = self.get_element_global_dofs(elem_idx).ravel()
            Fe = np.zeros(48)
            face_nodes = [4,5,6,7]
            for i, xi in enumerate(self.gauss_points):
                for j, eta in enumerate(self.gauss_points):
                    zeta = 1.0
                    w = self.gauss_weights[i] * self.gauss_weights[j]
                    NH, _, _, _ = self.get_hermite_shape_functions_and_derivatives(xi, eta, zeta)
                    J = self.get_hermite_jacobian(xi, eta, zeta, coords)
                    dx_dxi = J[0, :]
                    dx_deta = J[1, :]
                    cross = np.cross(dx_dxi, dx_deta)
                    detJs = np.linalg.norm(cross)
                    N_disp = np.zeros((3, 48))
                    for n in range(8):
                        col = 6 * n
                        N_disp[:, col:col+3] = NH[n] * np.eye(3)
                    Fe += N_disp.T @ surface_force * w * detJs
            F[dofs] += Fe
        return F

    def solve(self):
        K = self.build_global_K()
        F = self.build_force_vector()
        fixed_dof = []
        for node in self._nodes_min_z:
            node_dof = self.global_dofs[node]
            fixed_dof.extend(node_dof)
        free_dof = [d for d in range(K.shape[0]) if d not in fixed_dof]
        K_ff = K[np.ix_(free_dof, free_dof)]
        F_f = F[free_dof]
        U_f = spsolve(K_ff, F_f)
        U = np.zeros(K.shape[0])
        U[free_dof] = U_f
        return U

    def get_strain_at_center(self, elem_idx, U):
        coords = self.get_element_nodes(elem_idx)
        dofs = self.get_element_global_dofs(elem_idx).ravel()
        U_e = U[dofs]
        xi, eta, zeta = 0.0, 0.0, 0.0
        B, _, _, _ = self.build_B_matrix(xi, eta, zeta, coords)
        return B @ U_e
    
    def get_all_points(self):
        return self.grid_points
    
    def elastic_matrix(self):
        l = self.lambda_
        m = 2 * self.mu
        D = np.array([
            [l + m, l, l, 0, 0, 0],
            [l, l + m, l, 0, 0, 0],
            [l, l, l + m, 0, 0, 0],
            [0, 0, 0, self.mu, 0, 0],
            [0, 0, 0, 0, self.mu, 0],
            [0, 0, 0, 0, 0, self.mu]
        ])
        return D

    def get_interpolated_von_mises(self, x, y, z, U):
        if not hasattr(self, 'vms'):
            self.vms = np.zeros(self.get_num_elements())
            D = self.elastic_matrix()
            for elem_idx in range(self.get_num_elements()):
                strain = self.get_strain_at_center(elem_idx, U)
                stress = D @ strain
                sx, sy, sz, txy, tyz, txz = stress
                vm = np.sqrt(0.5 * ((sx - sy)**2 + (sy - sz)**2 + (sz - sx)**2) + 3 * (txy**2 + tyz**2 + txz**2))
                self.vms[elem_idx] = vm
        
        dx = self.Lx / (self.nx - 1)
        dy = self.Ly / (self.ny - 1)
        dz = self.Lz / (self.nz - 1)
        x_centers = np.linspace(-self.Lx/2 + dx/2, self.Lx/2 - dx/2, self.nx - 1)
        y_centers = np.linspace(-self.Ly/2 + dy/2, self.Ly/2 - dy/2, self.ny - 1)
        z_centers = np.linspace(dz/2, self.Lz - dz/2, self.nz - 1)
        
        vms_grid = self.vms.reshape((self.nx - 1, self.ny - 1, self.nz - 1))
        
        points = (x_centers, y_centers, z_centers)
        return interpn(points, vms_grid, (x, y, z), method='linear', bounds_error=False, fill_value=None)
    

    
