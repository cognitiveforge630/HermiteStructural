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
        self.D_stretch = self._build_D_stretch()
        self.D_curv = self.gamma_c * np.eye(9)
        self.gauss_points, self.gauss_weights = leggauss(3)

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

    def _build_D_stretch(self):
        D = np.zeros((9,9))
        for k in range(9):
            gamma_flat = np.zeros(9)
            gamma_flat[k] = 1
            gamma = gamma_flat.reshape(3,3)
            e = 0.5 * (gamma + gamma.T)
            a = 0.5 * (gamma - gamma.T)
            tr_e = np.trace(e)
            sigma_sym = self.lambda_ * tr_e * np.eye(3) + 2 * self.mu * e
            sigma_skew = 2 * self.mu_c * a
            sigma = sigma_sym + sigma_skew
            D[:,k] = sigma.ravel()
        return D

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

    def get_hermite_shapes_and_derivs(self, xi, eta, zeta):
        # Use element-node Hex8 translation shapes; the old cubic-only subset
        # over-stiffens cantilever bending by breaking the standard nodal basis.
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

    def build_global_K(self):
        num_dofs = self.grid_points.shape[0] * self.ndof_per_node
        K_data = []
        K_row = []
        K_col = []
        for elem_idx in range(self.get_num_elements()):
            coords = self.get_element_nodes(elem_idx)
            dofs = self.get_element_global_dofs(elem_idx).ravel()
            Ke = np.zeros((48, 48))
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        xi = self.gauss_points[i]
                        eta = self.gauss_points[j]
                        zeta = self.gauss_points[k]
                        w = self.gauss_weights[i] * self.gauss_weights[j] * self.gauss_weights[k]
                        dN_dn = self.hex8_shape_derivatives(xi, eta, zeta)
                        J = dN_dn @ coords
                        detJ = np.linalg.det(J)
                        J_inv = np.linalg.inv(J)
                        N_u, N_phi, dNu_dxi, dNu_deta, dNu_dzeta, dNp_dxi, dNp_deta, dNp_dzeta = self.get_hermite_shapes_and_derivs(xi, eta, zeta)
                        dNu_dx = J_inv[0, :] @ np.array([dNu_dxi, dNu_deta, dNu_dzeta])
                        dNu_dy = J_inv[1, :] @ np.array([dNu_dxi, dNu_deta, dNu_dzeta])
                        dNu_dz = J_inv[2, :] @ np.array([dNu_dxi, dNu_deta, dNu_dzeta])
                        dNp_dx = J_inv[0, :] @ np.array([dNp_dxi, dNp_deta, dNp_dzeta])
                        dNp_dy = J_inv[1, :] @ np.array([dNp_dxi, dNp_deta, dNp_dzeta])
                        dNp_dz = J_inv[2, :] @ np.array([dNp_dxi, dNp_deta, dNp_dzeta])
                        B_gamma = np.zeros((9, 48))
                        B_curv = np.zeros((9, 48))
                        for n in range(8):
                            col_trans = 6 * n
                            col_rot = 6 * n + 3
                            B_gamma[0, col_trans + 0] = dNu_dx[n]
                            B_gamma[1, col_trans + 0] = dNu_dy[n]
                            B_gamma[2, col_trans + 0] = dNu_dz[n]
                            B_gamma[3, col_trans + 1] = dNu_dx[n]
                            B_gamma[4, col_trans + 1] = dNu_dy[n]
                            B_gamma[5, col_trans + 1] = dNu_dz[n]
                            B_gamma[6, col_trans + 2] = dNu_dx[n]
                            B_gamma[7, col_trans + 2] = dNu_dy[n]
                            B_gamma[8, col_trans + 2] = dNu_dz[n]
                            B_gamma[1, col_rot + 2] = N_phi[n]
                            B_gamma[2, col_rot + 1] = -N_phi[n]
                            B_gamma[3, col_rot + 2] = -N_phi[n]
                            B_gamma[5, col_rot + 0] = N_phi[n]
                            B_gamma[6, col_rot + 1] = N_phi[n]
                            B_gamma[7, col_rot + 0] = -N_phi[n]
                            B_curv[0, col_rot + 0] = dNp_dx[n]
                            B_curv[1, col_rot + 0] = dNp_dy[n]
                            B_curv[2, col_rot + 0] = dNp_dz[n]
                            B_curv[3, col_rot + 1] = dNp_dx[n]
                            B_curv[4, col_rot + 1] = dNp_dy[n]
                            B_curv[5, col_rot + 1] = dNp_dz[n]
                            B_curv[6, col_rot + 2] = dNp_dx[n]
                            B_curv[7, col_rot + 2] = dNp_dy[n]
                            B_curv[8, col_rot + 2] = dNp_dz[n]
                        w_det = w * detJ
                        Ke += B_gamma.T @ self.D_stretch @ B_gamma * w_det
                        Ke += B_curv.T @ self.D_curv @ B_curv * w_det
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
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        xi = self.gauss_points[i]
                        eta = self.gauss_points[j]
                        zeta = self.gauss_points[k]
                        w = self.gauss_weights[i] * self.gauss_weights[j] * self.gauss_weights[k]
                        dN_dn = self.hex8_shape_derivatives(xi, eta, zeta)
                        J = dN_dn @ coords
                        detJ = np.linalg.det(J)
                        N_u, _, _, _, _, _, _, _ = self.get_hermite_shapes_and_derivs(xi, eta, zeta)
                        N_trans = np.zeros((3, 48))
                        for n in range(8):
                            col = 6*n
                            N_trans[:, col:col+3] = N_u[n] * np.eye(3)
                        Fe += N_trans.T @ body_force * w * detJ
            F[dofs] += Fe
        surface_force = np.array([0, 0, self.w])
        for elem_idx in self._elems_max_z:
            coords = self.get_element_nodes(elem_idx)
            dofs = self.get_element_global_dofs(elem_idx).ravel()
            Fe = np.zeros(48)
            face_nodes = [4,5,6,7]
            for i in range(3):
                for j in range(3):
                    xi = self.gauss_points[i]
                    eta = self.gauss_points[j]
                    zeta = 1.0
                    w = self.gauss_weights[i] * self.gauss_weights[j]
                    N_u, _, _, _, _, _, _, _ = self.get_hermite_shapes_and_derivs(xi, eta, zeta)
                    dN_dn = self.hex8_shape_derivatives(xi, eta, zeta)
                    J = dN_dn @ coords
                    dx_dxi = J[:,0]
                    dx_deta = J[:,1]
                    cross = np.cross(dx_dxi, dx_deta)
                    detJs = np.linalg.norm(cross)
                    N_trans = np.zeros((3, 48))
                    for nn in face_nodes:
                        col = 6*nn
                        N_trans[:, col:col+3] = N_u[nn] * np.eye(3)
                    Fe += N_trans.T @ surface_force * w * detJs
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
        U_trans = U_e.reshape(8,6)[:,:3]
        U_rot = U_e.reshape(8,6)[:,3:]
        xi, eta, zeta = 0.0, 0.0, 0.0
        N_u, N_phi, dNu_dxi, dNu_deta, dNu_dzeta, _, _, _ = self.get_hermite_shapes_and_derivs(xi, eta, zeta)
        dN_dn = self.hex8_shape_derivatives(xi, eta, zeta)
        J = dN_dn @ coords
        J_inv = np.linalg.inv(J)
        dNu_dx = dNu_dxi * J_inv[0,0] + dNu_deta * J_inv[0,1] + dNu_dzeta * J_inv[0,2]
        dNu_dy = dNu_dxi * J_inv[1,0] + dNu_deta * J_inv[1,1] + dNu_dzeta * J_inv[1,2]
        dNu_dz = dNu_dxi * J_inv[2,0] + dNu_deta * J_inv[2,1] + dNu_dzeta * J_inv[2,2]
        dNu_dxyz = np.vstack((dNu_dx, dNu_dy, dNu_dz))
        du_dxyz = U_trans.T @ dNu_dxyz.T
        phi = U_rot.T @ N_phi
        skew_phi = np.array([[0, -phi[2], phi[1]],
                             [phi[2], 0, -phi[0]],
                             [-phi[1], phi[0], 0]])
        gamma = du_dxyz - skew_phi
        e_sym = 0.5 * (gamma + gamma.T)
        strain_voigt = np.array([e_sym[0,0], e_sym[1,1], e_sym[2,2], e_sym[0,1], e_sym[1,2], e_sym[0,2]])
        return strain_voigt
    
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
    

    
