"""Consistent end-face traction and volume self-weight integration.

These loads are built using each element's actual nonlocal displacement
operator, not a lumped force approximation.
"""
import numpy as np
from core.elements.linear import RecoveredRotationHex8
from core.reference.beam_theory import P
LOAD_DOF = 1

def build_end_face_traction_load_vector(beam, total_load=P, load_dof=LOAD_DOF):
    force = np.zeros(beam.get_all_points().shape[0] * beam.ndof_per_node)
    face_area = beam.Lx * beam.Ly
    traction = total_load / face_area

    for elem_idx in beam.get_elements_at_max_z():
        coords = beam.get_element_nodes(elem_idx)
        nonlocal_stencil = isinstance(beam, RecoveredRotationHex8)
        dofs = (beam.element_transform(elem_idx)[0] if nonlocal_stencil
                else beam.get_element_global_dofs(elem_idx).ravel())
        Fe = np.zeros(len(dofs))

        for i, xi in enumerate(beam.gauss_points):
            for j, eta in enumerate(beam.gauss_points):
                zeta = 1.0
                weight = beam.gauss_weights[i] * beam.gauss_weights[j]
                if nonlocal_stencil:
                    _, N_disp = beam.displacement_operator(elem_idx, (xi, eta, zeta))
                else:
                    N_disp, _, _, _ = beam.get_hermite_displacement_matrices(
                        xi, eta, zeta, coords
                    )
                J = beam.get_hermite_jacobian(xi, eta, zeta, coords)
                dx_dxi = J[0, :]
                dx_deta = J[1, :]
                detJs = np.linalg.norm(np.cross(dx_dxi, dx_deta))

                unit_traction = np.zeros(3)
                unit_traction[load_dof] = traction
                Fe += N_disp.T @ unit_traction * weight * detJs

        force[dofs] += Fe

    return force, beam.get_nodes_at_max_z(), traction


def gravity_load(beam, density):
    """Consistent nonlocal body force: integrate N_disp.T b over each volume."""
    f = np.zeros(beam.get_all_points().shape[0] * beam.ndof_per_node)
    body = np.array([0.0, -density, 0.0])
    for e in range(beam.get_num_elements()):
        coords = beam.get_element_nodes(e)
        dofs, _ = beam.element_transform(e)
        fe = np.zeros(len(dofs))
        for i, xi in enumerate(beam.gauss_points):
            for j, eta in enumerate(beam.gauss_points):
                for k, zeta in enumerate(beam.gauss_points):
                    _, N = beam.displacement_operator(e, (xi, eta, zeta))
                    J = beam.get_hermite_jacobian(xi, eta, zeta, coords)
                    detJ = float(np.linalg.det(J))
                    if detJ <= 0:
                        raise ValueError("Invalid volume Jacobian")
                    fe += N.T @ body * detJ * (
                        beam.gauss_weights[i] * beam.gauss_weights[j]
                        * beam.gauss_weights[k])
        f[dofs] += fe
    return f
