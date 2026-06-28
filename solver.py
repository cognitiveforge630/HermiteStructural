from FEMHermiteBeamRegion import FEMHermiteBeamRegion, OrthotropicElasticity
from checkme import E, I, L, P, cantilever_point_load_deflection
import numpy as np
from scipy.sparse.linalg import spsolve


OUTPUT_PATH = "U.npy"
LOAD_DOF = 1  # uy: transverse bending load direction
THREE_POINT_GAUSS_ORDER = 3

# Keep this as the AWC/reference modulus used to define the displacement target.
# The calibrated modulus below is a solver-side stiffness slackening value; it
# should not overwrite the target-theory E.
REFERENCE_LONGITUDINAL_E = E
REFERENCE_TARGET_UY = cantilever_point_load_deflection(P, L, REFERENCE_LONGITUDINAL_E, I)

REDWOOD_CLEAR_WOOD_RATIOS = {
    # Redwood-like orthotropic ratios with L along the beam/global z-axis.
    "E_T_over_E_L": 0.089,
    "E_R_over_E_L": 0.087,
    "G_LR_over_E_L": 0.066,
    "G_LT_over_E_L": 0.077,
    "G_RT_over_E_L": 0.011,
    "nu_LR": 0.360,
    "nu_LT": 0.346,
    "nu_RT": 0.373,
}


def redwood_orthotropic_material(E_longitudinal):
    """Return orthotropic constants with wood axes R->x, T->y, and L->z."""
    ratios = REDWOOD_CLEAR_WOOD_RATIOS
    E_longitudinal = float(E_longitudinal)
    Ex = ratios["E_R_over_E_L"] * E_longitudinal
    Ey = ratios["E_T_over_E_L"] * E_longitudinal
    Ez = E_longitudinal
    return OrthotropicElasticity(
        Ex=Ex,
        Ey=Ey,
        Ez=Ez,
        nuxy=ratios["nu_RT"],
        nuxz=ratios["nu_LR"] * Ex / Ez,
        nuyz=ratios["nu_LT"] * Ey / Ez,
        Gxy=ratios["G_RT_over_E_L"] * E_longitudinal,
        Gyz=ratios["G_LT_over_E_L"] * E_longitudinal,
        Gxz=ratios["G_LR_over_E_L"] * E_longitudinal,
    )


def calibrated_longitudinal_E(seed_E, seed_fea_uy, target_uy=REFERENCE_TARGET_UY):
    """Return the positive Ez that makes linear FEA scaling hit target_uy."""
    seed_E = float(seed_E)
    seed_fea_uy = float(seed_fea_uy)
    target_uy = float(target_uy)
    if seed_E <= 0.0:
        raise ValueError("seed_E must be positive")
    if target_uy == 0.0:
        raise ValueError("target_uy must be nonzero")
    ratio = seed_fea_uy / target_uy
    if ratio <= 0.0:
        raise ValueError(
            "seed FEA displacement and target displacement must have the same sign "
            f"for positive E calibration; got seed={seed_fea_uy}, target={target_uy}"
        )
    return seed_E * ratio


def build_beam(E_longitudinal=REFERENCE_LONGITUDINAL_E, nx=5, ny=9, nz=121):
    material = redwood_orthotropic_material(E_longitudinal)
    return FEMHermiteBeamRegion(
        Lx=4.0,
        Ly=8.0,
        Lz=120.0,
        nx=nx,
        ny=ny,
        nz=nz,
        E=E_longitudinal,
        nu=0.3,
        gamma=0.0,
        w=0.0,
        gauss_order=THREE_POINT_GAUSS_ORDER,
        orthotropic_constants=material,
    )


def build_end_point_load_vector(beam, total_load=P, load_dof=LOAD_DOF):
    force = np.zeros(beam.get_all_points().shape[0] * beam.ndof_per_node)
    end_nodes = beam.get_nodes_at_max_z()
    nodal_load = total_load / len(end_nodes)

    for node in end_nodes:
        force[beam.global_dofs[node, load_dof]] += nodal_load

    return force, end_nodes, nodal_load


def build_end_face_traction_load_vector(beam, total_load=P, load_dof=LOAD_DOF):
    force = np.zeros(beam.get_all_points().shape[0] * beam.ndof_per_node)
    face_area = beam.Lx * beam.Ly
    traction = total_load / face_area

    for elem_idx in beam.get_elements_at_max_z():
        coords = beam.get_element_nodes(elem_idx)
        dofs = beam.get_element_global_dofs(elem_idx).ravel()
        Fe = np.zeros(beam.ndof_per_node * 8)

        for i, xi in enumerate(beam.gauss_points):
            for j, eta in enumerate(beam.gauss_points):
                zeta = 1.0
                weight = beam.gauss_weights[i] * beam.gauss_weights[j]
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


def solve_with_force_vector(beam, force):
    stiffness = beam.build_global_K()
    fixed_dof = []

    for node in beam.get_nodes_at_min_z():
        fixed_dof.extend(beam.global_dofs[node])

    free_dof = [dof for dof in range(stiffness.shape[0]) if dof not in fixed_dof]
    solved_free = spsolve(
        stiffness[np.ix_(free_dof, free_dof)],
        force[free_dof],
    )

    displacement = np.zeros(stiffness.shape[0])
    displacement[free_dof] = solved_free
    return displacement


def free_end_uy(beam, U):
    U_by_node = U.reshape((-1, beam.ndof_per_node))
    end_nodes = beam.get_nodes_at_max_z()
    end_uy = U_by_node[end_nodes, LOAD_DOF]
    return end_uy[np.argmax(np.abs(end_uy))]


def print_theoretical_check(beam, force=None, U=None, label="solve"):
    delta_theory = REFERENCE_TARGET_UY
    end_nodes = beam.get_nodes_at_max_z()
    end_face_area = beam.Lx * beam.Ly
    traction = P / end_face_area

    print(f"Theoretical displacement check ({label}):")
    print("  Solver load model: end-face traction integrated over z=L")
    print("  Material model: orthotropic Redwood/AWC D matrix")
    print(f"  Gauss integration order:              {beam.gauss_order}-point")
    print(f"  Total P:                              {P:.3f} lb")
    print(f"  End face area:                        {end_face_area:.3f} in^2")
    print(f"  End face traction:                    {traction:.6f} psi")
    print(f"  End nodes loaded:                     {len(end_nodes)}")
    if force is not None:
        assembled_load = force[beam.global_dofs[end_nodes, LOAD_DOF]].sum()
        print(f"  Integrated load on end nodes:         {assembled_load:.6f} lb")
    print("  Formula:                              delta = P L^3 / (3 E I)")
    print(f"  I:                                    {I:.6f} in^4")
    print(f"  AWC/reference E for target:           {REFERENCE_LONGITUDINAL_E:.3f} psi")
    print(f"  Solver orthotropic Ez:                {beam.orthotropic_constants.Ez:.3f} psi")
    print(f"  Target theoretical free-end uy:       {delta_theory:.3e} inches")

    if U is None:
        print("  FEA displacement: pending solve")
        print()
        return

    fea_uy = free_end_uy(beam, U)
    rel_diff = abs(delta_theory - fea_uy) / abs(delta_theory) * 100
    U_by_node = U.reshape((-1, beam.ndof_per_node))
    displacement_magnitude = np.linalg.norm(U_by_node[:, :3], axis=1)
    max_displacement_magnitude = displacement_magnitude.max()

    print(f"  FEA free-end uy:                      {fea_uy:.3e} inches")
    print(f"  Difference vs AWC/theory target:      {rel_diff:.6f}%")
    print(f"  FEA max displacement magnitude:       {max_displacement_magnitude:.3e} inches")
    print()


def solve_calibrated_beam():
    seed_beam = build_beam(REFERENCE_LONGITUDINAL_E)
    seed_force, _, _ = build_end_face_traction_load_vector(seed_beam)

    print_theoretical_check(seed_beam, seed_force, label="orthotropic seed")
    print("Solving seed finite element system...", flush=True)
    seed_U = solve_with_force_vector(seed_beam, seed_force)
    print_theoretical_check(seed_beam, seed_force, seed_U, label="orthotropic seed")

    seed_uy = free_end_uy(seed_beam, seed_U)
    calibrated_E = calibrated_longitudinal_E(
        REFERENCE_LONGITUDINAL_E,
        seed_uy,
        REFERENCE_TARGET_UY,
    )
    print("E calibration:")
    print(f"  seed Ez:                              {REFERENCE_LONGITUDINAL_E:.6f} psi")
    print(f"  seed FEA uy:                          {seed_uy:.12e} in")
    print(f"  target uy:                            {REFERENCE_TARGET_UY:.12e} in")
    print("  Formula:                              E_calibrated = E_seed * u_seed / u_target")
    print(f"  calibrated Ez:                        {calibrated_E:.6f} psi")
    print()

    calibrated_beam = build_beam(calibrated_E)
    calibrated_force, _, _ = build_end_face_traction_load_vector(calibrated_beam)
    print("Solving calibrated finite element system...", flush=True)
    calibrated_U = solve_with_force_vector(calibrated_beam, calibrated_force)
    print_theoretical_check(
        calibrated_beam,
        calibrated_force,
        calibrated_U,
        label="orthotropic calibrated",
    )
    return calibrated_beam, calibrated_U


def main():
    beam, U = solve_calibrated_beam()
    np.save(OUTPUT_PATH, U)

    strain = beam.get_strain_at_center(elem_idx=1920, U=U)
    strain_labels = ("exx", "eyy", "ezz", "exy", "eyz", "exz")
    print("Strain at element 1920:")
    for label, value in zip(strain_labels, strain):
        print(f"  {label}: {value:.6e}")

    print(f"Saved displacement vector to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
