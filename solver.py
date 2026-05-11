from FEMHermiteBeamRegion import FEMHermiteBeamRegion
from checkme import E, I, L, P, cantilever_point_load_deflection
import numpy as np
from scipy.sparse.linalg import spsolve


OUTPUT_PATH = "U.npy"
LOAD_DOF = 1  # uy: transverse bending load direction


def build_beam():
    return FEMHermiteBeamRegion(
        Lx=4.0,
        Ly=8.0,
        Lz=120.0,
        nx=5,
        ny=9,
        nz=121,
        E=2000000,
        nu=0.3,
        gamma=0.0,
        w=0.0,
    )


def build_end_point_load_vector(beam, total_load=P, load_dof=LOAD_DOF):
    force = np.zeros(beam.get_all_points().shape[0] * beam.ndof_per_node)
    end_nodes = beam.get_nodes_at_max_z()
    nodal_load = total_load / len(end_nodes)

    for node in end_nodes:
        force[beam.global_dofs[node, load_dof]] += nodal_load

    return force, end_nodes, nodal_load


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


def print_theoretical_check(beam, U=None):
    delta_theory = cantilever_point_load_deflection(P, L, E, I)
    end_nodes = beam.get_nodes_at_max_z()
    nodal_load = P / len(end_nodes)

    print("Theoretical displacement check:")
    print("  Solver load model: point load P distributed over z=L end nodes")
    print(f"  Total P:                              {P:.3f} lb")
    print(f"  End nodes loaded:                     {len(end_nodes)}")
    print(f"  Load per end node:                    {nodal_load:.6f} lb")
    print("  Formula:                              delta = P L^3 / (3 E I)")
    print(f"  I:                                    {I:.6f} in^4")
    print(f"  Theoretical free-end uy:              {delta_theory:.3e} inches")

    if U is None:
        print("  FEA displacement: pending solve")
        print()
        return

    U_by_node = U.reshape((-1, beam.ndof_per_node))
    end_uy = U_by_node[end_nodes, LOAD_DOF]
    fea_end_uy = end_uy[np.argmax(np.abs(end_uy))]
    rel_diff = abs(delta_theory - fea_end_uy) / abs(delta_theory) * 100
    displacement_magnitude = np.linalg.norm(U_by_node[:, :3], axis=1)
    max_displacement_magnitude = displacement_magnitude.max()

    print(f"  FEA free-end uy:                      {fea_end_uy:.3e} inches")
    print(f"  Difference vs theory:                 {rel_diff:.2f}%")
    print(f"  FEA max displacement magnitude:       {max_displacement_magnitude:.3e} inches")
    print()


def main():
    beam = build_beam()
    force, _, _ = build_end_point_load_vector(beam)

    print_theoretical_check(beam)
    print("Solving finite element system...", flush=True)
    U = solve_with_force_vector(beam, force)
    np.save(OUTPUT_PATH, U)
    print_theoretical_check(beam, U)

    strain = beam.get_strain_at_center(elem_idx=1920, U=U)
    strain_labels = ("exx", "eyy", "ezz", "exy", "eyz", "exz")
    print("Strain at element 1920:")
    for label, value in zip(strain_labels, strain):
        print(f"  {label}: {value:.6e}")

    print(f"Saved displacement vector to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
