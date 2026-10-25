from FEMHermiteBeamRegion import FEMHermiteBeamRegion, OrthotropicElasticity
from reference_hex8 import CauchyHex8Reference
from reference_hex27 import CauchyHex27Reference
from enriched_hex8 import AffineEnrichedHex8
from recovered_rotation_hex8 import RecoveredRotationHex8
from quadratic_rotation_hex8 import QuadraticRotationHex8
from cubic_rotation_hex8 import CubicRotationHex8
from checkme import E, I, L, P, cantilever_point_load_deflection
import numpy as np
from scipy.sparse.linalg import spsolve
from linear_solvers import solve_sparse_system


OUTPUT_PATH = "U.npy"
LOAD_DOF = 1  # uy: transverse bending load direction
THREE_POINT_GAUSS_ORDER = 3

# Keep this as the AWC/reference modulus used to define the displacement target.
# The calibrated modulus below is a solver-side stiffness slackening value; it
# should not overwrite the target-theory E.
REFERENCE_LONGITUDINAL_E = E
REFERENCE_TARGET_UY = cantilever_point_load_deflection(P, L, REFERENCE_LONGITUDINAL_E, I)

# Redwood specific weight; configurable independently of applied end-face traction.
DEFAULT_WEIGHT_DENSITY = 0.0174  # lbf/in^3; body force is in negative y
CROSS_SECTION_WIDTH_IN = 4.0
CROSS_SECTION_DEPTH_IN = 8.0
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


def build_beam(E_longitudinal=REFERENCE_LONGITUDINAL_E, nx=5, ny=9, nz=121, *,
               formulation='legacy_hermite'):
    if formulation not in ('legacy_hermite', 'cauchy_hex8_reference', 'cauchy_hex27_reference', 'affine_enriched_hex8', 'recovered_rotation_hermite', 'quadratic_rotation_hermite', 'cubic_rotation_hermite'):
        raise ValueError(f'Unknown formulation: {formulation}')
    element_class = {'legacy_hermite': FEMHermiteBeamRegion,
                     'cauchy_hex8_reference': CauchyHex8Reference,
                     'cauchy_hex27_reference': CauchyHex27Reference,
                     'affine_enriched_hex8': AffineEnrichedHex8,
                     'recovered_rotation_hermite': RecoveredRotationHex8,
                     'quadratic_rotation_hermite': QuadraticRotationHex8,
                     'cubic_rotation_hermite': CubicRotationHex8}[formulation]
    material = redwood_orthotropic_material(E_longitudinal)
    return element_class(
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
        gauss_order=(4 if formulation in ("affine_enriched_hex8", "quadratic_rotation_hermite", "cubic_rotation_hermite") else THREE_POINT_GAUSS_ORDER),
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


def solve_with_prescribed_displacement(beam, force, prescribed_translations=None,
                                       prescribed_rotations=None,*,stiffness=None,linear_solver='scipy'):
    """Solve with z=0 Dirichlet data; no virtual environment or calibration.

    With the experimental bubble-enriched Hex8, the last three nodal DOFs are
    *interior enrichment coefficients*, not boundary rotations. Consequently
    only actual displacement DOFs are fixed along the Dirichlet face.
    Existing legacy and Cauchy options retain their previous BC behavior.
    """
    stiffness = beam.build_global_K() if stiffness is None else stiffness
    force = np.asarray(force, dtype=float)
    if force.shape != (stiffness.shape[0],):
        raise ValueError('Force vector does not match global degrees of freedom')
    fixed = []
    u_known = np.zeros(stiffness.shape[0])
    for node in beam.get_nodes_at_min_z():
        nfixed = (3 if isinstance(beam, AffineEnrichedHex8)
                  else beam.ndof_per_node)
        ids = beam.global_dofs[node, :nfixed]
        fixed.extend(ids.tolist())
        if prescribed_translations is not None:
            values = np.asarray(prescribed_translations(beam.grid_points[node]), dtype=float)
            if values.shape != (3,):
                raise ValueError('Dirichlet field must return three translations')
            u_known[beam.global_dofs[node,:3]] = values
        if prescribed_rotations is not None and beam.ndof_per_node == 6 and not isinstance(beam, AffineEnrichedHex8):
            rotations = np.asarray(prescribed_rotations(beam.grid_points[node]), dtype=float)
            if rotations.shape != (3,):
                raise ValueError('Dirichlet rotation field must return three components')
            u_known[beam.global_dofs[node,3:6]] = rotations
    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(len(u_known)), fixed)
    if len(free):
        rhs = force[free] - stiffness[free][:,fixed] @ u_known[fixed]
        u_known[free],solve_info = solve_sparse_system(stiffness[free][:,free],rhs,linear_solver)
        beam.last_linear_solve=solve_info
    if not np.all(np.isfinite(u_known)):
        raise ValueError('Solve yielded nonfinite displacement')
    return u_known


def solve_with_force_vector(beam, force):
    if isinstance(beam, RecoveredRotationHex8):
        # Nonlocal recovered vertex derivatives make the full face trace depend
        # on some free neighbor DOFs. Fixing face nodes alone is insufficient.
        from face_dirichlet import solve_with_face_dirichlet
        U,_=solve_with_face_dirichlet(beam,force,lambda x:(0.,0.,0.),
                                      lambda x:(0.,0.,0.))
        return U
    return solve_with_prescribed_displacement(beam, force)


def free_end_uy(beam, U):
    U_by_node = U.reshape((-1, beam.ndof_per_node))
    end_nodes = beam.get_nodes_at_max_z()
    end_uy = U_by_node[end_nodes, LOAD_DOF]
    return end_uy[np.argmax(np.abs(end_uy))]


def print_theoretical_check(beam, force=None, U=None, label="solve", weight_density=DEFAULT_WEIGHT_DENSITY):
    """Report the SAME gravity-inclusive reference used by the FE load vector."""
    reference = awc_reference_components(weight_density)
    delta_theory = reference['total_uy_in']
    end_nodes = beam.get_nodes_at_max_z()
    end_face_area = beam.Lx * beam.Ly
    traction = P / end_face_area
    print(f"Theoretical displacement check ({label}):")
    print("  Load model: Gauss-integrated end traction + volume self-weight")
    print("  Material model: orthotropic Redwood/AWC D matrix")
    print(f"  Gauss integration order:              {beam.gauss_order}-point")
    print(f"  End-face force:                       {P:.6f} lbf")
    print(f"  End-face area:                        {end_face_area:.6f} in^2")
    print(f"  End-face traction:                    {traction:.6f} psi")
    print(f"  Self-weight force:                    {-reference['total_weight_lbf']:.6f} lbf")
    print(f"  End nodes:                            {len(end_nodes)}")
    if force is not None:
        translation_resultant = float(force[beam.global_dofs[:, LOAD_DOF]].sum())
        print(f"  Assembled y-translation resultant:    {translation_resultant:.6f} lbf")
        print("  Note: end-node subtotal is NOT the end-face traction resultant")
    print("  Formula: PL^3/(3EI) + wL^4/(8EI), signed")
    print(f"  I:                                    {I:.6f} in^4")
    print(f"  AWC reference / solver Ez:            {beam.orthotropic_constants.Ez:.3f} psi")
    print(f"  AWC end-load uy:                      {reference['end_load_uy_in']:.9f} in")
    print(f"  AWC gravity uy:                       {reference['self_weight_uy_in']:.9f} in")
    print(f"  AWC combined uy:                      {delta_theory:.9f} in")
    if U is None:
        print("  FEM displacement: pending solve\n")
        return
    fea_uy = free_end_uy(beam, U)
    rel_diff = abs(delta_theory-fea_uy)/abs(delta_theory)*100
    disp = np.linalg.norm(U.reshape((-1,beam.ndof_per_node))[:,:3],axis=1)
    print(f"  FEA free-end uy:                      {fea_uy:.9f} in")
    print(f"  Difference vs combined AWC:           {rel_diff:.6f}%")
    print(f"  FEA max displacement magnitude:       {disp.max():.9f} in\n")


def build_self_weight_load_vector(beam, weight_density=DEFAULT_WEIGHT_DENSITY):
    """Consistent body-force vector using the element's Hermite displacement operator.

    Orthotropic material orientation does not change the global downward-y force.
    weight_density is specific weight in lbf/in^3, not mass density.
    """
    if not np.isfinite(weight_density) or weight_density <= 0:
        raise ValueError('Tutorial 001 requires a positive finite Redwood specific weight')
    force = np.zeros(beam.get_all_points().shape[0] * beam.ndof_per_node)
    body = np.array([0., -weight_density, 0.])
    for e in range(beam.get_num_elements()):
        coords = beam.get_element_nodes(e)
        dofs = beam.get_element_global_dofs(e).ravel()
        fe = np.zeros(len(dofs))
        for i, xi in enumerate(beam.gauss_points):
            for j, eta in enumerate(beam.gauss_points):
                for k, zeta in enumerate(beam.gauss_points):
                    N, _, _, _ = beam.get_hermite_displacement_matrices(xi, eta, zeta, coords)
                    J = beam.get_hermite_jacobian(xi, eta, zeta, coords)
                    detJ = float(np.linalg.det(J))
                    if detJ <= 0: raise ValueError('Invalid Hermite volume Jacobian')
                    fe += N.T @ body * detJ * beam.gauss_weights[i] * beam.gauss_weights[j] * beam.gauss_weights[k]
        force[dofs] += fe
    return force


def awc_reference_components(weight_density=DEFAULT_WEIGHT_DENSITY):
    """Signed Euler-Bernoulli cantilever deflections (inches), identical load inputs.

    These are classical beam-theory references, not an independent certification.
    The transverse loading axis is global -y and the beam axis is global +z.
    """
    if not np.isfinite(weight_density) or weight_density <= 0:
        raise ValueError('Tutorial 001 requires positive finite self-weight')
    if not np.isclose(I, CROSS_SECTION_WIDTH_IN*CROSS_SECTION_DEPTH_IN**3/12):
        raise ValueError('The reference inertia does not match the FEM cross section')
    end = cantilever_point_load_deflection(P, L, REFERENCE_LONGITUDINAL_E, I)
    q = weight_density * CROSS_SECTION_WIDTH_IN * CROSS_SECTION_DEPTH_IN
    gravity = -q * L**4 / (8.0 * REFERENCE_LONGITUDINAL_E * I)
    return {'end_load_uy_in': float(end), 'self_weight_uy_in': float(gravity),
            'total_uy_in': float(end + gravity),
            'weight_per_length_lbf_in': float(q),
            'total_weight_lbf': float(q * L)}


def combined_reference_tip(weight_density=DEFAULT_WEIGHT_DENSITY):
    return awc_reference_components(weight_density)['total_uy_in']

def solve_referenced_beam(weight_density=DEFAULT_WEIGHT_DENSITY):
    """Run one Hermite solve with fixed, documented material constants."""
    awc_reference_components(weight_density)
    beam = build_beam(REFERENCE_LONGITUDINAL_E)
    traction_force, _, _ = build_end_face_traction_load_vector(beam)
    gravity_force = build_self_weight_load_vector(beam, weight_density)
    force = traction_force + gravity_force
    end_resultant = float(traction_force[beam.global_dofs[:,LOAD_DOF]].sum())
    gravity_resultant = float(gravity_force[beam.global_dofs[:,LOAD_DOF]].sum())
    expected_gravity = -weight_density * beam.Lx * beam.Ly * beam.Lz
    if not np.isclose(end_resultant, P, atol=1e-6, rtol=1e-9):
        raise RuntimeError(f'End traction resultant mismatch: {end_resultant} vs {P}')
    if not np.isclose(gravity_resultant, expected_gravity, atol=1e-6, rtol=1e-9):
        raise RuntimeError(f'Self-weight resultant mismatch: {gravity_resultant} vs {expected_gravity}')
    print(f'Gauss-integrated end-face force: {end_resultant:.6f} lbf')
    print(f'Gauss-integrated body force:     {gravity_resultant:.6f} lbf')
    print(f'Combined applied force:         {end_resultant + gravity_resultant:.6f} lbf')
    print_theoretical_check(beam, force, label='fixed orthotropic modulus',weight_density=weight_density)
    print('Solving Hermite finite element system...',flush=True)
    U = solve_with_force_vector(beam, force)
    print_theoretical_check(beam, force, U,label='fixed orthotropic modulus',weight_density=weight_density)
    return beam, U


def main():
    import argparse, json
    parser = argparse.ArgumentParser(description='Historical orthotropic Hermite solve with integrated traction and gravity')
    parser.add_argument('--weight-density',type=float,default=DEFAULT_WEIGHT_DENSITY,
                        help='Redwood specific weight, lbf/in^3; must be positive (gravity mandatory)')
    args = parser.parse_args()
    beam, U = solve_referenced_beam(args.weight_density)
    np.save(OUTPUT_PATH, U)
    import json
    report = {'formulation':'Historical Hermite Hex8 6 DOF/node',
              'material':'Orthotropic Redwood (fixed longitudinal modulus)',
              'weight_density_lbf_in3':args.weight_density,
              'self_weight_y_lbf':-args.weight_density*4*8*L,
              'end_traction_y_lbf':P,
              'awc_point_load_only_uy_in':REFERENCE_TARGET_UY,
              'awc_point_plus_gravity_uy_in':combined_reference_tip(args.weight_density),
              'awc_components':awc_reference_components(args.weight_density),
              'awc_reference_longitudinal_E_psi':REFERENCE_LONGITUDINAL_E,
              'fem_longitudinal_E_psi':beam.orthotropic_constants.Ez,
              'end_face_area_in2':4*8,
              'beam_length_in':L,
              'fem_tip_uy_in':free_end_uy(beam,U),
              'end_face':'Gauss-integrated traction',
              'body_force':'Gauss-integrated self-weight',
              'loading_validation':'positive self-weight plus end-face traction',
              'calibrated_to_reference':False}
    from pathlib import Path
    Path('hermite_result.json').write_text(json.dumps(report,indent=2))

    strain = beam.get_strain_at_center(elem_idx=1920, U=U)
    strain_labels = ("exx", "eyy", "ezz", "exy", "eyz", "exz")
    print("Strain at element 1920:")
    for label, value in zip(strain_labels, strain):
        print(f"  {label}: {value:.6e}")

    print(f"Saved displacement vector to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
