"""Element-level diagnostics for the Hermitian Hexa8 formulation.

These checks are intentionally small and local. They help separate shape
function/Jacobian/B-matrix problems from global cantilever modeling choices.
"""

import numpy as np

from FEMHermiteBeamRegion import FEMHermiteBeamRegion


NODE_NATURAL_COORDS = np.array(
    [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ]
)


def build_single_element():
    return FEMHermiteBeamRegion(
        Lx=4.0,
        Ly=8.0,
        Lz=1.0,
        nx=2,
        ny=2,
        nz=2,
        E=2_000_000.0,
        nu=0.3,
        gamma=0.0,
        w=0.0,
    )


def physical_point_from_mapping(beam, coords, xi, eta, zeta):
    N = beam.hex8_shape_functions(xi, eta, zeta)
    return N @ coords


def finite_difference_jacobian(beam, coords, xi, eta, zeta, eps=1e-6):
    center = np.array([xi, eta, zeta], dtype=float)
    jac = np.zeros((3, 3))

    for axis in range(3):
        step = np.zeros(3)
        step[axis] = eps
        plus = physical_point_from_mapping(beam, coords, *(center + step))
        minus = physical_point_from_mapping(beam, coords, *(center - step))
        jac[axis, :] = (plus - minus) / (2.0 * eps)

    return jac


def report(name, value, tolerance):
    status = "PASS" if value <= tolerance else "FAIL"
    print(f"{status:4} {name:<42} {value:.6e} <= {tolerance:.1e}")
    return status == "PASS"


def check_shape_identities(beam):
    print("\nShape function identities")
    passed = True

    max_nh_node_error = 0.0
    max_rh_node_value = 0.0
    max_rh_node_derivative_error = 0.0

    for node_idx, natural in enumerate(NODE_NATURAL_COORDS):
        NH, RH, _, dRH = beam.get_hermite_shape_functions_and_derivatives(*natural)
        expected = np.zeros(8)
        expected[node_idx] = 1.0
        max_nh_node_error = max(max_nh_node_error, np.max(np.abs(NH - expected)))
        max_rh_node_value = max(max_rh_node_value, np.max(np.abs(RH)))
        for axis in range(3):
            derivative_expected = np.zeros((3, 8))
            derivative_expected[axis, node_idx] = 1.0
            max_rh_node_derivative_error = max(
                max_rh_node_derivative_error,
                np.max(np.abs(dRH[axis] - derivative_expected)),
            )

    sample_points = [
        (0.0, 0.0, 0.0),
        (-0.3, 0.2, 0.7),
        (beam.gauss_points[0], beam.gauss_points[-1], beam.gauss_points[0]),
    ]
    max_partition_error = 0.0
    max_derivative_sum = 0.0
    for point in sample_points:
        NH, RH, dNH, dRH = beam.get_hermite_shape_functions_and_derivatives(*point)
        max_partition_error = max(max_partition_error, abs(np.sum(NH) - 1.0))
        max_derivative_sum = max(max_derivative_sum, np.max(np.abs(np.sum(dNH, axis=1))))

    passed &= report("NH nodal Kronecker delta", max_nh_node_error, 1e-12)
    passed &= report("RH nodal displacement contribution", max_rh_node_value, 1e-12)
    passed &= report("RH nodal first-derivative delta", max_rh_node_derivative_error, 1e-12)
    passed &= report("sum(NH) partition of unity", max_partition_error, 1e-12)
    passed &= report("sum(dNH) natural derivative", max_derivative_sum, 1e-12)

    return passed


def check_jacobian(beam, coords):
    print("\nJacobian consistency")
    passed = True

    max_error = 0.0
    for xi in beam.gauss_points:
        for eta in beam.gauss_points:
            for zeta in beam.gauss_points:
                analytic = beam.get_hermite_jacobian(xi, eta, zeta, coords)
                numeric = finite_difference_jacobian(beam, coords, xi, eta, zeta)
                max_error = max(max_error, np.max(np.abs(analytic - numeric)))

    passed &= report("analytic J vs finite difference J", max_error, 1e-8)
    return passed


def nodal_dofs_from_affine_field(coords, gradient, offset=None):
    if offset is None:
        offset = np.zeros(3)

    dofs = np.zeros(48)
    for node, xyz in enumerate(coords):
        displacement = offset + gradient @ xyz
        col = 6 * node
        dofs[col : col + 3] = displacement
    return dofs


def check_patch_tests(beam, coords):
    print("\nPatch tests")
    passed = True

    rigid_translation = nodal_dofs_from_affine_field(
        coords,
        gradient=np.zeros((3, 3)),
        offset=np.array([0.1, -0.2, 0.3]),
    )
    bending_slope = 1.0e-4
    bending = nodal_dofs_from_affine_field(
        coords,
        gradient=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, bending_slope],
                [0.0, 0.0, 0.0],
            ]
        ),
    )
    for node in range(8):
        bending[6 * node + 3] = -bending_slope
    expected_bending = np.array([0.0, 0.0, 0.0, 0.0, bending_slope, 0.0])

    max_rigid_strain = 0.0
    max_bending_error = 0.0
    for xi in beam.gauss_points:
        for eta in beam.gauss_points:
            for zeta in beam.gauss_points:
                B, _, _, _ = beam.build_B_matrix(xi, eta, zeta, coords)
                max_rigid_strain = max(
                    max_rigid_strain,
                    np.max(np.abs(B @ rigid_translation)),
                )
                max_bending_error = max(
                    max_bending_error,
                    np.max(np.abs((B @ bending) - expected_bending)),
                )

    passed &= report("constant translation strain", max_rigid_strain, 1e-12)
    passed &= report("uy = a*z bending/shear strain", max_bending_error, 1e-10)
    return passed


def main():
    beam = build_single_element()
    coords = beam.get_element_nodes(0)

    checks = [
        check_shape_identities(beam),
        check_jacobian(beam, coords),
        check_patch_tests(beam, coords),
    ]

    print("\nOverall:", "PASS" if all(checks) else "FAIL")


if __name__ == "__main__":
    main()
