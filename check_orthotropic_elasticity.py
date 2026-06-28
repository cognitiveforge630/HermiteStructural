"""Fast checks for the orthotropic elasticity matrix implementation."""

import numpy as np

from FEMHermiteBeamRegion import FEMHermiteBeamRegion, OrthotropicElasticity
from backups.redwood_orthotropic_study import (
    THREE_POINT_GAUSS_ORDER,
    build_beam,
    redwood_orthotropic_material,
)
from solver import (
    REFERENCE_LONGITUDINAL_E,
    REFERENCE_TARGET_UY,
    build_beam as build_solver_beam,
    calibrated_longitudinal_E,
    cantilever_point_load_deflection,
    I,
    L,
    P,
)


def isotropic_reference_matrix(E, nu):
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    return np.array(
        [
            [lam + 2.0 * mu, lam, lam, 0.0, 0.0, 0.0],
            [lam, lam + 2.0 * mu, lam, 0.0, 0.0, 0.0],
            [lam, lam, lam + 2.0 * mu, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, mu, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, mu, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, mu],
        ]
    )


def check_isotropic_equivalence():
    E = 2_000_000.0
    nu = 0.30
    G = E / (2.0 * (1.0 + nu))
    material = OrthotropicElasticity(
        Ex=E,
        Ey=E,
        Ez=E,
        nuxy=nu,
        nuxz=nu,
        nuyz=nu,
        Gxy=G,
        Gyz=G,
        Gxz=G,
    )
    actual = material.elastic_matrix()
    expected = isotropic_reference_matrix(E, nu)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-7)


def check_redwood_matrix():
    material = redwood_orthotropic_material(1_100_000.0)
    D = FEMHermiteBeamRegion.orthotropic_elastic_matrix(material)
    np.testing.assert_allclose(D, D.T, rtol=1e-12, atol=1e-9)
    eigenvalues = np.linalg.eigvalsh(D)
    if np.any(eigenvalues <= 0.0):
        raise AssertionError(f"Redwood orthotropic D is not positive definite: {eigenvalues}")


def check_three_point_gauss_defaults():
    beam = FEMHermiteBeamRegion(
        Lx=4.0,
        Ly=8.0,
        Lz=8.0,
        nx=2,
        ny=2,
        nz=2,
        E=2_000_000.0,
        nu=0.30,
        gamma=0.0,
        w=0.0,
    )
    if beam.gauss_order != THREE_POINT_GAUSS_ORDER:
        raise AssertionError(f"Default Gauss order should be 3, got {beam.gauss_order}")
    if len(beam.gauss_points) != THREE_POINT_GAUSS_ORDER:
        raise AssertionError("Expected exactly three Gauss points per direction")

    material = redwood_orthotropic_material(1_100_000.0)
    redwood_beam = build_beam(
        mesh=(2, 2, 2),
        gauss_order=THREE_POINT_GAUSS_ORDER,
        material=material,
        E_longitudinal=1_100_000.0,
    )
    if redwood_beam.gauss_order != THREE_POINT_GAUSS_ORDER:
        raise AssertionError("Redwood/AWC beam should use 3-point Gauss integration")


def check_solver_uses_orthotropic_three_point_beam():
    beam = build_solver_beam()
    if beam.gauss_order != THREE_POINT_GAUSS_ORDER:
        raise AssertionError(f"solver.py should use 3-point Gauss, got {beam.gauss_order}")
    if beam.orthotropic_constants is None:
        raise AssertionError("solver.py should build an orthotropic material, not isotropic D")
    if not isinstance(beam.orthotropic_constants, OrthotropicElasticity):
        raise AssertionError("solver.py should pass an OrthotropicElasticity instance")
    if abs(beam.orthotropic_constants.Ez - REFERENCE_LONGITUDINAL_E) > 1.0e-6:
        raise AssertionError("solver.py seed Ez should start at the AWC/reference E")
    D = beam.D
    np.testing.assert_allclose(D, D.T, rtol=1e-12, atol=1e-9)
    if np.any(np.linalg.eigvalsh(D) <= 0.0):
        raise AssertionError("solver.py orthotropic D is not positive definite")


def check_calibration_formula_keeps_reference_target_separate():
    seed_uy = 0.75 * REFERENCE_TARGET_UY
    calibrated_E = calibrated_longitudinal_E(REFERENCE_LONGITUDINAL_E, seed_uy)
    expected_E = 0.75 * REFERENCE_LONGITUDINAL_E
    if abs(calibrated_E - expected_E) > 1.0e-9:
        raise AssertionError("calibrated_longitudinal_E did not apply linear FEA scaling")

    calibrated_theory = cantilever_point_load_deflection(P, L, calibrated_E, I)
    if abs(calibrated_theory - REFERENCE_TARGET_UY) < 1.0e-6:
        raise AssertionError("calibrated solver E should not overwrite the AWC/reference target E")


def main():
    check_isotropic_equivalence()
    check_redwood_matrix()
    check_three_point_gauss_defaults()
    check_solver_uses_orthotropic_three_point_beam()
    check_calibration_formula_keeps_reference_target_separate()
    print("Orthotropic elasticity checks passed.")


if __name__ == "__main__":
    main()
