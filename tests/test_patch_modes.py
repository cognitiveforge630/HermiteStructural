"""Affine/rigid-body strains must be exact for each active recovery model."""
import unittest
import numpy as np
from core.solvers.beam import build_beam


class PatchModes(unittest.TestCase):
    def test_rigid_rotation_and_affine_deformation(self):
        rigid = np.array([[0., -.002, .001], [.002, 0., .003], [-.001, -.003, 0.]])
        affine = np.array([[.002, .001, -.002], [.004, -.003, .004], [.001, .006, .007]])
        for degree in ("linear", "quadratic", "cubic"):
            beam = build_beam(nx=3, ny=3, nz=5, formulation=f"{degree}_rotation_hermite")
            for name, A in (("rigid", rigid), ("affine", affine)):
                with self.subTest(degree=degree, field=name):
                    nodal = np.zeros((len(beam.grid_points), 6))
                    nodal[:, :3] = beam.grid_points @ A.T
                    nodal[:, 3:] = (
                        (A[2, 1] - A[1, 2]) / 2,
                        (A[0, 2] - A[2, 0]) / 2,
                        (A[1, 0] - A[0, 1]) / 2,
                    )
                    exact = np.array([
                        A[0, 0], A[1, 1], A[2, 2],
                        A[0, 1] + A[1, 0], A[1, 2] + A[2, 1], A[0, 2] + A[2, 0],
                    ])
                    for cell in (0, 3, 12, 15):
                        for point in ((0., 0., 0.), (-.66, .25, .53)):
                            strain = beam.strain_at(cell, nodal.ravel(), point)
                            np.testing.assert_allclose(strain, exact, rtol=0, atol=1e-11)
