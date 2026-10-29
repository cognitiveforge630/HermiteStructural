"""Protect transverse cubic-warping completeness from regression."""
import unittest
from core.solvers.beam import build_beam
from core.verification.shear_warping import audit_cubic_y_warping


class CubicTransverseWarpingTests(unittest.TestCase):
    def test_three_depth_nodes_explicitly_mark_incomplete(self):
        beam = build_beam(nx=3, ny=3, nz=9, formulation='cubic_rotation_hermite')
        report = audit_cubic_y_warping(beam)
        self.assertFalse(report['y_cubic_recovery_term_present'])
        self.assertGreater(report['relative_discrete_strain_L2'], 0.1)
        self.assertFalse(report['certified_error_bound'])

    def test_four_depth_nodes_reproduce_cubic_warping_exactly(self):
        beam = build_beam(nx=3, ny=4, nz=9, formulation='cubic_rotation_hermite')
        report = audit_cubic_y_warping(beam)
        self.assertTrue(report['y_cubic_recovery_term_present'])
        self.assertLess(report['relative_discrete_strain_L2'], 1e-11)
        self.assertLess(report['maximum_absolute_strain_difference'], 1e-10)
        self.assertTrue(report['samples_include_source_corners'])
        self.assertTrue(report['actual_source_B_and_T'])

    def test_quadratic_not_misreported_as_cubic(self):
        beam = build_beam(nx=3, ny=4, nz=9, formulation='quadratic_rotation_hermite')
        with self.assertRaises(TypeError):
            audit_cubic_y_warping(beam)

    def test_rejects_invalid_sample_positions(self):
        beam = build_beam(nx=3, ny=4, nz=9, formulation='cubic_rotation_hermite')
        with self.assertRaises(ValueError):
            audit_cubic_y_warping(beam, points=[(0, 0, 1.01)])

if __name__ == '__main__':
    unittest.main()
