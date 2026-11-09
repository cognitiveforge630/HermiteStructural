import unittest
import numpy as np
from core.solvers.beam import build_beam
from core.elements import RecoveredRotationHex8, QuadraticRotationHex8, CubicRotationHex8

class ElementSmokeTests(unittest.TestCase):
    def test_three_explicit_formulations(self):
        for degree, cls in (("linear",RecoveredRotationHex8), ("quadratic",QuadraticRotationHex8), ("cubic",CubicRotationHex8)):
            with self.subTest(degree=degree):
                b = build_beam(nx=3,ny=3,nz=5,formulation=f"{degree}_rotation_hermite")
                self.assertIsInstance(b,cls)
                self.assertEqual(b.ndof_per_node,6)
                ids, t = b.element_transform(0)
                self.assertEqual(ids.shape[0],t.shape[1])
                self.assertEqual(t.shape[0],b.primitive_dofs)
                self.assertTrue(np.all(np.isfinite(t)))
                N = b.displacement_operator(0,(0.,0.,0.))[1]
                self.assertEqual(N.shape,(3,len(ids)))
                self.assertTrue(np.all(np.isfinite(N)))

    def test_default_is_quadratic(self):
        self.assertIsInstance(build_beam(nx=3,ny=3,nz=5),QuadraticRotationHex8)

    def test_bad_formulation_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"Unknown active"):
            build_beam(nx=3,ny=3,nz=5,formulation="unsupported_formulation")

    def test_prescribed_material_invariant(self):
        for degree in ("linear","quadratic","cubic"):
            with self.subTest(degree=degree):
                b=build_beam(nx=3,ny=3,nz=5,formulation=f"{degree}_rotation_hermite")
                self.assertAlmostEqual(b.orthotropic_constants.Ez,1_100_000.)
