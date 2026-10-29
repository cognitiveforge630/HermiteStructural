import unittest
import numpy as np
from core.solvers.beam import build_beam
from core.solvers.loads import gravity_load, build_end_face_traction_load_vector
from core.reference.awc import awc_reference_components
from core.reference.beam_theory import P

class LoadTests(unittest.TestCase):
    def test_theory_contains_gravity(self):
        components=awc_reference_components()
        self.assertAlmostEqual(components["total_weight_lbf"],66.816)
        self.assertAlmostEqual(components["total_uy_in"],components["end_load_uy_in"]+components["self_weight_uy_in"])
        self.assertLess(components["total_uy_in"],0.)

    def test_integrated_resultants_for_each_recovery(self):
        for degree in ("linear","quadratic","cubic"):
            with self.subTest(degree=degree):
                beam=build_beam(nx=3,ny=3,nz=5,formulation=f"{degree}_rotation_hermite")
                ft,_,_=build_end_face_traction_load_vector(beam)
                fg=gravity_load(beam,0.0174)
                idx=beam.global_dofs[:,1]
                self.assertAlmostEqual(float(ft[idx].sum()),P,places=7)
                self.assertAlmostEqual(float(fg[idx].sum()),-66.816,places=6)
                self.assertTrue(np.all(np.isfinite(ft+fg)))
