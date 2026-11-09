"""Internal-face equilibrium jump regression checks."""
import unittest
import numpy as np
from core.solvers.beam import build_beam
from core.solvers.traction_jumps import source_traction_jumps


class TractionJumpsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beam=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')

    def test_affine_extension_has_no_internal_jump(self):
        beam=self.beam
        u=np.zeros(len(beam.grid_points)*6)
        u[beam.global_dofs[:,2]]=1e-4*beam.grid_points[:,2]
        result=source_traction_jumps(beam,u,order=3)
        for axis,value in result['per_face_normal_axis'].items():
            self.assertGreater(value['faces'],0)
            self.assertLess(value['maximum_face_traction_jump_rms_psi'],1e-7)
        self.assertFalse(result['certified_bound'])

    def test_random_dofs_also_have_continuous_stress_traction(self):
        # The recovered tensor Hermite field is numerically C1 on the
        # structured shared-vertex meshes. Zero jumps do NOT imply div sigma+b=0.
        beam=self.beam
        rng=np.random.default_rng(7)
        u=rng.normal(scale=1e-4,size=len(beam.grid_points)*6)
        jumps=source_traction_jumps(beam,u,order=3)
        for value in jumps['per_face_normal_axis'].values():
            self.assertLess(value['area_weighted_traction_jump_rms_psi'],1e-8)

    def test_reject_bad_input(self):
        with self.assertRaises(ValueError):source_traction_jumps(self.beam,np.zeros(1))
        with self.assertRaises(ValueError):source_traction_jumps(self.beam,np.zeros(len(self.beam.grid_points)*6),order=1)

if __name__=='__main__':unittest.main()
