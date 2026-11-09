"""Manufactured pointwise equilibrium test for the cubic source stress."""
import unittest
import numpy as np
from core.solvers.beam import build_beam
from core.solvers.strong_residual import point_residual,sample_volume_residual


class StrongResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beam=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')

    def test_affine_extension_has_zero_divergence(self):
        beam=self.beam
        u=np.zeros(len(beam.grid_points)*6)
        u[beam.global_dofs[:,2]]=1e-4*beam.grid_points[:,2]
        res=sample_volume_residual(beam,u,density=0.,gauss_order=2)
        self.assertLess(res['weighted_sample_rms_div_sigma_plus_b_lbf_in3'],1e-8)
        self.assertFalse(res['certified_bound'])

    def test_reject_centroid_only_stress_balance_diagnostic(self):
        beam=self.beam
        u=np.zeros(len(beam.grid_points)*6)
        with self.assertRaisesRegex(ValueError, 'center-only diagnostics prohibited'):
            sample_volume_residual(beam,u,gauss_order=1)
        with self.assertRaisesRegex(ValueError, 'Centroid-based stress diagnostics'):
            point_residual(beam,u,0,parent=(0.,0.,0.))

    def test_exact_quadratic_displacement_has_known_shear_divergence(self):
        beam=self.beam
        alpha=1e-5
        u=np.zeros(len(beam.grid_points)*6)
        u[beam.global_dofs[:,1]]=alpha*beam.grid_points[:,2]**2
        u[beam.global_dofs[:,3]]=-alpha*beam.grid_points[:,2] # half-curl_x
        # gamma_yz=du_y/dz=2 alpha z, hence (div sigma)_y=2G_yz alpha.
        expected=2*beam.D[4,4]*alpha
        r=point_residual(beam,u,4,(.1,-.3,.2),density=0.)
        np.testing.assert_allclose(r,[0,expected,0],rtol=1e-7,atol=1e-7)

    def test_nonzero_weight_without_stress_has_gravity_residual(self):
        beam=self.beam
        u=np.zeros(len(beam.grid_points)*6)
        r=point_residual(beam,u,0,density=.0174)
        np.testing.assert_allclose(r,[0,-.0174,0],atol=1e-12)
        with self.assertRaises(ValueError):
            point_residual(beam,u,0,parent=(1,0,0))

if __name__=='__main__':unittest.main()
