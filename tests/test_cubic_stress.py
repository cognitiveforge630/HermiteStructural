"""Cubic stress safety and bending profile research regressions (not certification)."""
import unittest
import numpy as np
from core.materials.design_limits import REDWOOD_BASIS, reference_load_screen
from core.solvers.beam import build_beam
from core.solvers.loads import gravity_load,build_end_face_traction_load_vector
from core.solvers.constraints import solve_with_face_dirichlet
from core.visualization.sampling import sample_volume
from core.visualization.cubic_stress import analytical_stress_at, cubic_stress_audit


class CubicStressTests(unittest.TestCase):
    def test_grade_values_and_reasonable_baseline(self):
        self.assertEqual(REDWOOD_BASIS.E_longitudinal_psi,1100000.)
        result=reference_load_screen(tip_force_lbf=-12.8,density_lbf_in3=.0174)
        self.assertEqual(result['status'],'PASS_PRELIMINARY')
        self.assertLess(result['checks']['bending']['demand'],200)
        self.assertLess(result['checks']['example_tip_serviceability']['utilization'],1)

    def test_rejects_unsafe_extreme_load(self):
        self.assertEqual(reference_load_screen(tip_force_lbf=-300.,density_lbf_in3=.0174)['status'],
                         'FAIL_PRELIMINARY')
        with self.assertRaises(ValueError):
            reference_load_screen(tip_force_lbf=3,density_lbf_in3=.0174)

    def test_beam_reference_is_distinct_and_signed(self):
        beam=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')
        sigma,tau=analytical_stress_at(np.array([[0,4,0],[0,-4,0],[0,0,60]]),beam)
        self.assertGreater(sigma[0],0)
        self.assertLess(sigma[1],0)
        self.assertAlmostEqual(sigma[0],-sigma[1])
        self.assertAlmostEqual(sigma[2],0)
        self.assertLess(tau[2],0)

    def test_cubic_stress_contours_against_beam_and_screen(self):
        b=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')
        f=build_end_face_traction_load_vector(b)[0]+gravity_load(b,.0174)
        U,_=solve_with_face_dirichlet(b,f,lambda x:(0.,0.,0.),lambda x:(0.,0.,0.))
        sample=sample_volume(b,U,subdivisions=2)
        before=sample.cell_fields['Sigma zz (psi)'].copy()
        audit=cubic_stress_audit(sample,b)
        np.testing.assert_array_equal(before,sample.cell_fields['Sigma zz (psi)'])
        self.assertTrue(np.max(abs(before-sample.cell_fields['Beam reference sigma zz (psi)']))>0)
        self.assertTrue(np.all(np.isfinite(sample.cell_fields['Cubic minus beam sigma zz (psi)'])))
        self.assertEqual(audit['status'],'PASS_PRELIMINARY_SAMPLE_SCREEN')
        self.assertFalse(audit['stress_distribution_certified'])
        self.assertEqual(audit['sample_subdivisions'],2)
        # This is only a baseline regression on bending; not rigorous 3D validation.
        self.assertLess(audit['quarter_to_threequarter_span_sigma_zz_relative_L2_pct'],15.)

if __name__=='__main__':
    unittest.main()
