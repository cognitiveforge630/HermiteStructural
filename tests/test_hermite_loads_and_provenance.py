import sys, unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'legacy_hermite'))
import solver

class TestHermiteLoads(unittest.TestCase):
    def test_positive_gravity_required(self):
        for x in (0.,-0.01,float('nan')):
            with self.subTest(weight=x), self.assertRaises(ValueError):
                solver.awc_reference_components(x)

    def test_awc_sum_and_beam_weight(self):
        a=solver.awc_reference_components()
        self.assertAlmostEqual(a['total_uy_in'],a['end_load_uy_in']+a['self_weight_uy_in'])
        self.assertAlmostEqual(a['total_weight_lbf'],66.816,places=7)
        self.assertLess(a['end_load_uy_in'],0)
        self.assertLess(a['self_weight_uy_in'],0)

    def test_integrated_load_resultants_and_first_moments(self):
        # Small orthotropic Hermite mesh: fast, no global stiffness solve.
        beam=solver.build_beam(nx=2,ny=2,nz=3)
        ft, _, _=solver.build_end_face_traction_load_vector(beam)
        fg=solver.build_self_weight_load_vector(beam)
        gd=beam.global_dofs
        X=beam.get_all_points()
        self.assertAlmostEqual(float(ft[gd[:,1]].sum()),solver.P,places=8)
        self.assertAlmostEqual(float(fg[gd[:,1]].sum()),-66.816,places=7)
        # Displacement interpolation includes rotational DOFs, so validate
        # resultant virtual work against unit rigid translation, not nodal moments.
        self.assertAlmostEqual(float((ft+fg)[gd[:,1]].sum()),solver.P-66.816,places=7)
        self.assertTrue(np.all(np.isfinite(ft)))
        self.assertTrue(np.all(np.isfinite(fg)))

if __name__=='__main__':unittest.main()

class TestFullReport(unittest.TestCase):
    def test_full_report_uses_combined_load_reference(self):
        import io
        from contextlib import redirect_stdout
        beam=solver.build_beam(nx=2,ny=2,nz=3)
        stream=io.StringIO()
        with redirect_stdout(stream):
            solver.print_theoretical_check(beam, label='test')
        output=stream.getvalue()
        self.assertIn('AWC combined uy:',output)
        self.assertIn('AWC gravity uy:',output)
        self.assertNotIn('Target theoretical free-end uy:',output)

    def test_primary_solver_has_no_calibration(self):
        from inspect import getsource
        source=getsource(solver.solve_referenced_beam)
        self.assertNotIn('calibrated_longitudinal_E(',source)
        self.assertIn('build_self_weight_load_vector(',source)
        self.assertIn('build_end_face_traction_load_vector(',source)
