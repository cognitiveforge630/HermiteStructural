"""Fast verification tests only: no production-size study or GUI."""
import unittest
import numpy as np
from tools.p_order_study import (
    build_plan, tensor_shapes, patch_checks, element_verification, local_nodes,
    baseline_compare, solve_study,
)

class TestPolynomialNodalHex(unittest.TestCase):
    def test_plan_sizes(self):
        q=build_plan(2,3,5,25)
        c=build_plan(3,3,5,25)
        self.assertEqual(q['interpolation_nodes_per_element'],27)
        self.assertEqual(c['interpolation_nodes_per_element'],64)
        self.assertEqual(q['free_and_fixed_translation_DOFs'],6615)
        self.assertEqual(c['free_and_fixed_translation_DOFs'],19929)
        self.assertEqual(q['physical_element_divisions'],c['physical_element_divisions'])

    def test_tensor_patch_and_rigid_modes(self):
        for p in (2,3):
            with self.subTest(order=p):
                self.assertTrue(patch_checks(p)['passed'])
                self.assertTrue(element_verification(p)['passed'])

    def test_shared_high_order_nodes(self):
        for p in (2,3):
            with self.subTest(order=p):
                dims=(2*p+1,p+1,2*p+1)
                left,_=local_nodes(p,0,0,0,dims)
                right,_=local_nodes(p,1,0,0,dims)
                a=left.reshape((p+1,p+1,p+1))[-1,:,:]
                b=right.reshape((p+1,p+1,p+1))[0,:,:]
                np.testing.assert_array_equal(a,b)

    def test_no_reduced_quadrature(self):
        with self.assertRaises(ValueError):
            build_plan(3,3,5,25,3)

    def test_small_solves_reactions_and_energy(self):
        # Only 16 cells per variant, well below the intended 192-cell comparison.
        for p in (2,3):
            with self.subTest(order=p):
                d=solve_study(p,3,3,5)
                self.assertLess(d['equilibrium']['force_defect_inf_lbf'],1e-5)
                self.assertLess(d['equilibrium']['moment_defect_inf_lbf_in'],1e-4)
                self.assertLess(d['equilibrium']['Clapeyron_work_energy_relative_defect'],1e-8)
                self.assertLess(d['equilibrium']['analytical_resultant_force_Fy_error_lbf'],1e-8)
                self.assertEqual(len(d['source_sections']),8)
                self.assertFalse(d['certified_error_bound'])

if __name__=='__main__':
    unittest.main()
