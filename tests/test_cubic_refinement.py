"""Offline refinement comparison tests; do not trigger an FE solve."""
import copy
from pathlib import Path
import unittest

from tools.cubic_refinement import read_case, compare

ROOT = Path(__file__).resolve().parents[1]


class CubicRefinementTests(unittest.TestCase):
    def test_shipped_cases_are_explicitly_sampling_confounded(self):
        coarse = read_case(ROOT/'samples/cubic_3x3x5.json')
        fine = read_case(ROOT/'samples/cubic_3x3x9.json')
        report = compare([coarse, fine])
        self.assertEqual(report['status'], 'RESEARCH_MESH_COMPARISON_NOT_CERTIFIED')
        self.assertFalse(report['certified_bound'])
        self.assertEqual([coarse['subdivisions'],fine['subdivisions']],[2,3])
        self.assertFalse(report['pairs'][0]['matched_sampling_subdivisions'])
        self.assertFalse(report['pairs'][0]['uniform_h_refinement'])
        self.assertAlmostEqual(fine['tip_reference_discrepancy_pct'],0.4651239155,places=7)
        self.assertAlmostEqual(fine['stress_reference_discrepancy_pct'],2.1625164051,places=7)

    def test_uniform_refinement_detected_from_element_lengths(self):
        a = read_case(ROOT/'samples/cubic_3x3x9.json')
        b = copy.deepcopy(a)
        b['mesh'] = [5,5,17]
        b['tip_uy_in'] = a['tip_uy_in'] + 0.0001
        b['result_file'] = 'synthetic_demo_finer_mesh.json'
        report = compare([a,b])
        self.assertTrue(report['pairs'][0]['uniform_h_refinement'])
        self.assertTrue(report['pairs'][0]['matched_sampling_subdivisions'])
        self.assertAlmostEqual(report['pairs'][0]['tip_change_pct_of_fine'],
                               100*abs(b['tip_uy_in']-a['tip_uy_in'])/abs(b['tip_uy_in']))

    def test_refuses_incompatible_or_reversed_cases(self):
        a = read_case(ROOT/'samples/cubic_3x3x9.json')
        b = copy.deepcopy(a)
        b['mesh'] = [5,5,17]
        b['end_force_lbf'] *= 2
        with self.assertRaisesRegex(ValueError,'end_force_lbf'):
            compare([a,b])
        with self.assertRaisesRegex(ValueError,'nondecreasing'):
            compare([a, read_case(ROOT/'samples/cubic_3x3x5.json')])


if __name__ == '__main__':
    unittest.main()
