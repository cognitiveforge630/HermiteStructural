"""Research convergence comparator, including fail-closed input validation."""
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from tools.section_refinement import read_case,compare


class SourceSectionComparatorTests(unittest.TestCase):
    def _fixture(self, folder, *, mesh=(3,3,5), tip=-.115, order=4, clamp='fixed', weight=.5568):
        folder=Path(folder)
        name='_'.join(map(str,mesh))
        data={
            'formulation':'cubic_rotation_hermite','node_mesh':list(mesh),
            'cubic_source_section_file':f'{name}_source_sections.json',
            'end_force_lbf':-12.8,'awc_reference':{'total_uy_in':-.116,'weight_per_length_lbf_in':weight},
            'fem_tip_uy_in':tip,'longitudinal_E_psi':1100000,
            'boundary_condition':clamp,'sampled_cubic_stress_screen_performed':False,
            'linear_solver_diagnostics':{'linear_solve':{'relative_backward_error':1e-13}},
        }
        report={'certified_error_bound':False,'status':'RESEARCH_SECTIONAL_DIAGNOSTIC_NOT_A_POSTERIORI',
                'visualization_subdivisions_affect_results':False,
                'cross_section_quadrature_order':order,'positions_parent_zeta':[0.],
                'sections':[{'z_in':0}]*(mesh[2]-1),
                'max_midspan_shear_y_difference_pct':7.,
                'max_midspan_moment_x_difference_pct':2.}
        (folder/f'{name}.json').write_text(json.dumps(data))
        (folder/f'{name}_source_sections.json').write_text(json.dumps(report))
        return folder/f'{name}.json'

    def test_compatible_cases_and_non_uniform_warning(self):
        with TemporaryDirectory() as d:
            first=read_case(self._fixture(d))
            second=read_case(self._fixture(d,mesh=(3,3,9),tip=-.116))
            comp=compare([first,second]);self.assertFalse(comp['certified_error_bound'])
            self.assertFalse(comp['pairs'][0]['uniform_refinement'])
            self.assertTrue(comp['pairs'][0]['tip_beam_discrepancy_decreased'])

    def test_rejects_quadrature_and_clamp_changes(self):
        with TemporaryDirectory() as d:
            first=read_case(self._fixture(d))
            changed=read_case(self._fixture(d,mesh=(3,3,9),order=5))
            with self.assertRaises(ValueError):compare([first,changed])
            changed=read_case(self._fixture(d,mesh=(3,3,9),clamp='free'))
            with self.assertRaises(ValueError):compare([first,changed])
            changed=read_case(self._fixture(d,mesh=(3,3,9),weight=1.))
            with self.assertRaises(ValueError):compare([first,changed])

if __name__=='__main__':unittest.main()
