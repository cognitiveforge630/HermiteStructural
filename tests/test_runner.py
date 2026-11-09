import unittest
import json
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch
from tutorials.cantilever import main

class RunnerTests(unittest.TestCase):
    def test_cli_rejects_bad_cubic_mesh(self):
        with patch("sys.argv",["cubic","--nz","3"]):
            with self.assertRaises(SystemExit) as exc:
                main(default_formulation="cubic_rotation_hermite")
            self.assertEqual(exc.exception.code,2)

    def test_cli_rejects_nonpositive_gravity(self):
        with patch("sys.argv",["quadratic","--weight-density","0"]):
            with self.assertRaises(SystemExit) as exc:
                main()
            self.assertEqual(exc.exception.code,2)

    def test_section_only_outputs_research_report_without_sample_screen(self):
        with TemporaryDirectory() as tmp:
            output=Path(tmp)/'cubic_3x3x5.json'
            with patch('sys.argv',['cubic','--nx','3','--ny','3',
                     '--nz','5','--section-only','--section-axial-positions','0','--out',str(output)]):
                result=main(default_formulation='cubic_rotation_hermite')
            report=json.loads(Path(result['cubic_source_section_file']).read_text())
            self.assertEqual(report['cross_section_quadrature_order'],4)
            self.assertFalse(report['certified_error_bound'])
            self.assertFalse(result['sampled_cubic_stress_screen_performed'])
            self.assertFalse(result['stress_distribution_certified'])
            self.assertAlmostEqual(report['max_midspan_shear_y_difference_pct'],18.8384380029,places=4)

    def test_optional_shear_balance_writes_unmodified_raw_diagnostic(self):
        with TemporaryDirectory() as tmp:
            output=Path(tmp)/'shear.json'
            with patch('sys.argv',['cubic','--nx','3','--ny','3','--nz','5',
                    '--section-only','--shear-balance','--out',str(output)]):
                result=main(default_formulation='cubic_rotation_hermite')
            path=Path(result['cubic_axial_shear_balance_file'])
            report=json.loads(path.read_text())
            self.assertFalse(report['certified_error_bound'])
            self.assertTrue(report['does_not_modify_solver_or_stress'])
            self.assertEqual(len(report['cells']),4)
            self.assertGreater(result['cubic_max_interior_sampled_shear_difference_pct'],0)
            self.assertIn('cubic_max_interior_axially_averaged_shear_difference_pct',result)

    def test_section_only_rejects_visual_screening_flags(self):
        with patch('sys.argv',['cubic','--section-only','--view']):
            with self.assertRaises(SystemExit) as exc:
                main(default_formulation='cubic_rotation_hermite')
            self.assertEqual(exc.exception.code,2)

    def test_section_positions_reject_interfaces_nonfinite_and_duplicates(self):
        for positions in (['-1'], ['nan'], ['0', '0']):
            with self.subTest(positions=positions), patch('sys.argv',[
                    'cubic','--section-axial-positions',*positions]):
                with self.assertRaises(SystemExit) as exc:
                    main(default_formulation='cubic_rotation_hermite')
                self.assertEqual(exc.exception.code,2)

    def test_multiple_section_positions_preserve_midpoint_result(self):
        with TemporaryDirectory() as tmp:
            output=Path(tmp)/'cubic_slices.json'
            with patch('sys.argv',['cubic','--nx','3','--ny','3',
                     '--nz','5','--section-only','--section-axial-positions',
                     '-0.5','0','0.5','--out',str(output)]):
                result=main(default_formulation='cubic_rotation_hermite')
            report=json.loads(Path(result['cubic_source_section_file']).read_text())
            self.assertEqual(report['positions_parent_zeta'],[-0.5,0.,0.5])
            self.assertEqual(len(report['sections']),12)
            midpoints=[s for s in report['sections'] if s['axial_parent_zeta']==0
                       and 30<=s['z_in']<=90]
            self.assertAlmostEqual(max(s['shear_y_reference_difference_pct']
                                     for s in midpoints),18.8384380029,places=4)
