"""Fast regression checks for restartable and non-centroid cubic stress export."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from core.solvers.beam import build_beam
from core.solvers.solution_archive import save_cubic_solution, load_cubic_solution
from core.visualization.corner_fields import prepare_corner_fields
from core.visualization.sampling import sample_volume
from core.visualization.shear_comparison import FE_SHEAR, BEAM_SHEAR
from tools.replay_cubic_fields import replay, summarize_stress


class CubicReplayFieldTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beam = build_beam(nx=3, ny=3, nz=5, formulation='cubic_rotation_hermite')
        cls.u = np.zeros(len(cls.beam.grid_points) * cls.beam.ndof_per_node)
        cls.u[1::6] = cls.beam.grid_points[:, 2] * 2.e-5

    def test_solution_archive_exact_round_trip_and_reject_summary(self):
        with tempfile.TemporaryDirectory() as td:
            state = Path(td) / 'solved.npz'
            save_cubic_solution(state, self.beam, self.u, end_force=-12.8,
                                weight_density=.0174, fem_tip_uy_in=-.005)
            beam, u, meta = load_cubic_solution(state)
            self.assertEqual((beam.nx, beam.ny, beam.nz), (3, 3, 5))
            np.testing.assert_array_equal(u, self.u)
            self.assertEqual(meta['weight_density_lbf_in3'], .0174)
            self.assertFalse(meta['stress_data_saved'])
            invalid = Path(td) / 'summary.npz'
            np.savez_compressed(invalid, tip=np.array([-.005]))
            with self.assertRaisesRegex(ValueError, 'Expected U and metadata'):
                load_cubic_solution(invalid)

    def test_saved_solution_rejects_changed_transverse_or_shear_modulus(self):
        """A matching longitudinal E is insufficient for stress replay provenance."""
        with tempfile.TemporaryDirectory() as td:
            original = Path(td) / 'original.npz'
            altered = Path(td) / 'altered.npz'
            save_cubic_solution(original, self.beam, self.u, end_force=-12.8,
                                weight_density=.0174)
            with np.load(original, allow_pickle=False) as state:
                u = np.array(state['U'], copy=True)
                meta = json.loads(str(state['metadata'].item()))
            self.assertEqual(np.asarray(meta['constitutive_matrix_psi']).shape, (6, 6))
            # Simulate a material source change: the replay must not silently
            # accept it even though the saved Ez and geometry are unchanged.
            meta['constitutive_matrix_psi'][3][3] *= 1.01
            np.savez_compressed(altered, U=u, metadata=np.array(json.dumps(meta)))
            with self.assertRaisesRegex(ValueError, 'constitutive matrix'):
                load_cubic_solution(altered)
            del meta['constitutive_matrix_psi']
            np.savez_compressed(altered, U=u, metadata=np.array(json.dumps(meta)))
            with self.assertRaisesRegex(ValueError, 'constitutive matrix'):
                load_cubic_solution(altered)

    def test_corner_stress_is_direct_source_B_not_cell_mean(self):
        sample = sample_volume(self.beam, self.u, subdivisions=1,
                               stress_sample_location='inset_corners', stress_inset_fraction=.95)
        data = prepare_corner_fields(sample, self.beam, end_force=-12.8, density=.0174)
        np.testing.assert_allclose(data.fields['Tau yz (psi)'],
                                   sample.stress_corner_values[:, :, 4].ravel())
        self.assertEqual(data.sample_location, 'inset_corners')
        self.assertEqual(len(data.points), 8 * len(sample.connectivity))
        self.assertTrue(np.any(np.abs(data.points - data.stress_evaluation_points) > 0))
        vm = data.fields['Von Mises (psi)']
        self.assertTrue(np.isfinite(vm).all())
        self.assertTrue((vm >= 0).all())
        # No inter-element FE stress interpolation or centroid substitution.
        self.assertEqual(data.fields['Sigma zz (psi)'].shape,
                         (len(sample.connectivity)*8,))
        report = summarize_stress(data, self.beam, {
            'weight_density_lbf_in3': .0174, 'end_force_lbf': -12.8,
            'node_mesh': [3,3,5],
        }, subdivisions=1, inset_fraction=.95)
        self.assertAlmostEqual(report['beam_maximum_shear_reference_psi'], 3.732, places=10)
        self.assertAlmostEqual(report['sampled_beam_tau_yz_maximum']['absolute_psi'],
                               max(abs(data.fields[BEAM_SHEAR])), places=12)
        self.assertAlmostEqual(report['sampled_fem_tau_yz_maximum']['absolute_psi'],
                               max(abs(data.fields[FE_SHEAR])), places=12)
        self.assertFalse(report['within_one_percent_claim'])

    def test_section_only_runner_can_save_replayable_dofs(self):
        import contextlib
        import io
        from unittest.mock import patch
        from tutorials.cantilever import main
        with tempfile.TemporaryDirectory() as td:
            result_file = Path(td) / 'result.json'
            state = Path(td) / 'complete_state.npz'
            args = ['cubic', '--nx', '3', '--ny', '3', '--nz', '5',
                    '--section-only', '--section-axial-positions', '0',
                    '--save-solution', str(state), '--out', str(result_file)]
            with patch('sys.argv', args), contextlib.redirect_stdout(io.StringIO()):
                result = main(default_formulation='cubic_rotation_hermite')
            self.assertTrue(state.is_file())
            self.assertEqual(result['solution_state_file'], str(state))
            beam, u, meta = load_cubic_solution(state)
            from core.solvers.results import free_end_uy
            self.assertAlmostEqual(free_end_uy(beam, u), result['fem_tip_uy_in'], places=10)
            self.assertAlmostEqual(meta['fem_tip_uy_in'], result['fem_tip_uy_in'], places=10)
            self.assertAlmostEqual(result['total_y_force_lbf'], -79.616, places=10)

    def test_replay_exports_all_point_corner_fields_vtu(self):
        try:
            import vtk
            from vtk.util.numpy_support import vtk_to_numpy
        except ImportError:
            self.skipTest('VTK unavailable')
        with tempfile.TemporaryDirectory() as td:
            state = Path(td) / 'state.npz'
            save_cubic_solution(state, self.beam, self.u, end_force=-12.8,
                                weight_density=.0174)
            prefix = Path(td) / 'replayed'
            report = replay(state, prefix, subdivisions=1, inset=.95)
            self.assertEqual(report['corner_samples'], 128)
            self.assertEqual(report['stress_sample_location'], 'inset_corners')
            self.assertFalse(report['certified_3d_stress_error_bound'])
            self.assertEqual(json.loads((Path(td)/'replayed_stress_comparison.json').read_text())['corner_samples'], 128)
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(str(Path(td)/'replayed_one_sided_corner_fields.vtu'))
            reader.Update()
            grid = reader.GetOutput()
            self.assertEqual(grid.GetNumberOfCells(), 16)
            self.assertEqual(grid.GetNumberOfPoints(), 128)
            self.assertIsNotNone(grid.GetPointData().GetArray('Displacement (in)'))
            for field in (BEAM_SHEAR, FE_SHEAR, 'Von Mises (psi)',
                          'Sigma zz (psi)', 'Longitudinal tension (psi)',
                          'Longitudinal compression (psi)',
                          'Stress evaluation Z (in)'):
                with self.subTest(field=field):
                    values = vtk_to_numpy(grid.GetPointData().GetArray(field))
                    self.assertEqual(len(values), 128)
                    self.assertTrue(np.isfinite(values).all())
            self.assertIsNone(grid.GetCellData().GetArray('Von Mises (psi)'))


if __name__ == '__main__':
    unittest.main()
