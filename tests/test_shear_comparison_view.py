"""No-GUI verification of signed two-field, one-sided shear visualization."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from core.solvers.beam import build_beam
from core.visualization.sampling import sample_volume
from core.visualization.cubic_stress import analytical_stress_at
from core.visualization.shear_comparison import (
    FE_SHEAR, BEAM_SHEAR, DIFFERENCE, prepare_shear_comparison,
    shared_signed_limits, export_shear_vtu,
)


class TwoPaneShearComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beam = build_beam(nx=3, ny=3, nz=5,
                              formulation='cubic_rotation_hermite')
        cls.U = np.zeros(len(cls.beam.grid_points)*6)
        cls.U[1::6] = 0.001 * cls.beam.grid_points[:, 2]

    def test_raw_fe_and_analytical_are_separate_original_corner_fields(self):
        sample = sample_volume(self.beam, self.U, subdivisions=1)
        data = prepare_shear_comparison(sample, self.beam)
        n = sample.connectivity.shape[0]
        self.assertEqual(data.points.shape, (n*8,3))
        self.assertEqual(data.connectivity.shape, (n,8))
        self.assertEqual(np.unique(data.connectivity).size, n*8)
        np.testing.assert_array_equal(data.connectivity.ravel(), np.arange(n*8))
        np.testing.assert_allclose(data.fields[FE_SHEAR], sample.stress_corner_values[:,:,4].ravel())
        _, independent_beam_tau = analytical_stress_at(sample.stress_corner_points.reshape(-1,3),self.beam)
        np.testing.assert_allclose(data.fields[BEAM_SHEAR], independent_beam_tau)
        np.testing.assert_allclose(data.fields[DIFFERENCE],data.fields[FE_SHEAR]-data.fields[BEAM_SHEAR])
        # The reference is zero at the free top/bottom surfaces; FE is neither
        # forcibly zeroed nor replaced by the beam profile.
        surface = np.isclose(np.abs(data.stress_evaluation_points[:,1]),self.beam.Ly/2)
        self.assertTrue(surface.any())
        np.testing.assert_allclose(data.fields[BEAM_SHEAR][surface],0,atol=1e-12)
        self.assertGreater(np.max(np.abs(data.fields[FE_SHEAR]-data.fields[BEAM_SHEAR])),1e-6)

    def test_inset_samples_are_labeled_and_do_not_move_display_geometry(self):
        sample = sample_volume(self.beam,self.U,subdivisions=1,
            stress_sample_location='inset_corners',stress_inset_fraction=0.95)
        data = prepare_shear_comparison(sample,self.beam)
        self.assertEqual(data.sample_location,'inset_corners')
        np.testing.assert_allclose(data.points,sample.points[sample.connectivity].reshape(-1,3))
        self.assertGreater(np.max(np.abs(data.points-data.stress_evaluation_points)),0)

    def test_shared_signed_color_range_and_finite_data(self):
        data=prepare_shear_comparison(sample_volume(self.beam,self.U,subdivisions=1),self.beam)
        a,b=shared_signed_limits(data)
        self.assertAlmostEqual(a,-b)
        self.assertGreaterEqual(b,np.max(np.abs(data.fields[FE_SHEAR])))
        self.assertGreaterEqual(b,np.max(np.abs(data.fields[BEAM_SHEAR])))

    def test_vtu_has_both_signed_fields_not_averaged(self):
        try:
            import vtk
            from vtk.util.numpy_support import vtk_to_numpy
        except ImportError:
            self.skipTest('VTK unavailable')
        data=prepare_shear_comparison(sample_volume(self.beam,self.U,subdivisions=1),self.beam)
        with tempfile.TemporaryDirectory() as td:
            path=export_shear_vtu(data,Path(td)/'dual_shear.vtu')
            r=vtk.vtkXMLUnstructuredGridReader()
            r.SetFileName(str(path));r.Update()
            grid=r.GetOutput()
            self.assertEqual(grid.GetNumberOfPoints(),len(data.points))
            self.assertEqual(grid.GetNumberOfCells(),len(data.connectivity))
            for key in (FE_SHEAR,BEAM_SHEAR,DIFFERENCE):
                array=grid.GetPointData().GetArray(key)
                self.assertIsNotNone(array,key)
                np.testing.assert_allclose(vtk_to_numpy(array),data.fields[key],rtol=1e-12)
            self.assertEqual(int(vtk_to_numpy(grid.GetFieldData().GetArray('Certified 3D error bound'))[0]),0)

    def test_tutorial_export_is_original_fe_corner_mesh(self):
        """Full solved integration, independent of whether PyVista is installed."""
        import contextlib
        import io
        import json
        import sys
        from unittest.mock import patch
        from tutorials.cantilever import main
        try:
            import vtk
        except ImportError:
            self.skipTest('VTK required for the export-only integration smoke test')
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)/'solved.json'
            argv=['cubic','--nx','3','--ny','3','--nz','5',
                  '--shear-export-only','--out',str(base)]
            with patch.object(sys,'argv',argv),contextlib.redirect_stdout(io.StringIO()):
                main(default_formulation='cubic_rotation_hermite')
            report=json.loads(base.read_text())
            vtupath=Path(report['cubic_shear_comparison_vtu'])
            self.assertTrue(vtupath.exists())
            r=vtk.vtkXMLUnstructuredGridReader()
            r.SetFileName(str(vtupath));r.Update();grid=r.GetOutput()
            self.assertEqual(grid.GetNumberOfCells(),(3-1)*(3-1)*(5-1))
            self.assertEqual(grid.GetNumberOfPoints(),grid.GetNumberOfCells()*8)
            self.assertIsNotNone(grid.GetPointData().GetArray(FE_SHEAR))
            self.assertIsNotNone(grid.GetPointData().GetArray(BEAM_SHEAR))
            self.assertFalse(report['cubic_shear_comparison_stress_sources']['certified_3d_error_bound'])

    def test_invalid_vtu_extension_rejected(self):
        data=prepare_shear_comparison(sample_volume(self.beam,self.U,subdivisions=1),self.beam)
        with self.assertRaises(ValueError):
            export_shear_vtu(data,'somewhere.png')

if __name__=='__main__':
    unittest.main()
