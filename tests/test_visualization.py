"""Verify visualization fields without requiring a GPU, desktop, or PyVista."""
import unittest
import importlib.util
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from core.solvers.beam import build_beam
from core.visualization.sampling import sample_volume, POINT_FIELDS, CELL_FIELDS


class VisualizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beams = {
            name: build_beam(nx=3,ny=3,nz=5,formulation=f"{name}_rotation_hermite")
            for name in ("linear","quadratic","cubic")
        }

    def test_sampled_mesh_and_linear_field(self):
        """All three recover exact affine displacement from identical DOFs."""
        for name,b in self.beams.items():
            with self.subTest(formulation=name):
                U=np.zeros((len(b.grid_points),6))
                U[:,1]=-0.002*b.grid_points[:,2]
                sample=sample_volume(b,U.ravel(),subdivisions=2)
                self.assertEqual(sample.connectivity.shape, (b.get_num_elements()*8,8))
                self.assertEqual(sample.points.shape[0],5*5*9)
                self.assertTrue(np.all((sample.connectivity>=0)&
                                       (sample.connectivity<len(sample.points))))
                np.testing.assert_allclose(sample.displacement[:,1],
                                           -0.002*sample.points[:,2],atol=1e-11)
                np.testing.assert_allclose(sample.point_fields["Displacement magnitude (in)"],
                                           np.linalg.norm(sample.displacement,axis=1))
                self.assertEqual(set(sample.cell_fields),set(CELL_FIELDS))
                self.assertEqual(set(sample.point_fields),set(POINT_FIELDS))
                self.assertTrue(np.all(sample.cell_fields["Von Mises (psi)"]>=0))

    def test_stresses_are_actual_element_strain(self):
        b=self.beams["quadratic"]
        U=np.zeros(len(b.grid_points)*6)
        U[1::6] = -0.002*b.grid_points[:,2]
        sampled=sample_volume(b,U,subdivisions=1)
        reference=np.array([b.elastic_matrix()@b.strain_at(0,U,p)
             for p in ((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
                       (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))])
        np.testing.assert_allclose(sampled.stress_corner_values[0],reference,rtol=1e-12,atol=1e-10)
        np.testing.assert_allclose(
            [sampled.cell_fields["Sigma xx (psi)"][0],
             sampled.cell_fields["Sigma yy (psi)"][0],
             sampled.cell_fields["Sigma zz (psi)"][0],
             sampled.cell_fields["Tau xy (psi)"][0],
             sampled.cell_fields["Tau yz (psi)"][0],
             sampled.cell_fields["Tau xz (psi)"][0]],reference.mean(axis=0),rtol=1e-12,atol=1e-10)

    @unittest.skipUnless(importlib.util.find_spec("pyvista"), "PyVista not installed in test environment")
    def test_pyvista_grid_and_vtu_round_trip(self):
        from core.visualization.pyvista_viewer import as_pyvista_grid, export_vtu
        import pyvista as pv
        b=self.beams["quadratic"]
        U=np.zeros(len(b.grid_points)*6)
        U[1::6]=-0.002*b.grid_points[:,2]
        sampled=sample_volume(b,U,subdivisions=2)
        grid=as_pyvista_grid(sampled)
        self.assertEqual(grid.n_cells,b.get_num_elements()*8)
        self.assertEqual(grid.n_points,len(sampled.points))
        with TemporaryDirectory() as d:
            path=Path(d)/"recovered_hermite.vtu"
            export_vtu(sampled,path)
            read=pv.read(str(path))
            self.assertEqual(read.n_cells,grid.n_cells)
            np.testing.assert_allclose(read.point_data["Displacement (in)"],sampled.displacement)
            np.testing.assert_allclose(read.cell_data["Von Mises (psi)"],sampled.cell_fields["Von Mises (psi)"])

    @unittest.skipUnless(
        importlib.util.find_spec("pyvista") and os.getenv("PYVISTA_OFFSCREEN_SMOKE")=="1",
        "Set PYVISTA_OFFSCREEN_SMOKE=1 on a machine with a working graphics driver")
    def test_offscreen_png(self):
        from core.visualization.pyvista_viewer import show_sample
        b=self.beams["linear"]
        U=np.zeros(len(b.grid_points)*6)
        U[1::6]=-0.002*b.grid_points[:,2]
        sampled=sample_volume(b,U,subdivisions=1)
        with TemporaryDirectory() as d:
            path=Path(d)/"screen.png"
            show_sample(sampled,formulation="Linear Hermite",off_screen=True,screenshot=path)
            self.assertGreater(path.stat().st_size,1000)

    def test_bad_input_rejected(self):
        b=self.beams["linear"]
        with self.assertRaises(ValueError):
            sample_volume(b,np.zeros(len(b.grid_points)*6),subdivisions=0)
        with self.assertRaises(ValueError):
            sample_volume(b,np.zeros(1))
        with self.assertRaises(ValueError):
            sample_volume(b,np.zeros(len(b.grid_points)*6),subdivisions=20,max_cells=10)

    @unittest.skipUnless(importlib.util.find_spec("vtk"), "VTK not installed")
    def test_saved_three_demo_vtu_files(self):
        import vtk
        root=Path(__file__).resolve().parents[1]
        for name in ("linear","quadratic","cubic"):
            with self.subTest(name=name):
                path=root/"samples"/f"{name}_3x3x5.vtu"
                self.assertTrue(path.is_file())
                r=vtk.vtkXMLUnstructuredGridReader()
                r.SetFileName(str(path))
                r.Update()
                data=r.GetOutput()
                self.assertEqual(data.GetNumberOfCells(),128)
                self.assertEqual(data.GetNumberOfPoints(),225)
                self.assertIsNotNone(data.GetPointData().GetArray("Displacement (in)"))
                self.assertIsNotNone(data.GetCellData().GetArray("Von Mises (psi)"))

    def test_hexa_connectivity_is_positive_orientation(self):
        b=self.beams["cubic"]
        sample=sample_volume(b,np.zeros(len(b.grid_points)*6),subdivisions=2)
        # VTK Hex ordering: x,y,z directed edges from corner 0 give det>0.
        pts=sample.points[sample.connectivity[0]]
        e=np.stack((pts[1]-pts[0],pts[3]-pts[0],pts[4]-pts[0]),axis=1)
        self.assertGreater(np.linalg.det(e),0)

if __name__=="__main__":
    unittest.main()
