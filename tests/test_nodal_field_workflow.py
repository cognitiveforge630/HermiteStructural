"""Fast regression tests for the field-export CLI plumbing and viewer field catalog."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from tools import nodal_field_viewer, p_order_study


class TestNodalFieldWorkflow(unittest.TestCase):
    def test_export_archive_from_small_solve(self):
        with tempfile.TemporaryDirectory() as directory:
            dest = Path(directory) / 'small.npz'
            result = p_order_study.solve_study(2, 2, 2, 3, field_output_path=dest)
            self.assertTrue(dest.is_file())
            self.assertIn('field_summary', result)
            with np.load(dest) as archive:
                self.assertEqual(archive['fe_stress'].shape, (16, 8, 6))
                self.assertEqual(str(archive['stress_sampling_scheme']), 'inset_subcell_corners_v2')
                self.assertEqual(archive['displacement'].shape[1], 3)

    def test_cli_passes_field_output(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'fields.npz'
            summary = Path(directory) / 'summary.json'
            with patch.object(p_order_study, 'solve_study', wraps=p_order_study.solve_study) as solver:
                p_order_study.main(['--mode','solve','--order','2','--nx','2','--ny','2','--nz','3',
                                    '--out',str(summary),'--field-output',str(out)])
                self.assertEqual(solver.call_args.kwargs['field_output_path'], out)
            self.assertTrue(out.is_file())

    def test_inset_points_never_at_source_or_subcell_centers(self):
        from core.visualization.nodal_convergence import inset_subcell_parent_points
        for order in (2,3):
            with self.subTest(order=order):
                points=inset_subcell_parent_points(order)
                self.assertEqual(points.shape,(order**3,8,3))
                self.assertTrue(np.all(np.abs(points)<1))
                self.assertFalse(np.any(np.all(np.isclose(points,0),axis=-1)))
                for n, sub in enumerate(np.ndindex((order,)*3)):
                    center=-1+(2*np.array(sub)+1)/order
                    self.assertTrue(np.all(np.abs(points[n]-center)>0))

    def test_old_centroid_archives_are_explicitly_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            old=Path(directory)/'old.npz'
            np.savez(old,fe_stress=np.zeros((1,6)))
            with np.load(old) as data:
                with self.assertRaisesRegex(ValueError,'CENTROID-STRESS ARCHIVE'):
                    nodal_field_viewer.require_inset_archive(data)

    def test_direct_stresses_at_saved_sample_points(self):
        from core.elements.nodal_hex import strain_matrix,tensor_shapes
        from core.mesh.nodal_hex import local_nodes
        from core.materials.redwood import redwood_orthotropic_material
        from core.reference.beam_theory import E
        for order in (2,3):
            with self.subTest(order=order), tempfile.TemporaryDirectory() as folder:
                path=Path(folder)/'field.npz'
                from core.mesh.nodal_hex import StructuredNodalHexMesh
                mesh=StructuredNodalHexMesh(order,(1,1,2),(4.,8.,120.),origin=(-2.,-4.,0.))
                # An arbitrary displacement field makes a constant stress check
                # inadequate: strain must be independently evaluated at each point.
                rng=np.random.default_rng(34)
                U=rng.standard_normal(mesh.global_dof_count)*1e-4
                D=redwood_orthotropic_material(E).elastic_matrix()
                from core.visualization.nodal_convergence import export_fields
                export_fields(path,mesh,U,D,.002, -12.8)
                _,dofs=local_nodes(order,0,0,0,mesh.node_shape)
                with np.load(path) as d:
                    self.assertEqual(d['fe_stress'].shape,(2*order**3,8,6))
                    self.assertEqual(d['stress_sample_coordinates'].shape,(2*order**3,8,3))
                    for sub,corner in ((0,0),(0,7),(order**3-1,3)):
                        xyz=d['stress_sample_coordinates'][sub,corner]
                        xi=2*(xyz-mesh.origin)/mesh.element_extent-1
                        _,grad=tensor_shapes(order,xi)
                        expected=D@strain_matrix(grad/(mesh.element_extent/2)[:,None])@U[dofs]
                        np.testing.assert_allclose(d['fe_stress'][sub,corner],expected,rtol=1e-11,atol=1e-9)

    def test_grid_uses_extreme_of_actual_corner_samples_not_centroid(self):
        import sys
        from types import SimpleNamespace
        from unittest.mock import patch
        from core.mesh.nodal_hex import StructuredNodalHexMesh
        from core.visualization.nodal_convergence import export_fields
        from core.materials.redwood import redwood_orthotropic_material
        from core.reference.beam_theory import E

        class FakeGrid:
            def __init__(self, _cells, _types, _points):
                self.point_data={};self.cell_data={};self.field_data={}

        fake=SimpleNamespace(UnstructuredGrid=FakeGrid,
                             CellType=SimpleNamespace(HEXAHEDRON=12))
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'fields.npz'
            mesh=StructuredNodalHexMesh(3,(1,1,2),(4.,8.,120.),origin=(-2.,-4.,0.))
            rng=np.random.default_rng(45)
            displacements=rng.normal(size=mesh.global_dof_count)*1e-4
            D=redwood_orthotropic_material(E).elastic_matrix()
            export_fields(path,mesh,displacements,D,.002,-12.8)
            with patch.dict(sys.modules, {'pyvista':fake}):
                grid=nodal_field_viewer.build_grid(path)
            with np.load(path) as saved:
                raw=saved['fe_stress']; vm=saved['fe_von_mises']
                self.assertEqual(len(grid.cell_data['Von Mises (psi)']),2*3**3)
                mapped=[np.ravel_multi_index((i//3,j//3,k//3,i%3,j%3,k%3),
                                            (1,1,2,3,3,3))
                        for i in range(3) for j in range(3) for k in range(6)]
                np.testing.assert_allclose(grid.cell_data['Von Mises (psi)'],vm[mapped].max(axis=1))
                np.testing.assert_allclose(grid.cell_data['Tension zz (psi)'],
                                           np.maximum(raw[mapped,:,2],0).max(axis=1))
                np.testing.assert_allclose(grid.cell_data['Compression zz (psi)'],
                                           np.maximum(-raw[mapped,:,2],0).max(axis=1))
                self.assertEqual(grid.field_data['Direct inset corner stress samples per subcell'][0],8)

    def test_named_fields(self):
        for name in ('Displacement magnitude (in)', 'Tension zz (psi)',
                     'Compression zz (psi)', 'Von Mises (psi)',
                     'Shear magnitude (psi)', 'Tau yz (psi)'):
            self.assertIn(name, nodal_field_viewer.FIELDS)


if __name__ == '__main__':
    unittest.main()
