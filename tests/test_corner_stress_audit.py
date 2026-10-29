"""Corner and inset-corner stress recovery, with source-side ownership."""
import unittest
import numpy as np
from core.solvers.beam import build_beam
from core.verification.nodal_stress import audit_source_vertex_stress
from core.visualization.sampling import sample_volume

class OriginalCornerStressTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.b=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')

    def test_exact_vertex_B_is_finite_and_nonsingular(self):
        b=self.b
        for e in (0,b.get_num_elements()//2,b.get_num_elements()-1):
            coords=b.get_element_nodes(e)
            for parent in ((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
                           (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)):
                B,det,_=b.primitive_B(*parent,coords)
                self.assertGreater(det,0.)
                self.assertTrue(np.isfinite(B).all())

    def test_exact_corners_recover_affine_stress_without_averaging(self):
        b=self.b
        u=np.zeros((len(b.grid_points),6))
        u[:,2]=0.001*b.grid_points[:,2]
        report=audit_source_vertex_stress(b,u.ravel(),end_force=-12.8,density=.0174)
        self.assertEqual(report['sample_count'],8*b.get_num_elements())
        self.assertEqual(report['sampling_location'],'exact_corners')
        self.assertTrue(report['no_cross_element_nodal_stress_averaging'])
        self.assertEqual(report['max_shared_vertex_element_side_spread_psi'],0.0)
        target=b.elastic_matrix()@np.array([0.,0.,.001,0.,0.,0.])
        for record in report['corner_samples']:
            self.assertIsInstance(record['source_node_index'],int)
            np.testing.assert_allclose(record['stress_psi'],target,atol=1e-8,rtol=1e-8)

    def test_inset_is_inside_source_cell_and_one_sided(self):
        b=self.b
        u=np.zeros(len(b.grid_points)*6)
        report=audit_source_vertex_stress(b,u,location='inset_corners',inset_fraction=.95,
                                          end_force=-12.8,density=.0174)
        self.assertEqual(report['sampling_location'],'inset_corners')
        for item in report['corner_samples']:
            self.assertIsNone(item['source_node_index'])
            self.assertTrue(all(abs(p)<1 for p in item['parent_coordinates']))
            self.assertTrue(np.allclose(item['stress_psi'],0))

    def test_sampling_has_eight_raw_noncenter_stresses_per_display_cell(self):
        b=self.b
        u=np.zeros(len(b.grid_points)*6)
        rng=np.random.default_rng(123)
        u[0::6]=rng.standard_normal(len(b.grid_points))*.001
        sampled=sample_volume(b,u,subdivisions=1)
        self.assertEqual(sampled.stress_sample_location,'exact_corners')
        self.assertEqual(sampled.stress_corner_values.shape,(b.get_num_elements(),8,6))
        np.testing.assert_allclose(sampled.stress_corner_points[0],b.get_element_nodes(0))
        self.assertTrue(np.isfinite(sampled.stress_corner_values).all())
        np.testing.assert_allclose(sampled.cell_fields['Sigma zz (psi)'][0],
                                   sampled.stress_corner_values[0,:,2].mean(),atol=1e-12)
        inset=sample_volume(b,u,subdivisions=1,stress_sample_location='inset_corners',
                            stress_inset_fraction=.9)
        self.assertFalse(np.allclose(inset.stress_corner_points[0],b.get_element_nodes(0)))
        self.assertEqual(inset.stress_inset_fraction,.9)
        with self.assertRaises(ValueError):
            sample_volume(b,u,stress_sample_location='center')
        with self.assertRaises(ValueError):
            sample_volume(b,u,stress_inset_fraction=1)

if __name__=='__main__': unittest.main()
