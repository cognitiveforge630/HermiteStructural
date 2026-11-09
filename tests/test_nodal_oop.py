"""Fast OOP mathematical regression and six-face physical loading tests."""
import json
from pathlib import Path
import unittest
import numpy as np
from scipy.sparse import coo_matrix

from core.elements.nodal_hex import (
    FACET_NAMES, Hex27Element, Hex64Element, element_type,
    face_shape_rule, gauss_rule, volume_shape_rule,
)
from core.materials.redwood import redwood_orthotropic_material
from core.mesh.nodal_hex import StructuredNodalHexMesh
from core.reference.beam_theory import E
from core.solvers.nodal_assembly import GlobalAssembler
from tests.reference_legacy_p_order import element_kernel as legacy_element_kernel
from tools.p_order_study import element_kernel, solve_study

D = redwood_orthotropic_material(E).elastic_matrix()


class TestNodalOOP(unittest.TestCase):
    def test_quadrature_rule_reused_and_immutable(self):
        a = gauss_rule(4)
        self.assertIs(a, gauss_rule(4))
        self.assertIs(face_shape_rule(3, 4, '+Z'), face_shape_rule(3, 4, '+Z'))
        self.assertEqual(len(face_shape_rule(3, 4, '+Z')), 16)
        self.assertEqual(len(volume_shape_rule(3, 4)), 64)
        with self.assertRaises(ValueError):
            a[0][0] = 0

    def test_element_and_facet_quadrature_getters(self):
        e = Hex64Element((2., 3., 5.), D)
        x3, w3 = e.get_gaussian_quadrature()
        x2, w2 = e.facets['-X'].get_gaussian_quadrature()
        self.assertIs(x3, e.get_gaussian_quadrature()[0])
        self.assertEqual(x3.shape, (64, 3))
        self.assertEqual(x2.shape, (16, 3))
        np.testing.assert_allclose(x2[:,0], -1)
        self.assertAlmostEqual(float(w3.sum()), 8.)
        self.assertAlmostEqual(float(w2.sum()), 4.)
        self.assertFalse(x3.flags.writeable)
        self.assertFalse(w2.flags.writeable)
        self.assertEqual(e.get_stiffness_matrix().flags.writeable, False)

    def test_corner_node_names_follow_xyz_convention(self):
        for p in (2, 3):
            with self.subTest(order=p):
                node_ids = np.arange((p+1)**3) + 1000
                e = element_type(p)((2, 4, 12), D, node_ids=node_ids)
                n = p+1
                self.assertEqual(e.corner_nodes['lower_back_left'], 1000)
                self.assertEqual(e.corner_nodes['upper_back_right'], 1000 + p*n*n + p*n)
                self.assertEqual(e.corner_nodes['upper_front_right'], 1000+n**3-1)
                self.assertEqual(len(e.corner_nodes), 8)
                self.assertEqual(set(e.facets), set(FACET_NAMES))
                for face in e.facets.values():
                    self.assertEqual(len(face.node_ids), n*n)
                    np.testing.assert_array_equal(face.node_ids,
                                                  node_ids[list(face.local_node_indices)])

    def test_legacy_element_kernel_equivalence(self):
        # Original method has an independent reference copy retained in tests.
        for p in (2, 3):
            with self.subTest(order=p):
                extent = (2., 4., 12.)
                old = legacy_element_kernel(p, extent, D, 0.02, -1.25, p+1)
                new = element_kernel(p, extent, D, 0.02, -1.25, p+1)
                for old_array, new_array in zip(old, new):
                    np.testing.assert_allclose(new_array, old_array, rtol=3e-12, atol=2e-7)

    def test_all_six_face_traction_resultants_and_moments(self):
        # Non-cubic affine geometry and nonsymmetric traction catch axis mixups.
        sizes = np.array([2., 3., 5.])
        start = np.array([3., -4., 7.])
        vector = np.array([1.25, -0.72, 0.55])
        for p in (2, 3):
            for name in FACET_NAMES:
                with self.subTest(order=p, face=name):
                    e = element_type(p)(sizes, D, origin=start)
                    f = e.facets[name]
                    self.assertAlmostEqual(np.linalg.norm(f.normal), 1.0)
                    result = f.integrate_traction(vector).reshape(-1, 3)
                    force = result.sum(axis=0)
                    area = f.area
                    np.testing.assert_allclose(force, area*vector, atol=2e-13, rtol=2e-13)
                    center = start + sizes/2
                    center[f.axis] += sizes[f.axis]/2 * f.normal[f.axis]
                    moment = np.cross(e.local_node_coordinates, result).sum(axis=0)
                    np.testing.assert_allclose(moment, np.cross(center, area*vector),
                                               atol=2e-12, rtol=1e-12)
                    opposite_name = ('-' if name[0]=='+' else '+') + name[1]
                    self.assertTrue(np.allclose(f.normal, -e.facets[opposite_name].normal))
                    self.assertEqual(np.count_nonzero(np.any(abs(result)>1e-12,axis=1)), (p+1)**2)

    def test_body_force_volume_and_moment(self):
        for p in (2, 3):
            e = element_type(p)((2., 3., 5.), D, origin=(3., -4., 7.))
            body = np.array([0., -0.2, 0.])
            nodal_forces = e.integrate_body_force(body).reshape(-1, 3)
            f = np.sum(nodal_forces, axis=0)
            np.testing.assert_allclose(f, 30*body, rtol=1e-13, atol=1e-13)
            center = e.origin + e.extent/2
            moment = np.cross(e.local_node_coordinates, nodal_forces).sum(axis=0)
            np.testing.assert_allclose(moment, np.cross(center, 30*body), rtol=1e-12, atol=1e-12)
            returned = e.get_force_vector()
            returned[:] = 100
            self.assertFalse(np.allclose(e.get_force_vector(), returned))
            e.clear_force_vector()
            np.testing.assert_array_equal(e.get_force_vector(), 0)

    def test_sparse_global_K_and_F_vs_frozen_legacy_COO(self):
        # Two adjacent cells: checks shared nodes and CSR accumulation.
        for p in (2, 3):
            with self.subTest(order=p):
                mesh = StructuredNodalHexMesh(p, (2, 1, 1), (4., 3., 5.))
                extent = mesh.element_extent
                old_K, old_fg, old_ft = legacy_element_kernel(p, extent, D, 0.02, -1.25, p+1)
                triplets_rows=[];triplets_cols=[];triplets_values=[]
                old_force=np.zeros(mesh.global_dof_count)
                element_count=0
                template=element_type(p)(extent,D)
                shared_K=template.get_stiffness_matrix()
                gravity=template.integrate_body_force((0., -0.02, 0.))
                assembler=GlobalAssembler(mesh.global_dof_count, batch_elements=1)
                for _,e in mesh.iter_elements(D,shared_stiffness=shared_K):
                    e.add_force_vector(gravity)
                    e.facets['+Z'].integrate_traction((0.,-1.25,0.))
                    assembler.add_element(e)
                    ids=e.dof_indices
                    triplets_rows.append(np.repeat(ids,len(ids)))
                    triplets_cols.append(np.tile(ids,len(ids)))
                    triplets_values.append(old_K.ravel())
                    old_force[ids] += old_fg+old_ft
                    element_count+=1
                self.assertEqual(element_count, 2)
                K,F=assembler.finalize()
                frozen=coo_matrix((np.concatenate(triplets_values),
                      (np.concatenate(triplets_rows),np.concatenate(triplets_cols))),
                        shape=K.shape).tocsr()
                np.testing.assert_allclose(K.toarray(),frozen.toarray(),rtol=3e-12,atol=3e-7)
                np.testing.assert_allclose(F,old_force,rtol=3e-12,atol=2e-13)
                self.assertLess(np.max(abs((K-K.T).data)), 1e-7)

    def test_small_solve_baseline_qoi(self):
        # Saved JSON was created from user-supplied code before editing.
        path = Path(__file__).resolve().parents[1] / 'docs/verification/legacy_baseline_3x3x5.json'
        baseline = json.loads(path.read_text())
        for p in (2, 3):
            with self.subTest(order=p):
                old = baseline[str(p)]
                new = solve_study(p, 3, 3, 5)
                for name in ('tip_max_abs_uy_in', 'tip_face_mean_uy_in',
                             'max_midspan_shear_y_difference_pct',
                             'max_midspan_moment_x_difference_pct'):
                    self.assertAlmostEqual(new[name], old[name], delta=2e-7)
                self.assertLess(new['equilibrium']['force_defect_inf_lbf'], 1e-6)
                self.assertLess(new['equilibrium']['moment_defect_inf_lbf_in'], 1e-5)
                self.assertFalse(new['certified_error_bound'])

    def test_reject_skewed_geometry_claims_and_wrong_facets(self):
        for p in (2,3):
            cls=element_type(p)
            with self.assertRaises(ValueError):
                cls((2,-1,5),D)
            with self.assertRaises(ValueError):
                cls((2,1,5),D,quadrature_order=p)
            e=cls((2,3,5),D)
            with self.assertRaises(ValueError):
                e.integrate_traction('FRONT',(0,-1,0))
            with self.assertRaises(ValueError):
                e.integrate_traction('+Z',(0,-1))

if __name__ == '__main__':
    unittest.main()
