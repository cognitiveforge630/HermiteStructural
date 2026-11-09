"""Fast, GUI-independent checks for all supported structured mesh previews."""
import unittest

import numpy as np

from tools.mesh_preview import coordinates, integration_locations, mesh_information


class TestMeshPreview(unittest.TestCase):
    def test_all_families(self):
        expected = {'hex8': 375, 'hermite': 375, 'hex27': 2205, 'hex64': 6643}
        for family, nodes in expected.items():
            with self.subTest(family=family):
                info = mesh_information(family, (3, 5, 25))
                self.assertEqual(info['physical_elements'], 192)
                self.assertEqual(info['global_geometric_nodes'], nodes)
                self.assertEqual(coordinates(info['node_shape']).shape, (nodes, 3))

    def test_gauss_locations_bounded_and_counted(self):
        divisions = (2, 4, 24)
        xyz = integration_locations(3, divisions, 4)
        self.assertEqual(xyz.shape, (192 * 64, 3))
        self.assertTrue(np.all(xyz[:, 0] >= -2))
        self.assertTrue(np.all(xyz[:, 0] <= 2))
        self.assertTrue(np.all(xyz[:, 1] >= -4))
        self.assertTrue(np.all(xyz[:, 1] <= 4))
        self.assertTrue(np.all(xyz[:, 2] >= 0))
        self.assertTrue(np.all(xyz[:, 2] <= 120))

    def test_bad_sizes(self):
        with self.assertRaises(ValueError):
            mesh_information('hex64', (1, 5, 25))


if __name__ == '__main__':
    unittest.main()
