import unittest
import numpy as np
from core.mesh.structured import structured_nodes, structured_hex8_cells

class MeshTests(unittest.TestCase):
    def test_geometry_and_topology(self):
        xyz=structured_nodes(3,3,5,4.,8.,120.)
        cells=structured_hex8_cells(3,3,5)
        self.assertEqual(xyz.shape,(45,3))
        self.assertEqual(cells.shape,(16,8))
        self.assertTrue(np.all(cells>=0))
        self.assertTrue(np.all(cells<len(xyz)))
        self.assertAlmostEqual(float(xyz[:,2].max()),120.)
