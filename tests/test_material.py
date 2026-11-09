import unittest
import numpy as np
from core.materials.redwood import redwood_orthotropic_material
from studies.orthotropic.comparison import E_LONGITUDINAL,redwood_constants,S_efunda_full,S_hermite

class MaterialTests(unittest.TestCase):
    def test_orthotropic_matrix(self):
        a=redwood_orthotropic_material(1_100_000.)
        self.assertAlmostEqual(a.Ez,1_100_000.)
        eig=np.linalg.eigvalsh(a.elastic_matrix())
        self.assertTrue(np.all(eig>0))

    def test_companion_compliance_agrees(self):
        constants=redwood_constants(E_LONGITUDINAL)
        np.testing.assert_allclose(S_efunda_full(constants),S_hermite(constants),atol=1e-18,rtol=0.)
