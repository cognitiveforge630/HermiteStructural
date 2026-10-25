import unittest
import numpy as np
from compare_efunda_vs_hermite_orthotropic import (
 E_LONGITUDINAL, redwood_constants, S_efunda_full, S_hermite
)

class OrthotropicEquivalenceTests(unittest.TestCase):
 def setUp(self):
  self.c = redwood_constants(E_LONGITUDINAL)
  self.S1 = S_efunda_full(self.c)
  self.S2 = S_hermite(self.c)

 def test_compliance_equivalence(self):
  np.testing.assert_allclose(self.S1, self.S2, rtol=0, atol=1e-18)

 def test_stiffness_equivalence(self):
  np.testing.assert_allclose(np.linalg.inv(self.S1), np.linalg.inv(self.S2), rtol=0, atol=1e-8)

 def test_symmetric_positive_definite(self):
  np.testing.assert_allclose(self.S1, self.S1.T, rtol=0, atol=1e-18)
  self.assertGreater(float(np.linalg.eigvalsh(self.S1).min()), 0)

 def test_wood_axis_mapping(self):
  self.assertEqual(self.c['Ez'], E_LONGITUDINAL)
  self.assertAlmostEqual(self.c['Ex'] / self.c['Ez'], .087)
  self.assertAlmostEqual(self.c['Ey'] / self.c['Ez'], .089)

 def test_poisson_reciprocity(self):
  for a,b in [('nuxy','nuyx'),('nuxz','nuzx'),('nuyz','nuzy')]:
   first,second=a[-2:],b[-2:]
   d={'x':self.c['Ex'],'y':self.c['Ey'],'z':self.c['Ez']}
   self.assertAlmostEqual(self.c[a]/d[first[0]],self.c[b]/d[second[0]])

 def test_breaking_reciprocity_is_detected(self):
  bad=dict(self.c)
  bad['nuyx'] *= 1.1
  self.assertGreater(np.max(np.abs(S_efunda_full(bad)-self.S2)),1e-9)

if __name__=='__main__': unittest.main()
