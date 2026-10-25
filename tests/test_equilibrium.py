import unittest
import numpy as np
from core.hex8 import solve
class EquilibriumTests(unittest.TestCase):
 def test_global_force_and_support_moment(self):
  xyz,cells,u,F,R=solve(nx=2,ny=2,nz=16,return_reactions=True)
  self.assertAlmostEqual(float(F[1::3].sum()),-100,places=8)
  self.assertAlmostEqual(float(R[:,1].sum()),100,places=5)
  loads=F.reshape(-1,3)
  moment=np.cross(xyz,loads+R).sum(axis=0)
  np.testing.assert_allclose(moment,0,atol=1e-5)
if __name__=='__main__':unittest.main()
