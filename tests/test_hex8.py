import unittest
import numpy as np
from core.hex8 import shape,B_matrix,stiffness,elasticity,mesh
class Hex8Tests(unittest.TestCase):
 def test_partition_unity(self):
  N,dN=shape((.2,-.3,.4));self.assertAlmostEqual(N.sum(),1);np.testing.assert_allclose(dN.sum(axis=0),0,atol=1e-14)
 def test_affine_constant_strain(self):
  xyz,cells=mesh(1,1,1);x=xyz[cells[0]]
  for A in [np.diag([1.,0.,0.]),np.diag([0.,1.,0.]),np.diag([0.,0.,1.]),np.array([[0,1,0],[1,0,0],[0,0,0]])]:
   B,_=B_matrix(x,(.1,-.2,.3));e=B@(x@A.T).ravel();expected=np.array([A[0,0],A[1,1],A[2,2],A[0,1]+A[1,0],A[1,2]+A[2,1],A[0,2]+A[2,0]])
   np.testing.assert_allclose(e,expected,atol=1e-13)
 def test_element_symmetric(self):
  x,c=mesh(1,1,1);K=stiffness(x[c[0]],elasticity(1e7,.3));np.testing.assert_allclose(K,K.T,atol=1e-7)
if __name__=='__main__':unittest.main()
