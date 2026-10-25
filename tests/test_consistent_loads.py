import unittest
import numpy as np
from core.hex8 import solve

class ConsistentLoadTests(unittest.TestCase):
    def test_integrated_traction_and_gravity_resultants(self):
        length, weight, force = 2.0, 0.0174, -10.0
        xyz, cells, u, F, reaction = solve(nx=1, ny=1, nz=2, L=length,
                                           tip_force=force, weight_density=weight,
                                           return_reactions=True)
        expected=force-weight*length
        self.assertAlmostEqual(F[1::3].sum(), expected, places=10)
        self.assertAlmostEqual(reaction[:,1].sum(), -expected, places=7)

    def test_switch_gravity_off(self):
        _,_,_,F = solve(nx=1,ny=1,nz=1,L=1,tip_force=-8,weight_density=0)
        self.assertAlmostEqual(F[1::3].sum(), -8.0)

if __name__ == '__main__':unittest.main()
