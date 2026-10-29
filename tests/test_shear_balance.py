"""Independent regression checks of axial shear diagnostic (not certification)."""
import unittest
from unittest.mock import patch
import numpy as np
from numpy.polynomial.legendre import leggauss
from core.solvers.beam import build_beam
from core.solvers.shear_balance import audit_axial_shear_balance, _differentiation_matrix
from core.reference.mechanics_of_materials import cantilever_fields


class AxialShearAuditTests(unittest.TestCase):
    def test_axial_moment_derivative_exact_for_cubic_polynomial(self):
        for order in (3,4,6,8):
            x,_=leggauss(order)
            D=_differentiation_matrix(x)
            for k in range(min(order,4)):
                np.testing.assert_allclose(D@x**k,0*x if k==0 else k*x**(k-1),atol=1e-12)
            np.testing.assert_allclose(D@np.ones(order),np.zeros(order),atol=1e-12)

    def test_exact_analytical_shear_moment_profile_has_no_mismatch(self):
        b=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')
        def exact_source(_b,_U,*,order,axial_positions,end_force,density):
            records=[]
            for iz in range(b.nz-1):
                for p in axial_positions:
                    z=(iz+(p+1)/2)*b.Lz/(b.nz-1)
                    d=cantilever_fields(z,length=b.Lz,E=1,I=1,
                                        P=end_force,w=-density*b.Lx*b.Ly)
                    records.append({'z_in':z,'axial_parent_zeta':p,
                                    'FE_shear_y_lbf':float(d['shear_y_lbf']),
                                    'FE_bending_moment_x_lbf_in':float(d['moment_x_lbf_in']),
                                    'beam_shear_y_lbf':float(d['shear_y_lbf']),
                                    'beam_bending_moment_x_lbf_in':float(d['moment_x_lbf_in'])})
            return {'sections':records}
        with patch('core.solvers.shear_balance.integrate_sections',side_effect=exact_source):
            for o in (4,6):
                result=audit_axial_shear_balance(b,np.zeros(len(b.grid_points)*6),
                    order=o,end_force=-12.8,density=.0174)
                self.assertLess(result['max_interior_sampled_shear_reference_difference_pct'],1e-10)
                self.assertLess(result['max_interior_axially_averaged_shear_reference_difference_pct'],1e-10)
                self.assertLess(result['max_interior_sampled_shear_moment_balance_difference_pct'],1e-9)
                self.assertFalse(result['certified_error_bound'])

    def test_axial_average_cannot_hide_pointwise_shear_error(self):
        b=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')
        def oscillatory(_b,_U,*,order,axial_positions,end_force,density):
            records=[]
            for iz in range(b.nz-1):
                for p in axial_positions:
                    z=(iz+(p+1)/2)*b.Lz/(b.nz-1)
                    exact=cantilever_fields(z,length=b.Lz,E=1,I=1,
                                            P=end_force,w=-density*b.Lx*b.Ly)
                    # Odd oscillation integrates to zero on each cell, but its
                    # maximum local error remains appreciable.
                    feV=float(exact['shear_y_lbf']+8*p)
                    records.append({'z_in':z,'axial_parent_zeta':p,
                        'FE_shear_y_lbf':feV,
                        'FE_bending_moment_x_lbf_in':float(exact['moment_x_lbf_in']),
                        'beam_shear_y_lbf':float(exact['shear_y_lbf']),
                        'beam_bending_moment_x_lbf_in':float(exact['moment_x_lbf_in'])})
            return {'sections':records}
        with patch('core.solvers.shear_balance.integrate_sections',side_effect=oscillatory):
            r=audit_axial_shear_balance(b,np.zeros(len(b.grid_points)*6),
                  end_force=-12.8,density=.0174)
        self.assertLess(r['max_interior_axially_averaged_shear_reference_difference_pct'],1e-10)
        self.assertGreater(r['max_interior_sampled_shear_reference_difference_pct'],5)
        self.assertFalse(r['does_not_modify_solver_or_stress'] is False)

    def test_real_source_zero_stress_and_fail_closed(self):
        b=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')
        zeros=np.zeros(len(b.grid_points)*6)
        result=audit_axial_shear_balance(b,zeros,order=4,end_force=0,density=0)
        self.assertEqual(len(result['cells']),4)
        self.assertEqual(result['max_interior_sampled_shear_reference_difference_pct'],0)
        self.assertEqual(result['max_interior_sampled_shear_moment_balance_difference_pct'],0)
        self.assertFalse(result['certified_error_bound'])
        with self.assertRaises(ValueError):
            audit_axial_shear_balance(b,zeros,order=2,end_force=0,density=0)

if __name__=='__main__': unittest.main()
