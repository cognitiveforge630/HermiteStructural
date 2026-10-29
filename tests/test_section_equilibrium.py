"""Source-element section quadrature regressions; no certification claims."""
import unittest
import numpy as np
from core.solvers.beam import build_beam
from core.solvers.section_equilibrium import integrate_sections


class SectionQuadratureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beam=build_beam(nx=3,ny=3,nz=5,formulation='cubic_rotation_hermite')

    def test_rejects_bad_quadrature_and_missing_dofs(self):
        b=self.beam
        zeros=np.zeros(len(b.grid_points)*6)
        with self.assertRaises(ValueError):
            integrate_sections(b,zeros,order=1)
        with self.assertRaises(ValueError):
            integrate_sections(b,zeros,axial_positions=(-1.,))
        with self.assertRaises(ValueError):
            integrate_sections(b,zeros[:-1])

    def test_exact_constant_strain_section_flux(self):
        """Affine extension is an element patch test, not a solved equilibrium test."""
        b=self.beam
        U=np.zeros(len(b.grid_points)*6)
        stretch=1e-4
        U[b.global_dofs[:,2]]=stretch*b.grid_points[:,2]
        data=integrate_sections(b,U,axial_positions=(-.6,0,.6),order=4,
                                end_force=0.,density=0.)
        expected=b.Lx*b.Ly*b.D[2,2]*stretch
        self.assertEqual(len(data['sections']),12)
        for section in data['sections']:
            self.assertAlmostEqual(section['section_area_in2'],b.Lx*b.Ly,places=9)
            self.assertAlmostEqual(section['FE_axial_force_lbf'],expected,places=8)
            self.assertAlmostEqual(section['FE_shear_y_lbf'],0.,places=8)
            self.assertAlmostEqual(section['FE_bending_moment_x_lbf_in'],0.,places=8)
        self.assertFalse(data['certified_error_bound'])
        self.assertFalse(data['visualization_subdivisions_affect_results'])


if __name__=='__main__': unittest.main()
