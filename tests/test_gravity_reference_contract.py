"""Default benchmark: independently check gravity in FEM loads and 1D theory.

Fast analytical/integration regression, NOT a 3D stress-convergence certificate.
"""
import unittest

import numpy as np

from core.elements.nodal_hex import element_type
from core.materials.design_limits import reference_load_screen
from core.reference.awc import DEFAULT_WEIGHT_DENSITY, awc_reference_components
from core.reference.beam_theory import B, E, H, I, L, P
from core.reference.mechanics_of_materials import cantilever_fields, rectangular_stresses


class TestGravityReferenceContract(unittest.TestCase):
    def test_theoretical_stress_and_displacement_include_both_loads(self):
        # Independent section formulas: weight density has units lbf/in^3.
        gamma = DEFAULT_WEIGHT_DENSITY
        area = B * H
        q = gamma * area
        F = -P
        weight = q * L
        moment_end = F * L
        moment_gravity = q * L**2 / 2
        bending = (moment_end + moment_gravity) * (H / 2) / I
        end_deflection = P * L**3 / (3 * E * I)
        weight_deflection = -q * L**4 / (8 * E * I)

        self.assertAlmostEqual(area, 32.0, places=12)
        self.assertAlmostEqual(I, 170.66666666666666, places=12)
        self.assertAlmostEqual(weight, 66.816, places=12)
        self.assertAlmostEqual(moment_end, 1536.0, places=12)
        self.assertAlmostEqual(moment_gravity, 4008.96, places=10)
        self.assertAlmostEqual(bending, 129.96, places=10)
        self.assertAlmostEqual(end_deflection, -0.03927272727272727, places=12)
        self.assertAlmostEqual(weight_deflection, -0.07687636363636363, places=12)

        awc = awc_reference_components(gamma, end_force=P)
        screen = reference_load_screen(tip_force_lbf=P, density_lbf_in3=gamma)
        fields = cantilever_fields(L, length=L, E=E, I=I, P=P, w=-q)
        stress = rectangular_stresses(H / 2, 0., width=B, depth=H,
                                       length=L, P=P, w=-q)
        self.assertAlmostEqual(awc['total_weight_lbf'], weight, places=12)
        self.assertAlmostEqual(awc['weight_per_length_lbf_in'], q, places=12)
        self.assertAlmostEqual(awc['end_load_uy_in'], end_deflection, places=12)
        self.assertAlmostEqual(awc['self_weight_uy_in'], weight_deflection, places=12)
        self.assertAlmostEqual(awc['total_uy_in'], end_deflection + weight_deflection, places=12)
        self.assertAlmostEqual(screen['checks']['bending']['demand'], bending, places=10)
        self.assertAlmostEqual(screen['checks']['example_tip_serviceability']['demand'],
                               -(end_deflection + weight_deflection), places=12)
        self.assertAlmostEqual(float(fields['uy_in']), awc['total_uy_in'], places=12)
        self.assertAlmostEqual(float(stress['sigma_zz_psi']), bending, places=10)
        self.assertAlmostEqual(screen['checks']['parallel_shear']['demand'], 3.732, places=12)

    def test_hex27_hex64_consistent_body_and_end_traction_loads(self):
        # Use a complete beam-shaped element to verify BOTH integrated resultants
        # and moments, separately, without solving a large stiffness system.
        gamma = DEFAULT_WEIGHT_DENSITY
        area = B * H
        weight = gamma * area * L
        for p in (2, 3):
            with self.subTest(order=p):
                elem = element_type(p)((B, H, L), np.eye(6),
                                       origin=(-B / 2, -H / 2, 0.))
                nodes = elem.local_node_coordinates
                gravity = elem.integrate_body_force((0., -gamma, 0.)).reshape(-1, 3)
                elem.clear_force_vector()
                traction = elem.facets['+Z'].integrate_traction(
                    (0., P / area, 0.)).reshape(-1, 3)
                np.testing.assert_allclose(gravity.sum(axis=0),
                                           (0., -weight, 0.), atol=1e-10)
                np.testing.assert_allclose(traction.sum(axis=0),
                                           (0., P, 0.), atol=1e-10)
                gravity_mx = np.cross(nodes, gravity)[:, 0].sum()
                end_mx = np.cross(nodes, traction)[:, 0].sum()
                self.assertAlmostEqual(gravity_mx, weight * L / 2, places=8)
                self.assertAlmostEqual(end_mx, -P * L, places=8)
                self.assertAlmostEqual(gravity_mx + end_mx,
                                       (weight * L / 2) - P * L, places=8)


if __name__ == '__main__':
    unittest.main()
