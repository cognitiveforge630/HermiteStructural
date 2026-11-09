"""Independent exact polynomial checks for Gaussian facet traction integration."""
import argparse
import itertools
import unittest

import numpy as np
from numpy.polynomial import Polynomial

from core.elements.nodal_hex import FACET_NAMES, element_type

SHOW_RESULTS = False


def polynomial_integrals(order):
    nodes = np.linspace(-1.0, 1.0, order + 1)
    result = []
    for i, xi in enumerate(nodes):
        poly = Polynomial([1.0])
        for j, xj in enumerate(nodes):
            if i != j:
                poly *= Polynomial([-xj, 1.0]) / (xi - xj)
        antiderivative = poly.integ()
        result.append(float(antiderivative(1) - antiderivative(-1)))
    return np.asarray(result)


def independent_facet_load(order, extent, face, traction):
    """Integrate tensor-product basis polynomials, without Gaussian quadrature."""
    integrals = polynomial_integrals(order)
    forces = np.zeros(((order + 1)**3, 3))
    axis = 'XYZ'.index(face[1])
    side = order if face[0] == '+' else 0
    tangents = [j for j in range(3) if j != axis]
    jacobian = float(np.prod(np.asarray(extent)[tangents]) / 4)
    for index, triple in enumerate(itertools.product(range(order + 1), repeat=3)):
        if triple[axis] == side:
            forces[index] = (jacobian * integrals[triple[tangents[0]]]
                             * integrals[triple[tangents[1]]] * traction)
    return forces.reshape(-1)


class TestIndependentFacetTraction(unittest.TestCase):
    def test_polynomial_vs_gaussian_all_faces_hex27_hex64(self):
        extent = np.array([2.3, 3.7, 5.1])
        origin = (3.0, -4.0, 7.0)
        traction = np.array([1.25, -0.72, 0.55])
        for order in (2, 3):
            for face in FACET_NAMES:
                with self.subTest(order=order, face=face):
                    element = element_type(order)(extent, np.eye(6), origin=origin)
                    actual = element.integrate_traction(face, traction)
                    expected = independent_facet_load(order, extent, face, traction)
                    delta = float(np.max(np.abs(actual - expected)))
                    relative = delta / max(float(np.max(np.abs(expected))), 1.0)
                    axis = 'XYZ'.index(face[1])
                    area = float(np.prod(np.delete(extent, axis)))
                    force = actual.reshape(-1, 3).sum(axis=0)
                    coords = element.local_node_coordinates
                    moment = np.cross(coords, actual.reshape(-1, 3)).sum(axis=0)
                    expected_moment = np.cross(coords, expected.reshape(-1, 3)).sum(axis=0)
                    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-12)
                    np.testing.assert_allclose(force, area * traction, rtol=2e-12, atol=2e-12)
                    np.testing.assert_allclose(moment, expected_moment, rtol=2e-12, atol=2e-12)
                    if SHOW_RESULTS:
                        print(f'Hex{(order+1)**3} {face}: max |delta Fe|={delta:.6e} lbf; '
                              f'relative={relative:.6e}; '
                              f'F={np.array2string(force, precision=9)} lbf; '
                              f'max |delta M|={np.max(np.abs(moment-expected_moment)):.6e} lbf-in',
                              flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--show-results', action='store_true',
                        help='Print all Hex27/Hex64 face-by-face discrepancies')
    SHOW_RESULTS = parser.parse_args().show_results
    unittest.main(argv=['test_facet_traction_independent', '-v'])
