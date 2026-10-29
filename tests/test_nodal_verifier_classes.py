"""Lightweight class API and legacy wrapper regression checks."""
import unittest
import numpy as np
from core.elements.nodal_hex import Hex27Element, Hex64Element
from core.verification.nodal_hex import PolynomialPatchVerifier, ElementMechanicsVerifier
from tools.p_order_study import patch_checks, element_verification

class VerifierClassesTest(unittest.TestCase):
    def test_both_orders(self):
        for cls in (Hex27Element, Hex64Element):
            with self.subTest(element=cls.__name__):
                element = cls((2., 4., 12.), np.eye(6))
                polynomial = PolynomialPatchVerifier().verify(element)
                mechanics = ElementMechanicsVerifier().verify(element)
                self.assertTrue(polynomial['passed'], polynomial)
                self.assertTrue(mechanics['passed'], mechanics)
                self.assertEqual(polynomial, patch_checks(element.order))
                self.assertEqual(mechanics, element_verification(element.order))

if __name__ == '__main__':
    unittest.main()
