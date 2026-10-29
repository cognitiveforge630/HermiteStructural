# Verification architecture

`core.verification.nodal_hex.PolynomialPatchVerifier` checks partition of unity, derivatives, affine displacements/gradients, and reproduction of `y**p` at three selected points on one rectangular cell. This is **not** a multi-element assembled patch test.

`core.verification.nodal_hex.ElementMechanicsVerifier` checks a single rectangular cell's gravity/traction load resultants, applied moments, six rigid-body displacement modes, constant-strain energy, and stiffness symmetry. Neither verifier certifies discretization error.

The original `tools.p_order_study.patch_checks` and `element_verification` call the new classes to preserve existing scripts and JSON output. Run `python -m tools.p_order_study --mode verify --order 3` or `python -m unittest tests.test_nodal_verifier_classes -v`.

The one-cell class currently constructs its reference test cell using the provided element's polynomial order; it deliberately does not reuse the supplied element's geometry or constitutive matrix. This keeps a fixed independent regression fixture, but is not a general verification of arbitrary elements. This affine implementation does not support skew or curved hexahedra.
