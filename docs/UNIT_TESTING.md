# Unit testing: independent traction integration

Run from the project root in PowerShell with your normal Python; no virtual environment.

**Only the new test (ordinary unittest output):**

```powershell
python -m unittest tests.test_facet_traction_independent.TestIndependentFacetTraction.test_polynomial_vs_gaussian_all_faces_hex27_hex64 -v
```

**Only the new test, printing the detailed numerical results:**

```powershell
python -m tests.test_facet_traction_independent --show-results
```

The `--show-results` flag belongs to the test module, not `python -m unittest`.
It prints one line for each of six faces and each of Hex27 and Hex64: maximum
absolute element-force error in lbf, relative difference, XYZ force resultant
in lbf, and maximum moment difference in lbf-in.

## Independent mathematical method

The production `integrate_traction` evaluates a two-dimensional Gaussian face
quadrature. The reference method instead constructs the cardinal Lagrange
polynomials using `numpy.polynomial.Polynomial`, integrates each analytically
on [-1,1], then multiplies the two tangential integrals and physical area
Jacobian. The reference does not call the production face quadrature helper.
The test compares **every entry** of the consistent nodal force vector,
checks total force against `face area * traction`, and checks the nodal first
moment. It covers all six faces, Hex27 and Hex64, nonsquare elements, offset
origin and a three-component traction vector. Each comparison uses
`rtol=2e-12, atol=2e-12` in consistent lbf / lbf-in units.

Scope: constant traction on axis-aligned affine hexahedra, not general curved
surfaces or spatially varying traction.

**Entire project suite:**

```powershell
python -m unittest discover -s tests -v
```
