# Tutorial001 Revision 4 — compliance audit

| Requirement | Status | Evidence / qualification |
|---|---|---|
| Hermite core source matches certified research source | PASS (byte identity for 5 core files) | `python tools/verify_source_provenance.py`; solver wrapper is changed; Hermite mathematics **not certified** |
| Gravity/self-weight in Hermite FEM | PASS (load-vector assembly tests) | `test_hermite_loads_and_provenance.py`, 0.0174 lbf/in³ acting global -y |
| End-face load integrated | PASS (load-vector assembly tests) | Same tests: 3-point Gauss 2D surface integration, −12.8 lbf resultant |
| Both loads in AWC-style beam reference | PASS (reference tests) | `awc_reference_components()`, signed EB tip deflections; same geometry and specific weight |
| Full-scale Hermite FEM result | NOT RUN IN THIS ENVIRONMENT | User to run `python run_legacy_hermite.py`; no fabricated numerical tip result |
| PyVista Windows view | NOT RUN IN THIS ENVIRONMENT | User to run `python view_legacy_hermite.py --hud` |
| Source-identical Cauchy tutorial | NOT YET | Optional Cauchy tutorial remains independent extracted code |
| Full legacy isolation | PARTIAL | `legacy_hermite/` still contains research files |

**Interpretation:** identical source is not certified Hermite accuracy. The Hermite workflow no longer calibrates the modulus against the AWC beam reference. Its displacement mismatch is reported without adjustment. This does not establish a rigorous certified error bound. The `AWC` reference here is an AWC-style cantilever Euler–Bernoulli formula, not a separately validated AWC design calculation. Gravity is always positive for Hermite. End force and gravity are both in negative global y.
