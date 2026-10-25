# Source provenance and remaining certification integration

- Historical `legacy_hermite/solver.py`, `FEMHermiteBeamRegion.py`, and related routines come from `HermiteStructural-Certified_PostBuild.zip`, preserved as historical research source. They are not a mathematically certified general Hermite formulation.
- `core/hex8.py` is an **independent tutorial extraction**, sharing geometry, shape functions, Voigt storage and numerical integration conventions with the research `reference_hex8.py`. **This is not yet a single shared source implementation.**
- The exact-arithmetic a posteriori certificate covers a narrowly defined manufactured Cauchy Hex8 example, not the Redwood cantilever nor legacy Hermite model.
- Integrating both projects behind one canonical Hex8 package and proving that all corresponding routines yield identical matrices remain release gates.
- Historical modulus calibration is explicitly disclosed; matching the calibrated theoretical displacement is not an independent validation.
- Density is specific weight (lbf/in³) and is applied downward along y. The beam reference includes the end force plus uniformly distributed gravitational loading.

## Rev4 audit

Byte-identical to `HermiteStructural-Certified_PostBuild.zip` for `FEMHermiteBeamRegion.py`, `FEMHexahedronRegion.py`, `FEMPointGrid.py`, `reference_hex8.py`, and `checkme.py` (pinned SHA-256 hashes in `HERMITE_SOURCE_SHA256.json`). `legacy_hermite/solver.py` is a modified orchestration layer. This is **source identity for the Hermite core**, not proof of general Hermite mathematical certification. The certified exact-arithmetic example is a Cauchy Hex8 manufactured case. The optional `core/hex8.py` remains independent source and is not unified.

Rev4 runs the orthotropic Hermite element with obligatory positive specific weight, volume-integrated self-weight, and end-surface Gauss-integrated traction. Beam reference uses both forces; the original historical **modulus calibration remains**, so low percentage difference is not independent evidence of agreement.
