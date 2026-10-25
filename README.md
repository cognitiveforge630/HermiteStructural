# HermiteStructural — Tutorial 001

**Start here:** [docs/RUN_COMMANDS_IN_ORDER.md](docs/RUN_COMMANDS_IN_ORDER.md).

This release is **Hermite-first**. The historical orthotropic Redwood cantilever integrates its end-face traction using Gauss quadrature and includes mandatory self-weight. The AWC-style Euler–Bernoulli deflection reference includes both terms. The longitudinal modulus is fixed at the documented reference value rather than calibrated to the predicted displacement. Numerical agreement is an independent displacement comparison, **not** a rigorous Hermite a posteriori certificate.

The original root commands referenced by the existing blog (`solver.py`, `checkme.py`, `check_orthotropic_elasticity.py`, `view.py`, and `view_end_face_traction_integration.py`) are compatibility launchers for `legacy_hermite/`. All newly exposed PyVista text defaults are 10.

- **Primary demonstration:** `python solver.py` then `python view.py --hud`.
- **Orthotropic matrix comparison:** `orthotropic_companion/` has its own VS Code workspace and commands.
- **Optional Cauchy reference:** `python tutorial_001.py` then `python view_cauchy.py --hud`. No need to read its source to complete the Hermite example.
- **Units:** in, lbf, psi. SciPy solver; no PyPardiso or virtual environment needed.

The historical Hermite source has known affine-patch limitations and is preserved for reproducibility. It must not be described as certified. See `docs/SOURCE_PROVENANCE_AND_LIMITATIONS.md`.
