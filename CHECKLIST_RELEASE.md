# HermiteStructural Tutorial 001 release checklist

## Mandatory source and physics acceptance

- [x] Historical Hermite element implementation matches five verified files from the previously supplied certified-build research archive (SHA-256 provenance script PASS). **This is source identity, not a rigorously certified Hermite bound.**
- [x] FEM includes specific weight (lbf/in³) in the negative global-y direction using 3D Gauss volume integration. Gravity is mandatory for the primary Hermite demonstration.
- [x] End-face load is an imposed traction integrated over z=L using 2D Gauss surface quadrature, with assembled resultant −12.800 lbf in a small-mesh test.
- [x] Beam-theory/AWC-style reference includes the *same* end-face resultant and gravitational weight (−66.816 lbf); 3EI and 8EI formulas use the appropriate section inertia.
- [x] FEM displacement and combined-load AWC reference are compared in console report and viewer, using the fixed input Ez (no auto-calibration).
- [x] Main workstation tests PASS (11 total); core provenance check PASS.
- [x] Orthotropic matrix companion packaged in `orthotropic_companion/` (6 comparison tests previously reported PASS).
- [ ] Production Hermite solve re-run with corrected reporting on Windows and results reviewed.
- [ ] PyVista HUD and screenshot visually confirmed on Windows.
- [ ] Complete repository source uploaded to GitHub and main branch updated.
- [ ] Rigorous Hermite a posteriori error certificate: **not claimed**; prior manufactured Cauchy certificate does not transfer to Hermite.

## Scientific interpretation

The original Hermite source is historical research code with known affine-mode limitations. It would be misleading to mark a general Hermite certification as complete. The Redwood specific weight of 0.0174 lbf/in³ and Ez of 2,000,000 psi are **demonstration input assumptions**, not universally valid properties of all Redwood lumber grades.
