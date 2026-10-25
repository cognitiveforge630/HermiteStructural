# Tutorial 001: staging status (2026-10-25)

**This branch has NOT received the complete release tree. Do not merge it.**

The corrected local release `HermiteStructural_Release_Corrected_20261025.zip` is available from the originating ChatGPT conversation, with `Publish-Tutorial001.ps1` for transferring both workspaces from a machine with Git authentication.

Verified in the release: 11 Tutorial 001 unit tests; 6 orthotropic companion tests; five-file SHA-256 historical Hermite source identity test. Primary Hermite uses Gauss-integrated end traction (-12.800 lbf), mandatory body-force self-weight (-66.816 lbf), a combined-load AWC-style Euler–Bernoulli reference, fixed rather than calibrated orthotropic Ez, and matching console/HUD comparison.

Outstanding: upload the full source tree; rerun the production-sized Hermite solve; verify the PyVista screenshot on Windows; then merge reviewed changes into main. Historical Hermite source identity is not a rigorous Hermite a posteriori certificate.
