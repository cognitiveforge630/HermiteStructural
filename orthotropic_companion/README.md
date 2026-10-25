# Tutorial 001 companion — Redwood orthotropic matrices

A standalone VS Code companion lesson preserving the supplied source script **unaltered**. It compares two *algebraic representations of the same compliance law*: a fully written reciprocal-Poisson-pair form (described as eFunda-style in the supplied script), and the symmetric compact form used by the historical HermiteStructural implementation. It does **not** establish independent empirical equivalence of two different material laws or validate the accuracy of the wood property data.

- Material orientation: **R → x, T → y, L → z**.
- Stress ordering: `xx, yy, zz, xy, yz, xz`.
- Strain shear entries are **engineering shear strains** (`gamma`), with compliance `1/G`.
- Elastic constants have units **psi**; compliance entries have **1/psi**.
- Published longitudinal modulus is **calibrated**, not independently measured.
- No PyPardiso and no virtual environment.

Start with [docs/RUN_COMMANDS_IN_ORDER.md](docs/RUN_COMMANDS_IN_ORDER.md). The original script and supplied expected console text are retained at the workspace root.

## Scope and caveats

Matching the two matrices proves a convention/implementation equivalence **for the supplied material parameters and reciprocity assumptions**. It does not prove that a particular FEM implementation assembles the same material matrix, that the Hermite basis is consistent, that the AWC displacement matches, or that the solution is certified.
