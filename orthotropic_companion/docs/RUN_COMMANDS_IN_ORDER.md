# Orthotropic matrix tutorial — commands in order

Run commands in **PowerShell from this workspace folder**. No virtual environment or PyPardiso is necessary.

```powershell
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python .\compare_efunda_vs_hermite_orthotropic.py
```

The last command displays the orthotropic constants, three reciprocal Poisson pairs, the two compliance matrices `S`, their difference, stiffness matrices `D = inverse(S)`, and their difference. The script uses the **published calibrated** longitudinal modulus by default; it does not independently determine that modulus. For the seed material, edit `E_LONGITUDINAL` to `2_000_000.0` psi.

## Connection to the two beam viewers

This is a **constitutive-matrix test**, not an FEA cantilever solve, an AWC deflection validation, or an a posteriori certificate. It supplies the same type of matrix `D` that a Cauchy Hex8 `B.T @ D @ B` formulation uses, and that the historical Hermite element uses with its own strain operator. Before using it in a beam run, verify the solver's actual component order and engineering-shear convention match those documented here.

Keep this workspace as an independent material-model tutorial; add a relative link to it from the Tutorial001 command guide and eventually from the blog. Do not replace the Cauchy or Hermite viewers with the matrix comparison script.
