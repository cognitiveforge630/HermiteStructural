# Reproducible nodal Hex27 / Hex64 convergence study

Run these commands **from the repository root** in PowerShell. They use the project Python installation directly; no virtual environment is required. This study uses the conventional C0 nodal Lagrange hexahedra, not the recovered-rotation Hermite solver. All physical sizes and loads are inherited unchanged from the existing Redwood benchmark: cross-section 4 × 8 in, length 120 in, downward end load −12.8 lbf, weight density 0.0174 lbf/in³, consistent end-face traction and body-force integration. The clamp fixes all three translations at z=0. **There is no calibration of E.**

The `--levels` values are **original vertex counts** `nx,ny,nz`, not element divisions. A 2,2,4 case has 1×1×3 actual hexahedra. Choose increasing levels deliberately and keep load, material, order, quadrature, and clamp identical when comparing h-refinement.

**SMOKE-TEST WARNING: Small cases (up to 16 source elements) provide unvalidated, rough stress samples only, NOT final stress values or structural design evidence.** See [the non-negotiable stress-sampling policy](NO_CENTROID_STRESS_POLICY.md).

## Order of operations: mesh → solve → report → field maps

### 1. Preview the mesh density (no solve)

```powershell
python -m tools.convergence_study --stage mesh --order 2 --levels "2,2,4;3,3,5;3,4,7"
```

Open `outputs/convergence_hex27/mesh_density.png`. This is a *longitudinal profile of source element edges*, not a 3D interactive mesh. It shows all levels together and gives element/DOF counts before you commit to a run. For Hex64 change `--order 2` to `--order 3`; output folder automatically becomes `convergence_hex64`.

### 2. Solve, with mathematical gates and cost protection

```powershell
python -m tools.convergence_study --stage solve --order 2 --levels "2,2,4;3,3,5;3,4,7"
```

Before solving the first case, the CLI runs the polynomial and one-element mechanics tests. Each level writes `hex27_NxNxN.json`, `hex27_NxNxN.npz` (raw nodal translations and direct inset-subcell-corner stresses), and `hex27_NxNxN_sections.csv`. It refuses **more than 8,000 displacement DOFs** unless you explicitly append `--accept-expensive-run`. Sparse factorization may need far more memory than the raw COO estimate. Start with Hex27 and small cases; Hex64 direct solves grow quickly. For costly runs, decide based on the mesh plan and available workstation RAM rather than uncritically adding that flag.

### 3. Generate an understandable convergence report

```powershell
python -m tools.convergence_study --stage report --order 2 --levels "2,2,4;3,3,5;3,4,7"
```

Read `outputs/convergence_hex27/report.md` in VS Code Markdown Preview (`Ctrl+Shift+V`). `summary.json` contains tabular data for automation. JSON per level includes global equilibrium residuals, reactions, work/strain-energy identity, section shear and moment departures from beam theory, reference tip deflection, and solution timing. Unsolved/skipped levels are marked rather than fabricated. Convergence percentage means difference between consecutive computed meshes, **not** a verified error versus the exact solution.

### 4. Plot displacement, tension, compression, von Mises, shear and analytical maps

```powershell
python -m tools.convergence_study --stage view --order 2 --levels "2,2,4;3,3,5;3,4,7"
```

Open `outputs/convergence_hex27/hex27_2x2x4_fields.png` and its counterparts. For each field, the **left side is FE and the right side is a classical Euler–Bernoulli beam-reference map**; matching rows have common color scales. The heatmaps display an x-near-zero longitudinal slice. For stresses, **one directly evaluated inset corner per visualization subcell** is used in each heatmap pixel; all eight corners are stored for extrema. **Never evaluate stress at centroids.** No nodal stress averaging or false global maximum claims. Tension means `max(sigma_zz,0)` and compression means `max(-sigma_zz,0)`. The final map row compares **absolute** FE and beam-approximate `tau_yz` shear stress (psi), not the section-integrated force. Von Mises uses the six engineering-Voigt Cauchy stress components, with shear order `[xy,yz,xz]`. Displacement displays FE total magnitude against the *absolute theoretical vertical deflection* (a 1D reference); these are not identical quantities. Theoretical stress uses flexure `sigma_zz=M(z)*y/I` and approximate parabolic transverse shear. **This reference is not a 3D elasticity solution** and excludes local support/traction stress effects. Pay particular attention to the fixed end.

### 5. One-command short workflow

```powershell
python -m tools.convergence_study --stage all --order 2 --levels "2,2,4;3,3,5;3,4,7"
```

Separate stages are preferable when a long factorization might be expensive.

### 6. Repeat for cubic Hex64

```powershell
python -m tools.convergence_study --stage mesh --order 3 --levels "2,2,4;2,3,5;3,3,5"
python -m tools.convergence_study --stage solve --order 3 --levels "2,2,4;2,3,5;3,3,5"
python -m tools.convergence_study --stage report --order 3 --levels "2,2,4;2,3,5;3,3,5"
python -m tools.convergence_study --stage view --order 3 --levels "2,2,4;2,3,5;3,3,5"
```

**Caution:** A Hex27 mesh and a Hex64 mesh with the same source-element divisions have *different numbers of global nodes and DOFs*. Comparison of orders at fixed geometry is a p-study, not strictly an h-refinement study.

## Mathematical scope and important limitations

- The reported results are research diagnostics, **not a posteriori certified error bounds** or engineering design approval.
- Positive stress convention: `sigma_zz > 0` is longitudinal tension; `sigma_zz < 0` is longitudinal compression. More general principal-stress tension/compression is not claimed here. Von Mises is a distortion-energy equivalent stress, not a direct wood strength criterion.
- The section resultants in each CSV come from the FE source element's stress at cross-sectional quadrature points. The viewer's inset-corner samples must **not** be substituted for Gaussian section integrals.
- Numerical convergence must be assessed separately for tip displacement, energies, shear and bending resultants, and local stress QoIs. Large peak stresses at a perfectly clamped edge may not converge meaningfully; investigate excluded-edge or percentile QoIs before interpreting them as material strengths.
- The beam theory's modulus is the project's nominal longitudinal Redwood E; no tuning is performed. Orthotropic 3D FE and Euler–Bernoulli 1D theories differ physically, even with perfectly converged meshes.
- No actual high-DOF convergence claim can be made until finer meshes are solved and compared. No Richardson extrapolation, validated effectivity index, or certified error estimator is included.

## Test the workflow

```powershell
python -m unittest discover -s tests -v
python -m unittest tests.test_convergence_study -v
```

`matplotlib` is included in `requirements.txt` for plotting; the numerical solver uses existing NumPy/SciPy dependencies. Viewer PNGs are noninteractive by design and can be opened in VS Code. For full 3D interactive plotting a separate PyVista viewer would require installation and an appropriate graphics environment.
