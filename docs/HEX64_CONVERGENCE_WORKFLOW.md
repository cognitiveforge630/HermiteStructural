# Hex64 mesh convergence: reproducible workflow

Run these commands from the repository root in PowerShell. Use the same mesh levels and loading for successive studies. The `nx,ny,nz` arguments are **legacy vertex counts**, so physical element divisions are `(nx-1, ny-1, nz-1)`.

```powershell
# 1. Preview increasingly dense physical elements without solving.
python -m tools.convergence_study --stage mesh --order 3 --levels "2,2,4;3,3,5;3,4,7;4,5,9" --out outputs/convergence_hex64

# 2. Verify cubic interpolation and element mechanics before solving.
python -m tools.p_order_study --mode verify --order 3

# 3. Solve modest Hex64 meshes first. No expensive-run override needed at or below 8,000 DOFs.
python -m tools.convergence_study --stage solve --order 3 --levels "2,2,4;3,3,5;3,4,7" --out outputs/convergence_hex64

# 4. Refresh the human-readable and machine-readable convergence summaries.
python -m tools.convergence_study --stage report --order 3 --levels "2,2,4;3,3,5;3,4,7" --out outputs/convergence_hex64

# 5. Render displacement, longitudinal tension/compression, von Mises, and |tau_yz| shear-stress comparisons.
python -m tools.convergence_study --stage view --order 3 --levels "2,2,4;3,3,5;3,4,7" --out outputs/convergence_hex64
```

The result files are `mesh_density.png`, `report.md`, `summary.json`, and per-mesh `hex64_*_fields.png`, `.npz`, `.json`, and `_sections.csv`. Missing heatmaps mean that mesh was not solved.

## Larger meshes (explicit opt-in)

The solver uses sparse direct elimination, and each Hex64 cell has 192 translation DOFs. Stiffness assembly and factorization can consume substantial memory. A study above 8,000 global DOFs is deliberately skipped unless `--accept-expensive-run` is supplied. For example, start with a *single* larger case before launching a batch:

```powershell
python -m tools.p_order_study --mode plan --order 3 --nx 4 --ny 5 --nz 9
python -m tools.convergence_study --stage solve --order 3 --levels "4,5,9" --out outputs/convergence_hex64 --accept-expensive-run
python -m tools.convergence_study --stage report --order 3 --levels "2,2,4;3,3,5;3,4,7;4,5,9" --out outputs/convergence_hex64
python -m tools.convergence_study --stage view --order 3 --levels "4,5,9" --out outputs/convergence_hex64
```

**Interpretation:** Hex64 is polynomial degree 3 (64 interpolation nodes per element); it does not automatically imply a dense physical mesh. Increase `nx,ny,nz` to perform *h-refinement* with Hex64 held fixed. Track changes between solved mesh levels in `report.md`, and do not treat agreement with the 1D beam reference as either an exact three-dimensional answer or an a posteriori certified error bound. Stresses are raw finite-element stresses; localized peaks near clamped boundary edges can be mesh-sensitive. The plotted beam reference is an approximate comparison, not a validated full 3D reference stress map.
