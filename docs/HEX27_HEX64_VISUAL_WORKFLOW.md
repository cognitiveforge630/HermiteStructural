# Hex27 / Hex64: mesh → solve → inspect results

**MANDATORY: Never evaluate/report stress at centroids.** Stress is sampled
directly at eight 95%-inset corners of each visualization subcell.
**16-cell examples are SMOKE TESTS ONLY; stress is rough, NOT VALIDATED and
NOT FOR STRUCTURAL DESIGN.** See [the project-wide policy](NO_CENTROID_STRESS_POLICY.md).


Run commands from the repository root in PowerShell. `--nx/--ny/--nz` are the **original Hermite vertex counts**; the actual physical element divisions are `(nx-1, ny-1, nz-1)`.

## 0. Install dependencies once

```powershell
python -m pip install -r requirements.txt
```

PyVista is required for interactive 3D viewing and VTU export. Mesh-density PNGs and 2D heatmaps use Matplotlib (included in `requirements.txt`). The 2D FE/beam heatmap now includes an absolute `tau_yz` shear-stress row, distinct from the integrated section-shear resultant.

For interactive wireframe previews for Hex8, Hermite geometry, Hex27 and Hex64, see [MESH_PREVIEW.md](MESH_PREVIEW.md).

## Open the included solved small-mesh samples (no solve needed)

The archive contains **both** conventional Hex27 and Hex64 3×3×5 vertex
mesh (16 physical cells) result JSON, source-section CSV and NPZ field archives.
**WARNING: SMOKE TEST ONLY — INACCURATE COARSE STRESSES; NOT VALIDATED, NOT FINAL, NOT FOR DESIGN.**
These are smoke examples, not the full 192-cell study:

```powershell
python -m tools.nodal_field_viewer samples/hex27_3x3x5_fields.npz --mesh-only
python -m tools.nodal_field_viewer samples/hex27_3x3x5_fields.npz --field "Von Mises (psi)"
python -m tools.nodal_field_viewer samples/hex64_3x3x5_fields.npz --field "Displacement magnitude (in)"
python -m tools.nodal_field_viewer samples/hex64_3x3x5_fields.npz --field "Tension zz (psi)"
python -m tools.nodal_field_viewer samples/hex64_3x3x5_fields.npz --field "Compression zz (psi)"
python -m tools.nodal_field_viewer samples/hex64_3x3x5_fields.npz --field "Von Mises (psi)"
python -m tools.nodal_field_viewer samples/hex64_3x3x5_fields.npz --field "Tau yz (psi)"
```

These scripts require PyVista. The 3D stress coloring uses the most extreme
directly evaluated inset-corner values among eight per visualization subcell.
These are **not certified or complete 3D stress maxima**. Old centroid-only
archives are rejected; rerun the solve to make current archives.

## 1. Hex27 example: 192 physical elements, 6,615 total DOFs

```powershell
# Plan (NO solve, NO visualization result)
python -m tools.p_order_study --mode plan --order 2 --nx 3 --ny 5 --nz 25

# Preview mesh density BEFORE solving (writes mesh_density.png)
python -m tools.convergence_study --stage mesh --order 2 --levels "3,5,25" --out outputs/hex27_3x5x25_preview

# Solve and save JSON, CSV, and an actual field archive (.npz)
python -m tools.p_order_study --mode solve --order 2 --nx 3 --ny 5 --nz 25 --field-output outputs/hex27_3x5x25.npz

# Copy source fields into the preview directory before the stage-view command;
# this workflow's --out is the directory used for its source field archive.
Copy-Item outputs/hex27_3x5x25.npz outputs/hex27_3x5x25_preview/hex27_3x5x25.npz
# FE/beam-reference heatmaps (displacement, longitudinal tension/compression and von Mises)
python -m tools.convergence_study --stage view --order 2 --levels "3,5,25" --out outputs/hex27_3x5x25_preview

# Interactive 3D views from the field archive
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --mesh-only
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --field "Displacement magnitude (in)"
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --field "Uy (in)"
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --field "Tension zz (psi)"
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --field "Compression zz (psi)"
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --field "Von Mises (psi)"
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --field "Shear magnitude (psi)"
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --field "Tau yz (psi)"

# Export all available fields to one ParaView-compatible VTU
python -m tools.nodal_field_viewer outputs/hex27_3x5x25.npz --vtu outputs/hex27_3x5x25.vtu
```

**Heatmap location:** `--stage view` requires the NPZ archive in its `--out` directory, with the filename `hex27_3x5x25.npz`. The preceding `Copy-Item` step puts it there. The integrated workflow below manages all output paths itself.

Alternatively, run the integrated convergence workflow which handles all paths itself:

```powershell
python -m tools.convergence_study --stage mesh --order 2 --levels "3,5,25" --out outputs/hex27_full
python -m tools.convergence_study --stage solve --order 2 --levels "3,5,25" --out outputs/hex27_full
python -m tools.convergence_study --stage report --order 2 --levels "3,5,25" --out outputs/hex27_full
python -m tools.convergence_study --stage view --order 2 --levels "3,5,25" --out outputs/hex27_full
python -m tools.nodal_field_viewer outputs/hex27_full/hex27_3x5x25.npz --field "Von Mises (psi)"
```

## 2. Hex64: same 192 physical elements, 19,929 DOFs

```powershell
python -m tools.p_order_study --mode plan --order 3 --nx 3 --ny 5 --nz 25
python -m tools.convergence_study --stage mesh --order 3 --levels "3,5,25" --out outputs/hex64_full
# CAUTION: large sparse direct factorization; estimate available RAM first.
python -m tools.convergence_study --stage solve --order 3 --levels "3,5,25" --out outputs/hex64_full --accept-expensive-run
python -m tools.convergence_study --stage report --order 3 --levels "3,5,25" --out outputs/hex64_full
python -m tools.convergence_study --stage view --order 3 --levels "3,5,25" --out outputs/hex64_full
python -m tools.nodal_field_viewer outputs/hex64_full/hex64_3x5x25.npz --mesh-only
python -m tools.nodal_field_viewer outputs/hex64_full/hex64_3x5x25.npz --field "Displacement magnitude (in)"
python -m tools.nodal_field_viewer outputs/hex64_full/hex64_3x5x25.npz --field "Tension zz (psi)"
python -m tools.nodal_field_viewer outputs/hex64_full/hex64_3x5x25.npz --field "Compression zz (psi)"
python -m tools.nodal_field_viewer outputs/hex64_full/hex64_3x5x25.npz --field "Von Mises (psi)"
python -m tools.nodal_field_viewer outputs/hex64_full/hex64_3x5x25.npz --field "Shear magnitude (psi)"
python -m tools.nodal_field_viewer outputs/hex64_full/hex64_3x5x25.npz --vtu outputs/hex64_full/hex64_3x5x25.vtu
```

**WARNING: This fast 16-element smoke solve is NOT validated for stress accuracy, and no reported stress values are appropriate for design.**

For a fast 16-element Hex64 smoke solve, substitute `--nx 3 --ny 3 --nz 5` in the `p_order_study` commands, or set `--levels "3,3,5"` in the convergence commands.

## Interpreting the views correctly

- **Displacement** comes from the actual solved displacement DOFs, reported in inches.
- **Tension / compression** mean positive / negative *longitudinal* normal stress `sigma_zz`, respectively. They are not principal-stress extrema.
- **Shear magnitude** is `sqrt(tau_xy²+tau_yz²+tau_xz²)`; use `Tau yz (psi)` for the primary beam shear component.
- **von Mises** is an equivalent stress measure, not a Redwood failure criterion.
- **Stress samples:** Eight independent 95%-inset corner evaluations of the original FE constitutive stress are made for **every visualization subcell**. The viewer uses their sampled extreme per subcell; there is **NO centroid-stress evaluation**, no smoothing across source elements, and no claim of full-field extrema. The PNG heatmap uses one directly sampled inset corner per subcell.
- The 2D paired beam-theory maps are approximations, **not** exact 3D solutions. No certified error bound is supplied.
- `samples/cubic_3x3x9.vtu` belongs to the **separate Hermite sample**, not these Hex27/Hex64 solves. The viewer `core.visualization.pyvista_viewer` is intended for the sampled Hermite export; use `tools.nodal_field_viewer` for these polynomial-order studies.

## Validation

```powershell
python -m unittest tests.test_p_order_study tests.test_nodal_field_workflow -v
```
