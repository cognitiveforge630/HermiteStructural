# Mesh preview for every supported structured beam family

**Run from the repository root in PowerShell.** Mesh previews do **not** assemble, solve, or export stresses. This viewer is for understanding physical cell density and geometric/interpolation points before a solve.

Install graphical dependency if needed (without a virtual environment):

```powershell
python -m pip install pyvista
```

## One command for all mesh families

```powershell
# Eight-node geometry for the existing Hermite formulations
python -m tools.mesh_preview --family hermite --nx 3 --ny 5 --nz 25

# Conventional Hex8 geometry (mesh-only; no Hex8 p-order solver in this CLI)
python -m tools.mesh_preview --family hex8 --nx 3 --ny 5 --nz 25

# Quadratic C0 Hex27
python -m tools.mesh_preview --family hex27 --nx 3 --ny 5 --nz 25

# Cubic C0 Hex64
python -m tools.mesh_preview --family hex64 --nx 3 --ny 5 --nz 25
```

All four previews show the **same 2 × 4 × 24 = 192 physical cells**. The difference is their interpolation lattice (and, for Hermite, separate enriched unknowns not drawn here):

| Family | Geometric nodes shown in full | Nodes per element |
| --- | ---: | ---: |
| Hermite Hex8 geometry | 375 | 8 |
| Hex8 geometry | 375 | 8 |
| Hex27 | 2,205 | 27 |
| Hex64 | 6,643 | 64 |

`--nx`, `--ny`, `--nz` specify *vertex counts*, not physical-element counts. For example `--nz 25` means 24 physical elements along the beam. The red `--nodes` markers are interpolation nodes, **not** Gaussian quadrature points. The orange `--gauss-points` markers are Gaussian positions (4 per axis for cubic by default). Both overlays are opt-in so the display remains clear on dense meshes.

```powershell
# Only cell wireframe (default)
python -m tools.mesh_preview --family hex64 --nx 3 --ny 5 --nz 25

# Show interpolation nodes, but not Gaussian points
python -m tools.mesh_preview --family hex64 --nx 3 --ny 5 --nz 25 --nodes

# Show actual Gaussian quadrature positions, but not nodes
python -m tools.mesh_preview --family hex64 --nx 3 --ny 5 --nz 25 --gauss-points

# Show both, explicitly using four quadrature points per direction
python -m tools.mesh_preview --family hex64 --nx 3 --ny 5 --nz 25 --nodes --gauss-points --quadrature-order 4

# Get counts without PyVista or any graphical window
python -m tools.mesh_preview --family hex64 --nx 3 --ny 5 --nz 25 --summary-only
```

**Compatibility:** The earlier command still works:

```powershell
python view_hex64_density.py --order 3 --nx 3 --ny 5 --nz 25
python view_hex64_density.py --order 2 --nx 3 --ny 5 --nz 25 --nodes
```

For FE solutions and stress/displacement viewers (a separate step), see [HEX27_HEX64_VISUAL_WORKFLOW.md](HEX27_HEX64_VISUAL_WORKFLOW.md). A Hermite mesh preview does not imply the Hex27/Hex64 solver will solve a Hermite enrichment model.

## Fast regression test

```powershell
python -m unittest tests.test_mesh_preview -v
```
