# PyVista viewer — linear, quadratic, and cubic Hermite

**The same viewer works for all three formulations.** It evaluates the **original recovered-rotation Hermite interpolation** when building the display mesh; it does not replace it with an ordinary first-order Hex8 solution. All quantities use inches, lbf, and psi.

## Installation

From the root of the workspace in Windows PowerShell (using your global Python, **no virtual environment**):

```powershell
python -m pip install -r requirements.txt
python -c "import pyvista as pv; print('PyVista', pv.__version__)"
```

PyVista requires VTK and a graphics-capable operating system for interactive rendering. Exporting VTU does **not** require opening a window. For unsupported remote/headless graphics, export a VTU and open it in ParaView on your workstation.

## Run one viewer per demonstration

The following three commands each solve the FE equations, launch an interactive PyVista window, and separately export a ParaView-compatible `.vtu` file:

```powershell
python -m tutorials.linear    --nx 3 --ny 3 --nz 5 --view --export-vtu
python -m tutorials.quadratic --nx 3 --ny 3 --nz 5 --view --export-vtu
python -m tutorials.cubic     --nx 3 --ny 3 --nz 5 --view --export-vtu
```

**Higher-resolution exploratory quadratic example** (the actual solve is more expensive; the historical comparison used a different modulus):

```powershell
python -m tutorials.quadratic --nx 4 --ny 5 --nz 17 --view --full-screen --export-vtu
```

The three `.vtu` files appear in `outputs/linear.vtu`, `outputs/quadratic.vtu`, `outputs/cubic.vtu`. Each run also produces a JSON text summary. The viewer is launched **after the solve**, with the original solution vector still in memory; results are never approximated from the JSON summary.

To render an image **without opening an interactive window** (if graphics off-screen support is available):

```powershell
python -m tutorials.quadratic --nx 3 --ny 3 --nz 5 --off-screen --screenshot outputs/quadratic.png
```

To save a PNG while also displaying the viewer, add `--view --screenshot` (or `--screenshot outputs/my_image.png`). To export without displaying: `--export-vtu` on its own. You can customize the initial field with `--field "Von Mises (psi)"` or `--field "Uy (in)"`; by default the viewer shows the positive magnitude `|Uy|`.

**Reuse an already exported mesh, without re-solving:**

```powershell
python -m core.visualization.pyvista_viewer outputs/quadratic.vtu --full-screen
python -m core.visualization.pyvista_viewer outputs/cubic.vtu --field "Von Mises (psi)"
```

## Precomputed sample VTUs shipped in this archive

You can inspect the three existing solve outputs without doing any FEM assembly:

```powershell
python -m core.visualization.pyvista_viewer samples/linear_3x3x5.vtu
python -m core.visualization.pyvista_viewer samples/quadratic_3x3x5.vtu
python -m core.visualization.pyvista_viewer samples/cubic_3x3x5.vtu
```

These three VTUs were exported from the original small-mesh solutions, with directly sampled recovered Hermite displacement (point arrays) and stress (cell arrays). Source-node markers and magnitude fields are preserved in the released VTU metadata. Their VTK XML structure and fields were round-trip checked with the VTK reader.

## Interactive keyboard controls

| Key | Action |
|:---|:---|
| `1` | Vertical displacement **magnitude** `|Uy|` (in), blue = least, red = greatest |
| `2` | Total displacement magnitude (in), blue = least, red = greatest |
| `U` | Signed `Uy` (negative downward; blue = negative, red = positive) |
| `A` | Absolute longitudinal stress `|σzz|` (blue = least, red = greatest) |
| `3` | Von Mises stress (psi) |
| `4` | Maximum principal stress (psi) |
| `5` | Minimum principal stress (psi) |
| `6` | Longitudinal stress `σzz` (psi) |
| `7`, `8`, `9` | `σxx`, `σyy`, `Uz` respectively |
| `0` | Classical beam-reference `σzz` (cubic samples only) |
| `B` | Cubic minus analytical longitudinal stress (cubic only) |
| `T` | Longitudinal tension utilization (cubic only) |
| `C` | Midspan sectional slice |
| `M` | Toggle yellow fixed-end source-node cubes and physical end labels |
| `D` | Toggle deformed/undeformed geometry |
| `+`, `-` | Increase/decrease deformation magnification by 25% |
| `E` | Toggle sampled-cell edges |
| `S` | Save screenshot to requested path, or `outputs/viewer_screenshot.png` |
| `R` | Reset camera |

Rotate, pan, and zoom with the standard PyVista mouse controls. The engineering HUD uses 10 pt lettering. The chosen deformation scale is displayed in the HUD and can be overridden with `--deform-scale 10`, for example. The default scale is set to visually expose deflection even in the overly stiff linear example. **Magnification changes the drawing only, never the physical solution.**

## Physical meaning of the fields

- **Point data**: `Displacement (in)` is a true, three-component vector recovered at each visualization point, with `Ux`, `Uy`, `Uz` and magnitude scalar arrays. Display coordinates are reference coordinates; the screen's deformation is applied at rendering only. In ParaView, use *Warp By Vector* on `Displacement (in)`.
- **Cell data**: von Mises, individual Cauchy stress components, and principal stresses are computed directly from the **nonlocal Hermite strain operator** at each visualization subcell's eight exact corners (or user-selected inset corners),
  with the arithmetic eight-point mean used **only for coloring**. Raw corner-side values
  are used in the stress audit. Stress is not silently averaged across element boundaries.
  The stress Voigt order is `xx, yy, zz, xy, yz, xz`. Stress is kept as cell data because stresses can be discontinuous at original element boundaries. Display refinement creates multiple sampling positions, not new FE solutions.
- **Visualization subdivisions**: `--subdivisions 2` (default) produces 8 display hexes per original element. Use `--subdivisions 3` for a smoother rendering at greater CPU/memory cost, or `--subdivisions 1` for a fast view. This parameter does **not** increase FE solution accuracy or establish convergence.
- The three formulations use *different recovered interpolation*, but the **same** mesh, loading, material, boundary conditions, viewer, and unit system for any like-for-like input dimensions.

## Limitations

A visualization is not proof of structural accuracy or an a posteriori bound. Stress is sampled at **eight individual display-cell corners** and the cell color is their arithmetic mean, not a guaranteed stress at any one location or a pointwise maximum. For original source vertices use the separate `--nodal-stress-audit` report. Stress is not integrated through a guaranteed equilibrated reconstruction. The linear formulation remains dramatically too stiff on the small benchmark. The cubic formulation remains experimental. Geometry interpolation within each small displayed Hex8 is linear; using more subdivisions improves the fidelity of the rendered recovered Hermite shape.

## Cubic stress study additions

Current cubic VTU includes **raw Hermite σzz**, independent **beam-reference σzz**, **cubic-minus-beam difference**, and longitudinal tension/compression and shear utilization. Press `6`, `0`, `b`, `t` respectively; `c` cuts a midspan section. Raw FE stress is *not overwritten*. The viewer deforms geometry for clarity; pass `--deform-scale 1` to disable magnification while retaining the true (small) displacement. See [grade design assumptions, benchmark caveats and book citations](CUBIC_STRESS_AND_REDWOOD.md).
