# Two synchronized transverse-shear meshes: classical beam vs cubic Hermite

**Research-grade visualization; not a certified three-dimensional error estimate.**

The cubic cantilever includes two simultaneous PyVista mesh views of signed `tau_yz` (psi):

1. **Left — classical mechanics of materials reference.** This is the usual first-order rectangular-beam `VQ/(It)` shear distribution, with the signed y-direction shear resultant from the applied end force and uniformly distributed self-weight. It is independent of the FE solution; the stress-cube picture represents the tensor components but **does not itself derive the parabolic beam-theory distribution**.
2. **Right — raw cubic Hermite stress.** This is `sigma = D @ (B_Hermite @ (T_element @ U_solved))`, with `tau_yz` in the fifth slot of the symmetric Cauchy stress/engineering shear convention `(xx, yy, zz, xy, yz, xz)`. This is computed from the FE displacement and rotation DOFs, *not* read from the beam-theory function or imposed by a stress fit.

Both views use **exactly the same physical geometry, one common signed blue–white–red scale, matching units and a linked camera**. The geometry is **undeformed**, to prevent deformation magnification from suggesting a different stress location. The wireframe depicts the **original source Hex8 cell geometry** (`subdivisions=1` for this particular paired view), while the source's cubic Hermite stress recovery still supplies the high-order stress. This is not a new eight-node displacement solve. The regular displacement/stress visualization can independently use multiple subcells for smoother geometry.

## Ready-to-open packaged 3×5×25 source mesh

The distribution also includes `samples/cubic_3x5x25_shear_comparison.vtu`, generated from an actual solved 3×5×25 cubic Hermite cantilever. **You can open it immediately without solving again**:

```powershell
python -m core.visualization.shear_comparison `
    samples/cubic_3x5x25_shear_comparison.vtu `
    --full-screen
```

See `samples/cubic_3x5x25_dual_shear_preview.png` for a static **VTK-rendered** illustration from this dataset. The interactive GUI is PyVista; the preview was rendered with the VTK backend in a test environment without PyVista installed.

The same saved-view and solved-view commands are also available in VS Code under
**Terminal → Run Task** as `View: Saved cubic tau_yz comparison (instant; no solve)`
and `Run: Cubic tau_yz side-by-side (3x5x25)`.

## Recommended solved example, on Windows PowerShell

Run from the repository root, with `python -m pip install -r requirements.txt` installed in your global Python interpreter:

```powershell
python -m tutorials.cubic `
    --nx 3 --ny 5 --nz 25 `
    --shear-comparison `
    --full-screen `
    --edges `
    --out outputs/cubic_3x5x25_shear_views.json
```

This opens one **interactive two-panel PyVista window**, and automatically writes:

- `outputs/cubic_3x5x25_shear_views.json` — solved tip, independent checks and reference metadata.
- `outputs/cubic_3x5x25_shear_views_shear_comparison.vtu` — both **one-sided corner-valued** signed shear fields and their difference, accessible in PyVista/ParaView.
- The regular source-section and research stress-audit JSON/CSV outputs from the cubic tutorial.

For an off-screen PNG (e.g. a figure for the README or a research report):

```powershell
python -m tutorials.cubic `
    --nx 3 --ny 5 --nz 25 `
    --shear-comparison `
    --off-screen `
    --edges `
    --out outputs/cubic_3x5x25_shear_views.json
```

This also writes `outputs/cubic_3x5x25_shear_views_shear_comparison.png` (1900 × 850 window settings; actual image pixel density depends on the rendering environment). Or specify an explicit output path with `--shear-screenshot outputs/my_dual_shear.png`.

For a **VTU-only run without launching PyVista**, use `--shear-export-only` instead of `--shear-comparison`:

```powershell
python -m tutorials.cubic `
    --nx 3 --ny 5 --nz 25 `
    --shear-export-only `
    --out outputs/cubic_3x5x25_shear_views.json
```

To reopen a saved dual-shear VTU, without a second expensive solve:

```powershell
python -m core.visualization.shear_comparison `
    outputs/cubic_3x5x25_shear_views_shear_comparison.vtu `
    --full-screen
```

To make an image from a saved file:

```powershell
python -m core.visualization.shear_comparison `
    outputs/cubic_3x5x25_shear_views_shear_comparison.vtu `
    --off-screen `
    --screenshot outputs/cubic_shear_side_by_side.png
```

## Exact mathematical definitions and provenance

Let `z` run from the fixed face (0) to the loaded end (`L`), `y` span the beam depth `h = Ly`, `A = Lx * Ly`, applied signed end load `P <= 0`, and positive weight density `rho` in lbf/in³. Classical rectangular-beam approximation:

```
V_y(z) = P - rho * A * (L - z)
tau_yz^beam(y,z) = (3/2) * V_y(z)/A * [1 - (2*y/h)^2]
```

This is the first-order beam approximation, **not an exact orthotropic three-dimensional elasticity solution**. In particular, free-surface tractions and three-dimensional warping cannot be judged solely from this reference.

For every original FE element `e`, `T_e` maps the global six-DOF solution onto the local Hermite primitive coefficients. At each selected source-element-side corner sample `(xi, eta, zeta)`:

```
strain = B_Hermite(xi,eta,zeta) @ (T_e @ U_e)
stress = D_orthotropic @ strain
tau_yz^FE = stress[4]
```

Only the reference shear uses `analytical_stress_at(...)`. FE shear comes **directly from the solved source-element field**. The difference array is `tau_yz^FE - tau_yz^beam` and is a *research comparison*, not an error bound.

The visualization uses `sample_volume(..., stress_sample_location="exact_corners")` by default. **All eight original FE element-corner stress evaluations per source Hex8 are retained separately**. A shared location may have more than one one-sided FE stress value. The comparison VTU deliberately duplicates the eight original FE corner points for each source cell so VTK **cannot falsely smooth one-sided discontinuities through shared point-data averaging**. The export's `Continuous nodal stress averaging` and `Certified 3D error bound` metadata are both `0`.

`--stress-sample-location inset_corners --stress-inset-fraction 0.95` is available for comparison. In inset mode the scalar samples are inside the visual corners, while the geometry remains a contiguous full-size display mesh. Therefore, **the display corner positions are approximate representations of those inset evaluation positions**; use the exact-corner default for a geographically precise visual comparison.

### Colors and interpretation

- Blue = more negative signed `tau_yz`, white ≈ 0, red = more positive.
- **The same numeric color bounds** cover both reference and FE so the colors can be compared visually.
- The two panes show `tau_yz` only, not von Mises stress, resultant sectional shear, magnitude of a shear vector, or the full stress tensor.
- Good agreement with first-order beam shear does not prove 3D equilibrium, strain-energy accuracy, or an a posteriori certified bound. The separate `--shear-balance`, `--equilibrium-diagnostics`, and `--nodal-stress-audit` remain important.

## Validation, source files, and remaining limitations

- FE sampling: `core/visualization/sampling.py`, computed from `beam.primitive_B` and orthotropic `D`.
- Reference: `core/visualization/cubic_stress.py::analytical_stress_at`.
- Independent two-pane export and PyVista viewer: `core/visualization/shear_comparison.py`.
- Tutorial wiring: `tutorials/cantilever.py` (`python -m tutorials.cubic ...`).
- Unit tests: `python -m unittest tests.test_shear_comparison_view -v`; whole suite: `python -m unittest discover -s tests -v`.

The **3×5×25** mesh is recommended here because earlier research tests found approximately **1.669% discrete interior original-source-corner shear L2 difference** from the beam reference. This is an unsigned sample norm; at 3×5×25 the least-squares fitted FE-to-beam amplitude is about 1.0066, while pointwise values occur both above and below the reference. Axial refinement on the same 3×5 cross-section lowers the discrete difference to 0.918% at 41 nodes and 0.912% at 49 nodes. These are beam-reference sample comparisons, not certified 3D errors; see the [transverse-warping benchmark](TRANSVERSE_CUBIC_SHEAR_RECOVERY_20261027.md) for the full sampling definition and limitations. Raw shear–moment-gradient balance and rigorous stress equilibration are still open issues. This plotting feature **does not alter the stiffness, loads, boundary conditions, B-matrices, constitutive law or FE solution**.
